"""
global_optimization.py

Global fixed-effects optimization for a 3-compartmental propofol PK model.
Parameterizes the Eleveld covariate equations with tunable theta coefficients,
replacing select exponential age-decay functions with Hill (sigmoidal) equations
to allow higher structural variability while maintaining smoothness.

Usage:
  python3 global_optimization.py --eval          # Baseline evaluation only
  python3 global_optimization.py --optimize      # Run optimizer
  python3 global_optimization.py --maxiter 100   # Run optimizer for N iterations
"""

import os
import sys
import time
import argparse
import math
import numpy as np
import pandas as pd
from scipy.optimize import minimize
import multiprocessing as mp
from functools import partial

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from core_solvers import run_algebraic_solver
from models.eleveld import EleveldModel
from models.eleveld_updated import EleveldUpdatedModel, THETA_DEFAULT as THETA_UPDATED_DEFAULT

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH  = os.path.join(SCRIPT_DIR, '..', 'data', 'supplementary_digital_content_1.txt')

# ═══════════════════════════════════════════════════════════════════════════════
# THETA VECTOR DEFINITION
# ═══════════════════════════════════════════════════════════════════════════════
#
# The theta vector encodes ALL tunable coefficients:
#
# --- Original Eleveld structural coefficients ---
# theta[0]  = Vc_base        (default 6.28)     - Central volume base (L)
# theta[1]  = Vc_W50         (default 33.6)     - Weight sigmoid midpoint for Vc
# theta[2]  = Vc_gamma       (default 1.0)      - Weight sigmoid steepness for Vc
# theta[3]  = V2_base        (default 25.5)     - V2 base volume (L)
# theta[4]  = V3_base        (default 273.0)    - V3 base volume (L)
# theta[5]  = CL_male        (default 1.79)     - CL1 base for males (L/min)
# theta[6]  = CL_female      (default 2.10)     - CL1 base for females (L/min)
# theta[7]  = CL_allom_exp   (default 0.75)     - Allometric exponent for CL1
# theta[8]  = CLmat_PMA50    (default 42.3)     - CL maturation PMA midpoint (weeks)
# theta[9]  = CLmat_gamma    (default 9.06)     - CL maturation steepness
# theta[10] = Q2_base        (default 1.75)     - Q2 base clearance (L/min)
# theta[11] = Q2_immat_coef  (default 1.3)      - Q2 immaturity amplification
# theta[12] = Q3_base        (default 1.11)     - Q3 base clearance (L/min)
# theta[13] = Q3mat_x50      (default 68.3)     - Q3 maturation midpoint (weeks PMA)
# theta[14] = Q3mat_gamma    (default 1.0)      - Q3 maturation steepness
#
# --- NEW: Hill sigmoid replacements for age-decay exponentials ---
# theta[15] = V2_age_top     (default 1.7235)   - V2 age Hill: top value
# theta[16] = V2_age_drop    (default 2.4810)   - V2 age Hill: drop magnitude
# theta[17] = V2_age_A50     (default 81.03)    - V2 age Hill: age midpoint
# theta[18] = V2_age_gamma   (default 1.056)    - V2 age Hill: steepness
# theta[19] = V3_age_top     (default 0.999)    - V3 age Hill: top value (opioid)
# theta[20] = V3_age_drop    (default 1.4798)   - V3 age Hill: drop magnitude
# theta[21] = V3_age_A50     (default 96.12)    - V3 age Hill: age midpoint
# theta[22] = V3_age_gamma   (default 1.043)    - V3 age Hill: steepness
#
# --- CL1 age decay (kept as exponential, tunable rate) ---
# theta[23] = CL1_age_rate   (default 0.00286)  - CL1 age decay rate
#
# --- V3 age decay for NON-opioid (Eleveld uses no decay; we keep = 1.0) ---
# (No extra theta needed)
#
# --- Residual error parameter ---
# theta[24] = sigma_prop     (default 0.3)      - Proportional error magnitude
#
# TOTAL: 25 parameters

N_THETA = 25

def get_default_theta():
    """Returns the Eleveld-equivalent default theta vector."""
    theta = np.zeros(N_THETA)
    # Original Eleveld
    theta[0]  = 6.28      # Vc_base
    theta[1]  = 33.6      # Vc_W50
    theta[2]  = 1.0       # Vc_gamma
    theta[3]  = 25.5      # V2_base
    theta[4]  = 273.0     # V3_base
    theta[5]  = 1.79      # CL_male
    theta[6]  = 2.10      # CL_female
    theta[7]  = 0.75      # CL_allom_exp
    theta[8]  = 42.3      # CLmat_PMA50
    theta[9]  = 9.06      # CLmat_gamma
    theta[10] = 1.75      # Q2_base
    theta[11] = 1.3       # Q2_immat_coef
    theta[12] = 1.11      # Q3_base
    theta[13] = 68.3      # Q3mat_x50
    theta[14] = 1.0       # Q3mat_gamma
    # Hill sigmoids for age-decay (fitted from analyze_population_ranges.py)
    theta[15] = 1.7235    # V2_age_top
    theta[16] = 2.4810    # V2_age_drop
    theta[17] = 81.03     # V2_age_A50
    theta[18] = 1.056     # V2_age_gamma
    theta[19] = 0.999     # V3_age_top
    theta[20] = 1.4798    # V3_age_drop
    theta[21] = 96.12     # V3_age_A50
    theta[22] = 1.043     # V3_age_gamma
    # CL1 age decay rate
    theta[23] = 0.00286   # CL1_age_rate
    # Residual error
    theta[24] = 0.3       # sigma_prop
    return theta

# Theta parameter names for display
THETA_NAMES = [
    'Vc_base', 'Vc_W50', 'Vc_gamma',
    'V2_base', 'V3_base',
    'CL_male', 'CL_female', 'CL_allom_exp',
    'CLmat_PMA50', 'CLmat_gamma',
    'Q2_base', 'Q2_immat_coef', 'Q3_base', 'Q3mat_x50', 'Q3mat_gamma',
    'V2_age_top', 'V2_age_drop', 'V2_age_A50', 'V2_age_gamma',
    'V3_age_top', 'V3_age_drop', 'V3_age_A50', 'V3_age_gamma',
    'CL1_age_rate', 'sigma_prop',
]

# ═══════════════════════════════════════════════════════════════════════════════
# PARAMETERIZED ELEVELD MODEL
# ═══════════════════════════════════════════════════════════════════════════════

def _fsigmoid(x, y, z):
    """Hill sigmoid: x^z / (x^z + y^z)"""
    return x**z / (x**z + y**z)

def _hill_falling(x, top, drop, x50, gamma):
    """Falling Hill: top - drop * x^gamma / (x50^gamma + x^gamma)"""
    return top - drop * x**gamma / (x50**gamma + x**gamma)

def _fffm(weight, height, age, m1f2):
    """Fat-Free Mass (Janmahasatian/Al-Sallami) — NOT reparameterized."""
    bmi = weight / (height / 100.0) ** 2
    if m1f2 == 1:  # male
        factor = 0.88 + (1 - 0.88) / (1 + (age / 13.4)**(-12.7))
        return factor * (9270 * weight) / (6680 + 216 * bmi)
    else:  # female
        factor = 1.11 + (1 - 1.11) / (1 + (age / 7.1)**(-1.1))
        return factor * (9270 * weight) / (8780 + 244 * bmi)

# Fixed reference FFM (35yo, 70kg, 170cm, male)
FFMREF = (0.88 + (1 - 0.88) / (1 + (35 / 13.4)**(-12.7))) * ((9270 * 70) / (6680 + 216 * 24.22145))

def compute_params(theta, weight, age, height, m1f2, opioid):
    """
    Compute the 3-compartment PK parameters for a single patient
    using the parameterized Eleveld equations with Hill sigmoid replacements.
    
    Returns dict with keys: vc, k10, k12, k13, k21, k31, ke0
    """
    toweeks = 52.1429
    pma = age * toweeks + 40
    ffm = _fffm(weight, height, age, m1f2)
    
    # ── Vc (Central Volume) ──────────────────────────────────────────────────
    vc = theta[0] * _fsigmoid(weight, theta[1], theta[2]) / _fsigmoid(70, theta[1], theta[2])
    
    # ── V2 (Rapid Peripheral) ────────────────────────────────────────────────
    # Weight scaling: linear (same as Eleveld)
    # Age scaling: Hill falling sigmoid (REPLACES exponential)
    v2_age_factor = _hill_falling(age, theta[15], theta[16], theta[17], theta[18])
    v2 = theta[3] * (weight / 70) * v2_age_factor
    v2ref = theta[3]
    
    # ── V3 (Slow Peripheral) ─────────────────────────────────────────────────
    if opioid:
        v3_age_factor = _hill_falling(age, theta[19], theta[20], theta[21], theta[22])
    else:
        v3_age_factor = 1.0
    v3 = theta[4] * ffm / FFMREF * v3_age_factor
    v3ref = theta[4]
    
    # ── CL1 (Metabolic Clearance) ────────────────────────────────────────────
    mat_factor = _fsigmoid(pma, theta[8], theta[9]) / _fsigmoid(35 * toweeks + 40, theta[8], theta[9])
    
    if m1f2 == 1:  # male
        base_cl = theta[5]
    else:          # female
        base_cl = theta[6]
    
    cl1 = base_cl * (weight / 70)**theta[7] * mat_factor
    if opioid:
        cl1 *= math.exp(-theta[23] * age)
    
    # Safeguard against negative volumes before computing rate constants
    if vc <= 0 or v2 <= 0 or v3 <= 0 or cl1 <= 0:
        return {'vc': -1} # Invalid flag

    # ── CL2 / Q2 (Rapid Inter-compartmental Clearance) ──────────────────────
    q3mat_age = _fsigmoid(age * toweeks + 40, theta[13], theta[14])
    cl2 = theta[10] * (v2 / v2ref)**0.75 * (1 + theta[11] * (1 - q3mat_age))
    
    # ── CL3 / Q3 (Slow Inter-compartmental Clearance) ───────────────────────
    q3mat_ref = _fsigmoid(35 * toweeks + 40, theta[13], theta[14])
    cl3 = theta[12] * (v3 / v3ref)**0.75 * (q3mat_age / q3mat_ref)
        
    k10 = cl1 / vc
    k12 = cl2 / vc
    k13 = cl3 / vc
    k21 = cl2 / v2
    k31 = cl3 / v3
    ke0 = 0.146 * (weight / 70)**(-0.25)  # not optimized (doesn't affect Cp)
    
    return {
        'vc': vc, 'k10': k10, 'k12': k12, 'k13': k13,
        'k21': k21, 'k31': k31, 'ke0': ke0,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════════

def load_patient_data(data_path, max_patients=None):
    """
    Load and pre-process the dataset.
    Returns a list of dicts, one per patient, with:
      - covariates (weight, age, height, m1f2, opioid)
      - infusion_bounds: list of (start, end, rate) tuples
      - t_obs, c_obs: observation times and concentrations
    """
    df = pd.read_csv(data_path, sep=r'\s+')
    
    patient_ids = df['ID'].unique()
    if max_patients is not None:
        patient_ids = patient_ids[:max_patients]
    
    patients = []
    skipped = 0
    
    for pid in patient_ids:
        p = df[df['ID'] == pid].sort_values('TIME')
        first = p.iloc[0]
        
        # Covariates
        weight = first['WGT']
        age    = first['AGE']
        height = first['HGT']
        m1f2   = int(first['M1F2'])
        opioid = (first['TECH'] == 2)
        
        # Dosing events → infusion bounds
        dosing = p[p['EVID'] == 1]
        infusion_bounds = []
        for _, row in dosing.iterrows():
            t_start = row['TIME']
            if row['RATE'] > 0:
                duration = row['AMT'] / row['RATE']
                infusion_bounds.append((t_start, t_start + duration, row['RATE']))
            elif row['AMT'] > 0 and row['RATE'] == 0:
                infusion_bounds.append((t_start, t_start + 1.0, row['AMT'] / 1.0))
        
        # Observations
        obs = p[(p['EVID'] == 0) & (p['DV'] > 0)]
        if len(obs) == 0:
            skipped += 1
            continue
        
        t_obs = obs['TIME'].values.astype(float)
        c_obs = np.exp(obs['DV'].values.astype(float))   # DV is log(Cp)
        max_time = t_obs.max()
        
        # Pre-compute dense time grid for algebraic solver
        n_points = max(int(max_time * 10), 100)
        t_eval = np.linspace(0, max_time + 10, n_points)
        obs_idx = np.searchsorted(t_eval, t_obs)
        obs_idx = np.clip(obs_idx, 0, len(t_eval) - 1)
        
        patients.append({
            'pid': int(pid),
            'weight': weight, 'age': age, 'height': height,
            'm1f2': m1f2, 'opioid': opioid,
            'infusion_bounds': infusion_bounds,
            't_obs': t_obs, 'c_obs': c_obs,
            't_eval': t_eval, 'obs_idx': obs_idx,
        })
    
    if skipped > 0:
        print(f"  Skipped {skipped} patients with no valid observations.")
    
    return patients


def make_infusion_func(infusion_bounds):
    """Create infusion rate function from pre-computed bounds."""
    def infusion_rate(t):
        rate = 0.0
        for start, end, r in infusion_bounds:
            if start <= t < end:
                rate += r
        return rate
    return infusion_rate


# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL COST FUNCTION
# ═══════════════════════════════════════════════════════════════════════════════

def _evaluate_single_patient(pat, theta):
    """Helper for multiprocessing: evaluates cost for a single patient."""
    try:
        params = compute_params(
            theta,
            pat['weight'], pat['age'], pat['height'],
            pat['m1f2'], pat['opioid']
        )
        
        # Sanity: ensure no negative/zero volumes or clearances
        if params['vc'] <= 0 or params.get('k10', 0) <= 0:
            return pat['pid'], 1e6, len(pat['c_obs'])
        
        infusion_func = make_infusion_func(pat['infusion_bounds'])
        cp_pred, _ = run_algebraic_solver(params, infusion_func, pat['t_eval'])
        cp_at_obs = cp_pred[pat['obs_idx']]
        
        # Symmetric Log-Residual Error proxy wrapper
        # This perfectly penalizes over-predictions and under-predictions symmetrically
        c_obs_safe = np.maximum(pat['c_obs'], 1e-4)
        c_pred_safe = np.maximum(cp_at_obs, 1e-4)
        residuals = np.log(c_obs_safe) - np.log(c_pred_safe)
        cost_p = np.sum(residuals**2)
        if np.isnan(cost_p) or np.isinf(cost_p):
            return pat['pid'], 1e6, len(pat['c_obs'])
        return pat['pid'], cost_p, len(pat['c_obs'])
        
    except Exception:
        return pat['pid'], 1e6, len(pat['c_obs'])


def global_cost(theta, patients, return_per_patient=False, pool=None):
    """
    Evaluate the global cost across all patients.
    Supports multiprocessing via Pool.
    """
    if pool is not None:
        eval_func = partial(_evaluate_single_patient, theta=theta)
        results = pool.map(eval_func, patients)
    else:
        results = [_evaluate_single_patient(pat, theta) for pat in patients]
        
    total_cost = 0.0
    total_obs = 0
    per_patient = [] if return_per_patient else None
    
    for pid, cost_p, n_obs in results:
        total_cost += cost_p
        total_obs += n_obs
        if return_per_patient:
            rmse = np.sqrt(cost_p / n_obs) if n_obs > 0 else float('inf')
            per_patient.append({'pid': pid, 'prop_rmse': rmse, 'cost': cost_p, 'n_obs': n_obs})
            
    if return_per_patient:
        return total_cost, total_obs, per_patient
    return total_cost


def global_cost_eleveld(patients):
    """
    Evaluate the same proportional-error cost using the ORIGINAL EleveldModel class.
    This serves as a verification reference.
    """
    total_cost = 0.0
    total_obs = 0
    
    for pat in patients:
        gender = 0 if pat['m1f2'] == 1 else 1
        model = EleveldModel(
            weight=pat['weight'], age=pat['age'], height=pat['height'],
            gender=gender, opioid=pat['opioid']
        )
        params = model.get_parameters()
        
        infusion_func = make_infusion_func(pat['infusion_bounds'])
        cp_pred, _ = run_algebraic_solver(params, infusion_func, pat['t_eval'])
        cp_at_obs = cp_pred[pat['obs_idx']]
        
        c_obs_safe = np.maximum(pat['c_obs'], 1e-4)
        c_pred_safe = np.maximum(cp_at_obs, 1e-4)
        residuals = np.log(c_obs_safe) - np.log(c_pred_safe)
        cost_p = np.sum(residuals**2)
        
        total_cost += cost_p
        total_obs += len(pat['c_obs'])
    
    return total_cost, total_obs


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description='Global PK model optimization.')
    parser.add_argument('--eval', action='store_true', help='Evaluate baseline only')
    parser.add_argument('--optimize', action='store_true', help='Run optimization')
    parser.add_argument('--maxiter', type=int, default=50, help='Max optimizer iterations')
    parser.add_argument('--patients', type=int, default=None, help='Limit number of patients (for testing)')
    parser.add_argument('--load-theta', type=str, default=None,
                        help='Path to a .npy file with a previous optimized theta vector (warm-start)')
    args = parser.parse_args()
    
    print("=" * 70)
    print("GLOBAL PK MODEL OPTIMIZATION")
    print("=" * 70)
    
    # ── Load data ────────────────────────────────────────────────────────────
    print(f"\nLoading patient data from {DATA_PATH}...")
    t0 = time.time()
    patients = load_patient_data(DATA_PATH, max_patients=args.patients)
    t_load = time.time() - t0
    print(f"  Loaded {len(patients)} patients in {t_load:.2f}s")
    
    total_obs = sum(len(p['c_obs']) for p in patients)
    print(f"  Total observations: {total_obs}")
    
    # ── Evaluate with original EleveldModel (reference) ──────────────────────
    print(f"\n--- Reference: Original EleveldModel class ---")
    t0 = time.time()
    eleveld_cost, eleveld_nobs = global_cost_eleveld(patients)
    t_eleveld = time.time() - t0
    eleveld_rmse = np.sqrt(eleveld_cost / eleveld_nobs)
    print(f"  Cost (proportional SSE): {eleveld_cost:.4f}")
    print(f"  Global proportional RMSE: {eleveld_rmse:.6f}")
    print(f"  Eval time: {t_eleveld:.3f}s  ({t_eleveld/len(patients)*1000:.2f} ms/patient)")
    
    # ── Evaluate with parameterized model (should match) ─────────────────────
    theta0 = get_default_theta()

    # ── Warm-start: load a previously saved optimized theta ───────────────────
    if args.load_theta:
        theta_path = os.path.abspath(args.load_theta)
        if not os.path.isfile(theta_path):
            print(f"  ⚠️  --load-theta: file not found: {theta_path}")
        else:
            loaded = np.load(theta_path)
            if len(loaded) == N_THETA:
                theta0 = loaded
                print(f"  ✅ Warm-starting from saved theta: {theta_path}")
            else:
                print(f"  ⚠️  --load-theta: vector length {len(loaded)} ≠ {N_THETA}, ignoring.")

    # Initialize Multiprocessing Pool
    n_cores = mp.cpu_count()
    print(f"  Initializing multiprocessing pool with {n_cores} cores...")
    pool = mp.Pool(processes=n_cores)
    
    print(f"\n--- Parameterized model (current theta) ---")
    t0 = time.time()
    param_cost, param_nobs, per_patient = global_cost(theta0, patients, return_per_patient=True, pool=pool)
    t_param = time.time() - t0
    param_rmse = np.sqrt(param_cost / param_nobs)
    print(f"  Cost (proportional SSE): {param_cost:.4f}")
    print(f"  Global proportional RMSE: {param_rmse:.6f}")
    print(f"  Eval time: {t_param:.3f}s  ({t_param/len(patients)*1000:.2f} ms/patient)")
    
    # ── Comparison ───────────────────────────────────────────────────────────
    cost_diff = abs(param_cost - eleveld_cost) / eleveld_cost * 100
    print(f"\n--- Comparison ---")
    print(f"  Cost difference: {abs(param_cost - eleveld_cost):.4f} ({cost_diff:.4f}%)")
    if cost_diff < 1.0:
        print("  ✅ Parameterized model matches Eleveld within 1%")
    else:
        print(f"  ⚠️  Difference is {cost_diff:.2f}% — expected due to Hill sigmoid replacement of exp()")
        print("      The Hill age functions are close but not identical to the exponentials.")
    
    # ── Per-patient statistics ───────────────────────────────────────────────
    rmses = [p['prop_rmse'] for p in per_patient]
    print(f"\n--- Per-patient proportional RMSE ---")
    print(f"  Min:    {min(rmses):.4f}")
    print(f"  Median: {np.median(rmses):.4f}")
    print(f"  Mean:   {np.mean(rmses):.4f}")
    print(f"  P90:    {np.percentile(rmses, 90):.4f}")
    print(f"  Max:    {max(rmses):.4f}")
    
    # ── Worst 5 patients ─────────────────────────────────────────────────────
    sorted_pp = sorted(per_patient, key=lambda x: x['prop_rmse'], reverse=True)
    print(f"\n  Worst 5 patients:")
    for pp in sorted_pp[:5]:
        print(f"    PID {pp['pid']:4d}: prop_RMSE={pp['prop_rmse']:.4f} ({pp['n_obs']} obs)")
    
    # ── Optimization ─────────────────────────────────────────────────────────
    if args.optimize:
        print(f"\n{'='*70}")
        print(f"STARTING OPTIMIZATION (maxiter={args.maxiter})")
        print(f"{'='*70}")
        print(f"  Theta vector: {N_THETA} parameters")
        print(f"  Cost function: {len(patients)} patients, {total_obs} observations")
        
        # Bounds: physiological constraints
        lower = np.array([
            1.0, 5.0, 0.1,       # Vc_base, Vc_W50, Vc_gamma
            5.0, 50.0,            # V2_base, V3_base
            0.5, 0.5, 0.5,       # CL_male, CL_female, CL_allom_exp
            10.0, 1.0,           # CLmat_PMA50, CLmat_gamma
            0.5, 0.1, 0.2,       # Q2_base, Q2_immat_coef, Q3_base
            10.0, 0.1,           # Q3mat_x50, Q3mat_gamma
            0.5, 0.1, 10.0, 0.1, # V2_age: top, drop, A50, gamma
            0.3, 0.1, 10.0, 0.1, # V3_age: top, drop, A50, gamma
            0.0001,              # CL1_age_rate
            0.01,                # sigma_prop
        ])
        upper = np.array([
            20.0, 100.0, 5.0,    # Vc_base, Vc_W50, Vc_gamma
            100.0, 1000.0,       # V2_base, V3_base
            5.0, 5.0, 1.5,       # CL_male, CL_female, CL_allom_exp
            200.0, 20.0,         # CLmat_PMA50, CLmat_gamma
            5.0, 5.0, 5.0,       # Q2_base, Q2_immat_coef, Q3_base
            200.0, 5.0,          # Q3mat_x50, Q3mat_gamma
            5.0, 10.0, 200.0, 10.0,  # V2_age
            3.0, 5.0, 200.0, 10.0,   # V3_age
            0.05,                # CL1_age_rate
            1.0,                 # sigma_prop
        ])
        bounds = list(zip(lower, upper))
        
        iter_count = [0]
        best_cost = [param_cost]
        checkpoint_path = os.path.join(SCRIPT_DIR, 'optimized_theta_checkpoint.npy')
        
        def callback(xk):
            iter_count[0] += 1
            cost = global_cost(xk, patients, pool=pool)
            if cost < best_cost[0]:
                best_cost[0] = cost
                improvement = (1 - cost / param_cost) * 100
                print(f"  iter {iter_count[0]:4d}: cost={cost:.4f}  improvement={improvement:+.3f}% (saved)")
                np.save(checkpoint_path, xk)
            else:
                print(f"  iter {iter_count[0]:4d}: cost={cost:.4f}")
        
        t0 = time.time()
        try:
            result = minimize(
                global_cost,
                theta0,
                args=(patients, False, pool),
                method='L-BFGS-B',
                bounds=bounds,
                options={'maxiter': args.maxiter, 'disp': True, 'ftol': 1e-8},
                callback=callback,
            )
            opt_x = result.x
            success_str = str(result.success)
            msg_str = result.message
            nit = result.nit
            nfev = result.nfev
        except KeyboardInterrupt:
            print("\n  [!] Optimization interrupted by user (KeyboardInterrupt).")
            print(f"  [!] Recovering last checkpoint from: {checkpoint_path}")
            if os.path.exists(checkpoint_path):
                opt_x = np.load(checkpoint_path)
            else:
                opt_x = theta0
            success_str = "Interrupted"
            msg_str = "User aborted. Best parameters restored from checkpoint."
            nit = iter_count[0]
            nfev = "Unknown"
        t_opt = time.time() - t0
        
        print(f"\n--- Optimization Result ---")
        print(f"  Success: {success_str}")
        print(f"  Message: {msg_str}")
        print(f"  Iterations: {nit}")
        print(f"  Function evals: {nfev}")
        if isinstance(nfev, int) and nfev > 0:
            print(f"  Time: {t_opt:.1f}s ({t_opt/nfev:.2f}s per function eval)")
        else:
            print(f"  Time: {t_opt:.1f}s")
        
        opt_cost = global_cost(opt_x, patients, pool=pool)
        opt_pseudo_rmse = np.sqrt(opt_cost / total_obs)
        improvement = (1 - opt_cost / eleveld_cost) * 100
        
        print(f"\n  Eleveld Log-Residual cost:   {eleveld_cost:.4f}  (RMSE={eleveld_rmse:.6f})")
        print(f"  Optimized Log-Residual cost: {opt_cost:.4f}  (RMSE={opt_pseudo_rmse:.6f})")
        print(f"  Improvement:                 {improvement:+.3f}%")
        
        print(f"\n--- Theta comparison (default → optimized) ---")
        for i in range(N_THETA):
            changed = '  ←' if abs(opt_x[i] - theta0[i]) / (abs(theta0[i]) + 1e-10) > 0.01 else ''
            print(f"  θ[{i:2d}] {THETA_NAMES[i]:16s}: {theta0[i]:10.4f} → {opt_x[i]:10.4f}{changed}")
            
        theta_out = os.path.join(SCRIPT_DIR, 'optimized_theta.npy')
        print(f"\n  Saving finalized parameters to: {theta_out}")
        np.save(theta_out, opt_x)
        print(f"  To warm-start next run use: --load-theta {theta_out}")
    
    pool.close()
    pool.join()
    
    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")


if __name__ == '__main__':
    main()
