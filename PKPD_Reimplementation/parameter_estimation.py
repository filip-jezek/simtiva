import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, Bounds

# Ensure we can import from local modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from core_solvers import run_algebraic_solver
from models.eleveld import EleveldModel

def get_infusion_func(dosing_rows):
    infusion_bounds = []
    for _, row in dosing_rows.iterrows():
        t_start = row['TIME']
        if row['RATE'] > 0:
            duration = row['AMT'] / row['RATE']
            infusion_bounds.append((t_start, t_start + duration, row['RATE']))
        elif row['AMT'] > 0 and row['RATE'] == 0:
            infusion_bounds.append((t_start, t_start + 1.0, row['AMT'] / 1.0))
            
    def infusion_rate(t):
        rate = 0.0
        for start, end, r in infusion_bounds:
            if start <= t < end:
                rate += r
        return rate
    return infusion_rate

class PKParameterEstimator:
    def __init__(self, t_obs, c_obs, infusion_func, max_time):
        """
        Estimates the parameters of a 3-compartmental model.
        t_obs: array of observation times (minutes)
        c_obs: array of observed plasma concentrations (mg/L)
        infusion_func: function giving infusion rate at time t
        max_time: maximum simulation time needed
        """
        self.t_obs = np.array(t_obs)
        self.c_obs = np.array(c_obs)
        self.infusion_func = infusion_func
        self.max_time = max_time
        
        # Dense time vector for simulation
        self.t_eval = np.linspace(0, self.max_time + 10, int((self.max_time + 10) * 10))
        
        # Cache observation indices to avoid searchsorted on every iteration
        self.obs_idx = np.searchsorted(self.t_eval, self.t_obs)

    def _params_to_k(self, p):
        # p is [V1, V2, V3, CL1, CL2, CL3, KE0]
        # Or just [V1, V2, V3, CL1, Q2, Q3] if we ignore KE0 for Cp optimization
        v1, v2, v3, cl1, q2, q3 = p[:6]
        
        k10 = cl1 / v1
        k12 = q2 / v1
        k21 = q2 / v2
        k13 = q3 / v1
        k31 = q3 / v3
        
        # Default Eleveld ke0 roughly 0.146 for an average person, but it doesn't affect Cp
        ke0 = 0.146 if len(p) == 6 else p[6]

        return {
            'vc': v1,
            'k10': k10,
            'k12': k12,
            'k13': k13,
            'k21': k21,
            'k31': k31,
            'ke0': ke0
        }

    def simulate(self, param_dict):
        cp, _ = run_algebraic_solver(param_dict, self.infusion_func, self.t_eval)
        return self.t_eval, cp

    def objective_least_squares(self, p):
        """Standard least squares objective (Log error)"""
        params = self._params_to_k(p)
        cp_pred, _ = run_algebraic_solver(params, self.infusion_func, self.t_eval)
        
        # Extract predictions at observation times
        cp_at_obs = cp_pred[self.obs_idx]
        
        # We use Log prediction error to give equal weight to high/low concentrations
        # Add small epsilon to avoid log(0)
        c_pred_log = np.log(np.maximum(cp_at_obs, 1e-6))
        c_obs_log = np.log(np.maximum(self.c_obs, 1e-6))
        
        residuals = c_pred_log - c_obs_log
        sse = np.sum(residuals**2)
        return sse

    def objective_map(self, p, prior_means, prior_vars):
        """Maximum A Posteriori objective (adding regularization to population parameters)"""
        sse = self.objective_least_squares(p)
        
        # Add penalty for deviating from prior (Empirical Bayes / MAP)
        penalty = 0
        for i in range(len(p)):
            # Assuming log-normal distribution for PK parameters
            penalty += ((np.log(p[i]) - np.log(prior_means[i]))**2) / (2 * prior_vars[i])
        
        return sse + penalty

    def fit_minimze(self, initial_guess=None, bounds=None, method='local'):
        """
        Fit using scipy.optimize.minimize
        Default bounds prevent unphysiological negative volumes or clearances
        method: 'local' (L-BFGS-B) or 'map' (Requires prior population initial guess)
        """
        if initial_guess is None:
            # Baseline generic adult starting point
            initial_guess = [10.0, 30.0, 200.0, 2.0, 2.0, 1.0]
            
        if bounds is None:
            # Min and Max bounds for [V1, V2, V3, CL1, Q2, Q3]
            lower = [0.1, 1.0, 10.0, 0.1, 0.1, 0.1]
            upper = [100.0, 200.0, 1000.0, 10.0, 15.0, 5.0]
            bounds = Bounds(lower, upper)
            
        print(f"Starting Scipy Optimization ({method})...")
        res = minimize(
            self.objective_least_squares, 
            initial_guess, 
            method='L-BFGS-B', 
            bounds=bounds,
            options={'disp': True, 'maxiter': 500}
        )
        return res

def main():
    parser = argparse.ArgumentParser(description='Estimate PK parameters for a patient.')
    parser.add_argument('--pid', type=int, default=1, help='Patient ID to estimate')
    parser.add_argument('--plot', action='store_true', default=True, help='Plot fit and save')
    args = parser.parse_args()

    # Load data
    dataset_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'supplementary_digital_content_1.txt')
    if not os.path.exists(dataset_path):
        print(f"Dataset not found at {dataset_path}")
        return
        
    df = pd.read_csv(dataset_path, sep=r'\s+')
    p1 = df[df['ID'] == args.pid].sort_values('TIME').copy()
    
    if len(p1) == 0:
        print(f"No data for patient ID {args.pid}")
        return

    # Patient Covariates
    first_row = p1.iloc[0]
    weight = first_row['WGT']
    age = first_row['AGE']
    height = first_row['HGT']
    gender = 0 if first_row['M1F2'] == 1 else 1 
    opioid = first_row['TECH'] == 2

    dosing_rows = p1[p1['EVID'] == 1]
    obs_rows = p1[(p1['EVID'] == 0) & (p1['DV'] > 0)]
    
    if len(obs_rows) < 6:
        print(f"Not enough observations ({len(obs_rows)}) to fit a 6-parameter model fully unconstrained reliably.")
        print("Will attempt anyway, but expect potential overfitting.")

    infusion_func = get_infusion_func(dosing_rows)
    t_obs = obs_rows['TIME'].values
    c_obs = np.exp(obs_rows['DV'].values)  # DV is log(Cp)
    max_time = t_obs.max() if len(t_obs) > 0 else dosing_rows['TIME'].max() + 60

    print(f"--- Processing Patient {args.pid} ---")
    print(f"Covariates: {weight}kg, {age:.0f}yo, {'Female' if gender==1 else 'Male'}")
    print(f"Observations: {len(t_obs)}")
    
    # 1. Baseline Population Model (Eleveld)
    pop_model = EleveldModel(weight=weight, age=age, height=height, gender=gender, opioid=opioid)
    pop_params = pop_model.get_parameters()
    
    # Convert Population k to raw parameters for initial guess
    v1_pop = pop_params['vc']
    cl1_pop = pop_params['k10'] * v1_pop
    q2_pop = pop_params['k12'] * v1_pop
    v2_pop = q2_pop / pop_params['k21']
    q3_pop = pop_params['k13'] * v1_pop
    v3_pop = q3_pop / pop_params['k31']
    
    initial_guess = [v1_pop, v2_pop, v3_pop, cl1_pop, q2_pop, q3_pop]
    
    estimator = PKParameterEstimator(t_obs, c_obs, infusion_func, max_time)
    
    # Calculate baseline RMSE
    base_sse = estimator.objective_least_squares(initial_guess)
    base_rmse = np.sqrt(base_sse / len(t_obs))
    print(f"\n[Population Baseline] Eleveld Log(SSE): {base_sse:.4f}, RMSE on log scale: {base_rmse:.4f}")
    
    # 2. Estimate Best Fit Individual Parameters
    res = estimator.fit_minimze(initial_guess=initial_guess)
    
    if not res.success:
        print("Optimization failed to converge.")
        print(res.message)
    
    opt_p = res.x
    opt_sse = res.fun
    opt_rmse = np.sqrt(opt_sse / len(t_obs))
    print(f"\n[Optimized Individual] Log(SSE): {opt_sse:.4f}, RMSE on log scale: {opt_rmse:.4f}")
    
    print("\nParameter Comparison (Pop vs Opt):")
    param_names = ['V1 (L)', 'V2 (L)', 'V3 (L)', 'CL1 (L/min)', 'Q2 (L/min)', 'Q3 (L/min)']
    for i, name in enumerate(param_names):
        print(f"  {name:12s} : {initial_guess[i]:8.3f} --> {opt_p[i]:8.3f}")

    if args.plot:
        t_eval = estimator.t_eval
        
        # Simulate both
        _, cp_pop = run_algebraic_solver(pop_params, infusion_func, t_eval)
        opt_params_dict = estimator._params_to_k(opt_p)
        _, cp_opt = run_algebraic_solver(opt_params_dict, infusion_func, t_eval)
        
        plt.figure(figsize=(10, 6))
        plt.plot(t_eval, cp_pop, 'b--', label='Eleveld Pop', alpha=0.7)
        plt.plot(t_eval, cp_opt, 'r-', label='Estimated Indiv', linewidth=2)
        plt.scatter(t_obs, c_obs, c='k', marker='X', s=80, label='Observed (Cp)')
        
        plt.xlabel('Time (min)')
        plt.ylabel('Concentration (mg/L)')
        plt.title(f'Patient {args.pid} Parameter Estimation Fit')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.yscale('log')
        plt.tight_layout()
        
        out_name = f'estimation_pid_{args.pid}.png'
        plt.savefig(out_name, dpi=150)
        print(f"\nSaved plot to {out_name}")

if __name__ == '__main__':
    main()
