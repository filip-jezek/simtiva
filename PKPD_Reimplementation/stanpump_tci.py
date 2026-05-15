"""
STANPUMP-based TCI (Target-Controlled Infusion) engine for propofol.

Implements the Effect-Site Concentration Targeting (CET) algorithm from STANPUMP,
as used in the SimTIVA JavaScript application (pharmacology.js).

Algorithm references:
  - Shafer SL, Gregg KM. Algorithms to rapidly achieve and maintain stable drug
    concentrations at the site of drug effect with a computer-controlled infusion pump.
    J Pharmacokinet Biopharm 1992;20(2):147-69.
  - SimTIVA pharmacology.js: cube(), calculate_udfs(), virtual_model(), find_peak(),
    deliver_cet_real(), deliver_cpt()

Units convention:
  - Rate constants (k10, k12 …): /min in pk_params, converted to /s internally
  - Volumes (vc): litres
  - Concentrations (Cp, Ce): mcg/mL  (1 mg/L == 1 mcg/mL)
  - Drug amounts: mg
  - Infusion rate: mg/s internally; converted to ml/h for output
  - Time: seconds
"""

import math
from typing import Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# 1.  Characteristic-polynomial solver  (cube.c translated from STANPUMP)
# ---------------------------------------------------------------------------

def cube(k10: float, k12: float, k21: float, k13: float, k31: float) -> List[float]:
    """
    Solve X³ + a2·X² + a1·X + a0 = 0 for a 3-compartment model.

    All rate constants in /s.  Returns the three real roots sorted ascending
    (most-negative first), which become λ1 ≤ λ2 ≤ λ3 (all negative for a
    physiologically valid model).
    """
    toradian = math.asin(1.0) * 2.0 / 180.0

    a0 = k10 * k21 * k31
    a1 = k10 * k31 + k21 * k31 + k21 * k13 + k10 * k21 + k31 * k12
    a2 = k10 + k12 + k13 + k21 + k31

    # Depressed cubic  y³ + p·y + q = 0
    p = a1 - (a2 * a2 / 3.0)
    q = (2.0 * a2 ** 3 / 27.0) - (a1 * a2 / 3.0) + a0

    r1 = math.sqrt(-(p ** 3) / 27.0)
    phi = (-q / 2.0) / r1
    phi = max(-1.0, min(1.0, phi))          # guard against numerical rounding
    phi = math.acos(phi) / 3.0
    r1 = 2.0 * math.exp(math.log(r1) / 3.0)

    roots = [
        -(math.cos(phi)                      * r1 - a2 / 3.0),
        -(math.cos(phi + 120.0 * toradian)   * r1 - a2 / 3.0),
        -(math.cos(phi + 240.0 * toradian)   * r1 - a2 / 3.0),
    ]
    roots.sort()
    return roots   # [λ1, λ2, λ3] ascending (all negative)


# ---------------------------------------------------------------------------
# 2.  Precompute TCI coefficients  (calculate_udfs translated from STANPUMP)
# ---------------------------------------------------------------------------

def calculate_udfs(pk_params: Dict[str, float], delta_seconds: int = 1) -> dict:
    """
    Precompute unit-disposition functions and exponential coefficients.

    Parameters
    ----------
    pk_params : dict
        {vc, k10, k12, k13, k21, k31, ke0} — rate constants in /min, vc in L.
    delta_seconds : int
        Duration of a bolus "unit infusion" used to define e_udf (default 1 s,
        matching SimTIVA's delta_seconds = 1).

    Returns
    -------
    dict
        lam1..lam4   : characteristic roots in /s (lam4 = ke0/60)
        l1..l4       : per-second decay factors  exp(-lam_i)
        p_coef       : plasma macro-constants  [pc1, pc2, pc3]
        e_coef       : effect-site macro-constants  [ec1, ec2, ec3, ec4]
        p_udf        : list[301] — plasma conc from 1 mg/s at each second 0..300
        e_udf        : dict[int → float] — Ce from delta_seconds of 1 mg/s,
                       then off, at each second up to peak_time+200
        peak_time    : seconds to peak Ce after a delta_seconds-long unit bolus
    """
    # /min → /s
    k10 = pk_params['k10'] / 60.0
    k12 = pk_params['k12'] / 60.0
    k13 = pk_params['k13'] / 60.0
    k21 = pk_params['k21'] / 60.0
    k31 = pk_params['k31'] / 60.0
    k41 = pk_params['ke0'] / 60.0
    vc  = pk_params['vc']

    # --- Characteristic roots ---
    lam1, lam2, lam3 = cube(k10, k12, k21, k13, k31)
    lam4 = k41

    # --- Plasma macro-constants (STANPUMP p_coef) ---
    pc1 = ((k21 - lam1) * (k31 - lam1)
           / ((lam1 - lam2) * (lam1 - lam3) * vc * lam1))
    pc2 = ((k21 - lam2) * (k31 - lam2)
           / ((lam2 - lam1) * (lam2 - lam3) * vc * lam2))
    pc3 = ((k21 - lam3) * (k31 - lam3)
           / ((lam3 - lam2) * (lam3 - lam1) * vc * lam3))
    p_coef = [pc1, pc2, pc3]

    # --- Effect-site macro-constants (STANPUMP e_coef) ---
    ec1 = pc1 / (k41 - lam1) * k41
    ec2 = pc2 / (k41 - lam2) * k41
    ec3 = pc3 / (k41 - lam3) * k41
    ec4 = ((k41 - k21) * (k41 - k31)
           / ((lam1 - k41) * (lam2 - k41) * (lam3 - k41) * vc))
    e_coef = [ec1, ec2, ec3, ec4]

    # --- Per-second decay factors ---
    l1 = math.exp(-lam1)
    l2 = math.exp(-lam2)
    l3 = math.exp(-lam3)
    l4 = math.exp(-lam4)

    # --- p_udf: plasma conc from 1 mg/s continuous infusion, each second ---
    p_udf = [0.0] * 301
    t1 = t2 = t3 = 0.0
    for i in range(1, 301):
        t1 = t1 * l1 + pc1 * (1.0 - l1)
        t2 = t2 * l2 + pc2 * (1.0 - l2)
        t3 = t3 * l3 + pc3 * (1.0 - l3)
        p_udf[i] = t1 + t2 + t3

    # --- e_udf: Ce from delta_seconds of 1 mg/s, then off ---
    # Phase 1 — infusion on for delta_seconds
    e_udf: Dict[int, float] = {0: 0.0}
    t1 = t2 = t3 = t4 = 0.0
    for i in range(1, delta_seconds + 1):
        t1 = t1 * l1 + ec1 * (1.0 - l1)
        t2 = t2 * l2 + ec2 * (1.0 - l2)
        t3 = t3 * l3 + ec3 * (1.0 - l3)
        t4 = t4 * l4 + ec4 * (1.0 - l4)
        e_udf[i] = t1 + t2 + t3 + t4

    # Phase 2 — infusion off; track until peak then 200 extra seconds
    prior = e_udf[delta_seconds]
    peak_time = delta_seconds
    found_peak = False
    for i in range(delta_seconds + 1, 5000):
        t1 *= l1
        t2 *= l2
        t3 *= l3
        t4 *= l4
        val = t1 + t2 + t3 + t4
        e_udf[i] = val
        if not found_peak:
            if prior >= val:
                peak_time = i - 1
                found_peak = True
            prior = val
        elif i > peak_time + 200:
            break

    return {
        'lam1': lam1, 'lam2': lam2, 'lam3': lam3, 'lam4': lam4,
        'l1': l1, 'l2': l2, 'l3': l3, 'l4': l4,
        'p_coef': p_coef,
        'e_coef': e_coef,
        'p_udf': p_udf,
        'e_udf': e_udf,
        'peak_time': peak_time,
        'delta_seconds': delta_seconds,
        'vc': vc,
        'k10_s': k10, 'k12_s': k12, 'k13_s': k13,
        'k21_s': k21, 'k31_s': k31, 'k41_s': k41,
    }


# ---------------------------------------------------------------------------
# 3.  Virtual model: predict concentration t seconds ahead
# ---------------------------------------------------------------------------

def virtual_model_cp(ps: List[float], t: int, tci: dict) -> float:
    """Predict plasma conc (mcg/mL) after t seconds from exponential states."""
    lam = [tci['lam1'], tci['lam2'], tci['lam3']]
    return sum(ps[i] * (0.0 if lam[i] * t > 100 else math.exp(-lam[i] * t))
               for i in range(3))


def virtual_model_ce(es: List[float], t: int, tci: dict) -> float:
    """Predict effect-site conc (mcg/mL) after t seconds from exponential states."""
    lam = [tci['lam1'], tci['lam2'], tci['lam3'], tci['lam4']]
    return sum(es[i] * (0.0 if lam[i] * t > 100 else math.exp(-lam[i] * t))
               for i in range(4))


# ---------------------------------------------------------------------------
# 4.  find_peak  (translated from find_peak() in pharmacology.js)
# ---------------------------------------------------------------------------

def find_peak_time(
    current_time: int,
    rate: float,
    es: List[float],
    tci: dict,
) -> int:
    """
    Hill-climb to find the time at which Ce is maximum given current effect-site
    state *es* plus one unit-bolus period at *rate* mg/s.

    Translated from find_peak() in pharmacology.js.
    """
    e_udf = tci['e_udf']
    delta_s = tci['delta_seconds']
    current_time = max(delta_s, current_time)

    def ce_at(t: int) -> float:
        t = max(delta_s, t)
        eu = e_udf.get(t, 0.0)
        return virtual_model_ce(es, t, tci) + eu * rate

    current = ce_at(current_time)
    earlier = ce_at(current_time - 1)
    later   = ce_at(current_time + 1)

    for _ in range(1000):
        if current >= earlier and current >= later:
            break
        if current < earlier:
            if current_time <= delta_s:
                return current_time
            current_time -= 1
            later   = current
            current = earlier
            earlier = ce_at(current_time - 1)
        else:
            current_time += 1
            earlier = current
            current = later
            later   = ce_at(current_time + 1)

    return current_time


# ---------------------------------------------------------------------------
# 5.  Main simulation: STANPUMP CET
# ---------------------------------------------------------------------------

def simulate_stanpump_cet(
    pk_params: Dict[str, float],
    ce_target: float,
    duration_s: int,
    syringe_capacity_mg: float = 500.0,
    infusate_conc_mg_ml: float = 10.0,
    syringe_change_durations: Optional[List[int]] = None,
    delta_seconds: int = 1,
    max_rate_ml_h: float = 0.0,
    extend_until_ce: float = 0.0,
    rate_step_ml_h: float = 0.0,
) -> dict:
    """
    Simulate propofol TCI using the STANPUMP effect-site targeting algorithm.

    Parameters
    ----------
    pk_params : dict
        {vc, k10, k12, k13, k21, k31, ke0} — /min and litres.
    ce_target : float
        Target effect-site concentration (mcg/mL).
    duration_s : int
        Total simulation duration in seconds.
    syringe_capacity_mg : float
        Propofol per syringe (mg).  Default 500 mg (50 mL × 10 mg/mL).
    infusate_conc_mg_ml : float
        Propofol concentration in syringe (mg/mL).  Default 10.
    syringe_change_durations : list[int] or None
        Durations of successive syringe changes (seconds), applied cyclically.
        [180, 20] → 1st=180 s, 2nd=20 s, 3rd=180 s, …
    delta_seconds : int
        Bolus unit length and CPT lookahead (seconds).  SimTIVA uses 1.
    max_rate_ml_h : float
        Maximum infusion rate (ml/h); 0 = unlimited.
    extend_until_ce : float
        After duration_s, continue with rate = 0 until Ce falls below this
        threshold (mcg/mL).  0 = no extension.  Capped at 7200 extra seconds.
    rate_step_ml_h : float
        Quantization step for CPT maintenance rates (ml/h); 0 = continuous.
        E.g. 0.1 → rates floored to nearest 0.1 ml/h; below 0.1 rounds to 0.

    Returns
    -------
    dict with keys:
        time_s, time_min, rate_mgs, rate_mlh,
        cp, ce, dose_mg,
        total_dose_mg, syringe_changes, n_changes,
        ce_target, initial_bolus_mg, peak_time_s, tci,
        wakeup_time_s, wakeup_time_min,
        dose_5min_mg, dose_10min_mg, dose_30min_mg, dose_60min_mg, dose_end_mg,
        peak_rate_mlh, auc_rate_mlh_min, stabilization_time_min,
        rate_change_freq, zero_infusion_dur_s, infusion_cv
    """
    if syringe_change_durations is None:
        syringe_change_durations = [20, 180]

    def change_dur(n: int) -> int:
        return syringe_change_durations[n % len(syringe_change_durations)]

    # ------------------------------------------------------------------
    # Precompute coefficients
    # ------------------------------------------------------------------
    tci = calculate_udfs(pk_params, delta_seconds)

    l1, l2, l3, l4 = tci['l1'], tci['l2'], tci['l3'], tci['l4']
    p_coef = tci['p_coef']
    e_coef = tci['e_coef']
    p_udf  = tci['p_udf']
    e_udf  = tci['e_udf']
    peak_time = tci['peak_time']

    max_rate_mgs = (max_rate_ml_h * infusate_conc_mg_ml / 3600.0
                    if max_rate_ml_h > 0 else 0.0)
    rate_step_mgs = (rate_step_ml_h * infusate_conc_mg_ml / 3600.0
                     if rate_step_ml_h > 0 else 0.0)

    # ------------------------------------------------------------------
    # CET induction: compute initial bolus from zero state
    # (STANPUMP lines 2002-2020 / deliver_cet_real logic)
    # ------------------------------------------------------------------
    es0 = [0.0, 0.0, 0.0, 0.0]   # zero initial state

    temp_peak  = peak_time
    trial_rate = ce_target / e_udf[temp_peak]           # first estimate
    temp_peak  = find_peak_time(temp_peak, trial_rate, es0, tci)

    min_dif = ce_target * 0.0001
    for _ in range(100):
        predicted_ce = virtual_model_ce(es0, temp_peak, tci)
        trial_rate   = (ce_target - predicted_ce) / e_udf[temp_peak]
        new_peak     = find_peak_time(temp_peak, trial_rate, es0, tci)
        achieved_ce  = predicted_ce + e_udf[temp_peak] * trial_rate
        if abs(achieved_ce - ce_target) <= min_dif:
            temp_peak = new_peak
            break
        temp_peak = new_peak

    cet_bolus_mg = math.ceil(trial_rate)    # mg (1-second bolus at trial_rate mg/s)

    # Bolus duration: if max_rate set, administer over multiple seconds
    if max_rate_mgs > 0 and cet_bolus_mg / max_rate_mgs > 1.0:
        bolus_duration_s = math.ceil(cet_bolus_mg / max_rate_mgs)
        bolus_rate_mgs   = cet_bolus_mg / bolus_duration_s
        # Recompute Ce-peak time: simulate rate-limited bolus, then hill-climb to peak
        _ps = [0.0, 0.0, 0.0]
        _es = [0.0, 0.0, 0.0, 0.0]
        for _ in range(bolus_duration_s):
            _ps[0] = _ps[0] * l1 + p_coef[0] * bolus_rate_mgs * (1.0 - l1)
            _ps[1] = _ps[1] * l2 + p_coef[1] * bolus_rate_mgs * (1.0 - l2)
            _ps[2] = _ps[2] * l3 + p_coef[2] * bolus_rate_mgs * (1.0 - l3)
            _es[0] = _es[0] * l1 + e_coef[0] * bolus_rate_mgs * (1.0 - l1)
            _es[1] = _es[1] * l2 + e_coef[1] * bolus_rate_mgs * (1.0 - l2)
            _es[2] = _es[2] * l3 + e_coef[2] * bolus_rate_mgs * (1.0 - l3)
            _es[3] = _es[3] * l4 + e_coef[3] * bolus_rate_mgs * (1.0 - l4)
        prior = sum(_es)
        cpt_start_s = bolus_duration_s
        for _dt in range(1, 600):
            _val = virtual_model_ce(_es, _dt, tci)
            if _val < prior:
                cpt_start_s = bolus_duration_s + _dt - 1
                break
            prior = _val
        else:
            cpt_start_s = bolus_duration_s + 300
    else:
        bolus_duration_s = 1
        bolus_rate_mgs   = float(cet_bolus_mg)
        # Pause phase ends at temp_peak; CPT starts there
        cpt_start_s = temp_peak

    # ------------------------------------------------------------------
    # Output arrays  (index = second)
    # ------------------------------------------------------------------
    N = duration_s + 1
    time_arr = np.arange(N, dtype=float)
    rate_arr = np.zeros(N)
    cp_arr   = np.zeros(N)
    ce_arr   = np.zeros(N)
    dose_arr = np.zeros(N)

    # ------------------------------------------------------------------
    # Exponential state vectors
    # ------------------------------------------------------------------
    ps = [0.0, 0.0, 0.0]       # plasma
    es = [0.0, 0.0, 0.0, 0.0]  # effect site

    # ------------------------------------------------------------------
    # Syringe & phase state
    # ------------------------------------------------------------------
    syringe_mg    = syringe_capacity_mg
    n_changes     = 0
    change_cd     = 0   # countdown seconds in current change
    syringe_events: List[Tuple[int, int]] = []   # (time_s, duration_s)

    total_dose = 0.0
    cpt_rate   = 0.0   # current CPT maintenance rate (mg/s)

    # ------------------------------------------------------------------
    # Run second-by-second
    # ------------------------------------------------------------------
    for t in range(duration_s):
        # Record state at start of this second
        cp_arr[t]   = sum(ps)
        ce_arr[t]   = sum(es)
        rate_arr[t] = 0.0          # will be filled below
        dose_arr[t] = total_dose

        # --- Determine infusion rate for this second ---
        if change_cd > 0:
            # Syringe change in progress: pump off
            r = 0.0
            change_cd -= 1
            if change_cd == 0:
                # Change complete — refill and force CPT recompute next second
                syringe_mg = syringe_capacity_mg
                cpt_rate   = 0.0   # will be recomputed

        elif t < bolus_duration_s:
            # CET induction bolus
            r = bolus_rate_mgs

        elif t < cpt_start_s:
            # Post-bolus pause while Ce is still rising
            r = 0.0

        else:
            # CPT maintenance: compute required rate every delta_seconds
            if t == cpt_start_s or (t - cpt_start_s) % delta_seconds == 0:
                predicted_cp = virtual_model_cp(ps, delta_seconds, tci)
                if ce_target > predicted_cp:
                    cpt_rate = (ce_target - predicted_cp) / p_udf[delta_seconds]
                else:
                    cpt_rate = 0.0
                if max_rate_mgs > 0:
                    cpt_rate = min(cpt_rate, max_rate_mgs)
                if rate_step_mgs > 0 and cpt_rate > 0:
                    cpt_rate = math.floor(cpt_rate / rate_step_mgs) * rate_step_mgs
            r = cpt_rate

        # --- Syringe constraint ---
        if r > 0 and change_cd == 0:
            if r > syringe_mg:
                # Use what remains, then trigger syringe change next second
                r = syringe_mg
                syringe_mg = 0.0
                dur = change_dur(n_changes)
                change_cd = dur
                n_changes += 1
                syringe_events.append((t + 1, dur))
            else:
                syringe_mg -= r

        # Store this second's rate
        rate_arr[t] = r
        total_dose += r

        # --- Advance exponential states by 1 second ---
        ps[0] = ps[0] * l1 + p_coef[0] * r * (1.0 - l1)
        ps[1] = ps[1] * l2 + p_coef[1] * r * (1.0 - l2)
        ps[2] = ps[2] * l3 + p_coef[2] * r * (1.0 - l3)

        es[0] = es[0] * l1 + e_coef[0] * r * (1.0 - l1)
        es[1] = es[1] * l2 + e_coef[1] * r * (1.0 - l2)
        es[2] = es[2] * l3 + e_coef[2] * r * (1.0 - l3)
        es[3] = es[3] * l4 + e_coef[3] * r * (1.0 - l4)

    # Final timepoint
    cp_arr[duration_s]   = sum(ps)
    ce_arr[duration_s]   = sum(es)
    rate_arr[duration_s] = 0.0
    dose_arr[duration_s] = total_dose

    # ------------------------------------------------------------------
    # Washout extension: run with rate = 0 until Ce < extend_until_ce
    # ------------------------------------------------------------------
    wakeup_time_s = None
    if extend_until_ce > 0.0 and sum(es) > extend_until_ce:
        ext_t, ext_r, ext_cp, ext_ce, ext_dose = [], [], [], [], []
        wps = list(ps)
        wes = list(es)
        for _w in range(1, 7201):
            wps[0] *= l1;  wps[1] *= l2;  wps[2] *= l3
            wes[0] *= l1;  wes[1] *= l2;  wes[2] *= l3;  wes[3] *= l4
            _t = duration_s + _w
            ext_t.append(float(_t))
            ext_r.append(0.0)
            ext_cp.append(sum(wps))
            ext_ce.append(sum(wes))
            ext_dose.append(total_dose)
            if sum(wes) < extend_until_ce:
                wakeup_time_s = _t
                break
        time_arr = np.concatenate([time_arr, np.array(ext_t)])
        rate_arr = np.concatenate([rate_arr, np.zeros(len(ext_t))])
        cp_arr   = np.concatenate([cp_arr,   np.array(ext_cp)])
        ce_arr   = np.concatenate([ce_arr,   np.array(ext_ce)])
        dose_arr = np.concatenate([dose_arr, np.array(ext_dose)])

    # ------------------------------------------------------------------
    # Endpoint metrics
    # ------------------------------------------------------------------
    def _dose_at(t: int) -> float:
        return float(dose_arr[min(t, len(dose_arr) - 1)])

    dose_5min_mg  = _dose_at(300)
    dose_10min_mg = _dose_at(600)
    dose_30min_mg = _dose_at(1800)
    dose_60min_mg = _dose_at(3600)
    dose_end_mg   = _dose_at(duration_s)

    rate_mlh_proc  = rate_arr[:duration_s + 1] * 3600.0 / infusate_conc_mg_ml
    time_min_proc  = time_arr[:duration_s + 1] / 60.0
    peak_rate_mlh_val = float(np.max(rate_mlh_proc))
    auc_rate_mlh_min  = float(np.trapezoid(rate_mlh_proc, time_min_proc))

    # Stabilization: first 5-min window after maintenance starts with CV < 0.15
    maint_rate_mlh = rate_arr[cpt_start_s:duration_s + 1] * 3600.0 / infusate_conc_mg_ml
    stab_min = None
    for _i in range(max(0, len(maint_rate_mlh) - 300)):
        chunk = maint_rate_mlh[_i:_i + 300]
        _m = np.mean(chunk)
        if _m > 0 and np.std(chunk) / _m < 0.15:
            stab_min = (cpt_start_s + _i) / 60.0
            break

    # Secondary: rate-change frequency, zero-infusion duration, CV
    maint_seq = rate_arr[cpt_start_s:duration_s] * 3600.0 / infusate_conc_mg_ml
    _n_maint_changes = int(np.sum(np.diff(maint_seq) != 0))
    maint_dur_min = (duration_s - cpt_start_s) / 60.0
    rate_chg_freq = _n_maint_changes / maint_dur_min if maint_dur_min > 0 else 0.0
    zero_dur_s    = int(np.sum(maint_seq == 0))
    _nz = maint_seq[maint_seq > 0]
    infusion_cv   = float(np.std(_nz) / np.mean(_nz)) if len(_nz) > 0 else 0.0

    return {
        'time_s':         time_arr,
        'time_min':       time_arr / 60.0,
        'rate_mgs':       rate_arr,
        'rate_mlh':       rate_arr * 3600.0 / infusate_conc_mg_ml,
        'cp':             cp_arr,
        'ce':             ce_arr,
        'dose_mg':        dose_arr,
        'total_dose_mg':  total_dose,
        'syringe_changes': syringe_events,
        'n_changes':      n_changes,
        'ce_target':      ce_target,
        'initial_bolus_mg':    cet_bolus_mg,
        'peak_time_s':         cpt_start_s,
        'tci':                 tci,
        'wakeup_time_s':       wakeup_time_s,
        'wakeup_time_min':     wakeup_time_s / 60.0 if wakeup_time_s is not None else None,
        'wakeup_ce_threshold': extend_until_ce,
        # Primary endpoints
        'dose_5min_mg':           dose_5min_mg,
        'dose_10min_mg':          dose_10min_mg,
        'dose_30min_mg':          dose_30min_mg,
        'dose_60min_mg':          dose_60min_mg,
        'dose_end_mg':            dose_end_mg,
        'peak_rate_mlh':          peak_rate_mlh_val,
        'auc_rate_mlh_min':       auc_rate_mlh_min,
        'stabilization_time_min': stab_min,
        # Secondary endpoints
        'rate_change_freq':    rate_chg_freq,
        'zero_infusion_dur_s': zero_dur_s,
        'infusion_cv':         infusion_cv,
    }
