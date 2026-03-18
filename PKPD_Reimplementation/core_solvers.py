import numpy as np
from scipy.integrate import solve_ivp

# --- A. Differential Equations Numerical Solver ---
def run_ode_solver(params, infusion_rate_func, t_eval):
    # Unwrap PK params
    k10, k12, k13 = params['k10'], params['k12'], params['k13']
    k21, k31, ke0 = params['k21'], params['k31'], params['ke0']
    vc = params['vc']

    def pk_model_ode(t, y):
        x1, x2, x3, ce = y
        u = infusion_rate_func(t)
        
        dx1_dt = -(k10 + k12 + k13) * x1 + k21 * x2 + k31 * x3 + u
        dx2_dt = k12 * x1 - k21 * x2
        dx3_dt = k13 * x1 - k31 * x3
        
        cp = x1 / vc
        dce_dt = ke0 * (cp - ce)
        
        return [dx1_dt, dx2_dt, dx3_dt, dce_dt]

    y0 = [0.0, 0.0, 0.0, 0.0]
    span = [t_eval[0], t_eval[-1]]
    
    sol = solve_ivp(pk_model_ode, span, y0, t_eval=t_eval, method='Radau', rtol=1e-6, atol=1e-8)
    
    cp_ode = sol.y[0] / vc
    ce_ode = sol.y[3]
    return cp_ode, ce_ode


# --- B. SimTIVA Exact Algebraic Analytical Solver ---
def run_algebraic_solver(params, infusion_rate_func, t_eval):
    k10, k12, k13 = params['k10'], params['k12'], params['k13']
    k21, k31, ke0 = params['k21'], params['k31'], params['ke0']
    vc = params['vc']

    # Roots of characteristic polynomial (cube() translated accurately)
    a0 = k10 * k21 * k31
    a1 = k10 * k31 + k21 * k31 + k21 * k13 + k10 * k21 + k31 * k12
    a2 = k10 + k12 + k13 + k21 + k31

    p = a1 - (a2**2 / 3.0)
    q = (2 * a2**3 / 27.0) - (a1 * a2 / 3.0) + a0
    r1 = np.sqrt(-(p**3) / 27.0)
    phi = (-q / 2.0) / r1
    phi = np.clip(phi, -1, 1)
    phi = np.arccos(phi) / 3.0
    r1 = 2.0 * np.exp(np.log(r1) / 3.0)

    toradian = np.arcsin(1.0) * 2.0 / 180.0
    roots = np.zeros(3)
    roots[0] = -(np.cos(phi) * r1 - a2 / 3.0)
    roots[1] = -(np.cos(phi + 120.0 * toradian) * r1 - a2 / 3.0)
    roots[2] = -(np.cos(phi + 240.0 * toradian) * r1 - a2 / 3.0)

    # Sort roots ascending (most negative first) 
    roots = np.sort(roots)

    # Lambda assignment (positive magnitudes)
    lam1 = roots[0]
    lam2 = roots[1]
    lam3 = roots[2]
    lam4 = ke0 

    # Coefficients logic array 
    p_coef = np.zeros(3)
    p_coef[0] = (k21 - lam1) * (k31 - lam1) / ((lam1 - lam2) * (lam1 - lam3) * vc * lam1)
    p_coef[1] = (k21 - lam2) * (k31 - lam2) / ((lam2 - lam1) * (lam2 - lam3) * vc * lam2)
    p_coef[2] = (k21 - lam3) * (k31 - lam3) / ((lam3 - lam2) * (lam3 - lam1) * vc * lam3)

    e_coef = np.zeros(4)
    e_coef[0] = p_coef[0] / (ke0 - lam1) * ke0
    e_coef[1] = p_coef[1] / (ke0 - lam2) * ke0
    e_coef[2] = p_coef[2] / (ke0 - lam3) * ke0
    e_coef[3] = (ke0 - k21) * (ke0 - k31) / ((lam1 - ke0) * (lam2 - ke0) * (lam3 - ke0) * vc)

    # Simulation Sequence
    dt = t_eval[1] - t_eval[0]
    cp_alg = []
    ce_alg = []
    p_state = np.zeros(3)
    e_state = np.zeros(4)

    dec_l1 = np.exp(-lam1 * dt)
    dec_l2 = np.exp(-lam2 * dt)
    dec_l3 = np.exp(-lam3 * dt)
    dec_l4 = np.exp(-lam4 * dt)

    for t in t_eval:
        u = infusion_rate_func(t)
        
        # Store snapshot BEFORE integration addition (aligns cleanly with ODE evaluation) 
        cp_alg.append(np.sum(p_state))
        ce_alg.append(np.sum(e_state))
        
        # Exponential State Recursion
        p_state[0] = p_state[0] * dec_l1 + p_coef[0] * u * (1 - dec_l1)
        p_state[1] = p_state[1] * dec_l2 + p_coef[1] * u * (1 - dec_l2)
        p_state[2] = p_state[2] * dec_l3 + p_coef[2] * u * (1 - dec_l3)
        
        e_state[0] = e_state[0] * dec_l1 + e_coef[0] * u * (1 - dec_l1)
        e_state[1] = e_state[1] * dec_l2 + e_coef[1] * u * (1 - dec_l2)
        e_state[2] = e_state[2] * dec_l3 + e_coef[2] * u * (1 - dec_l3)
        e_state[3] = e_state[3] * dec_l4 + e_coef[3] * u * (1 - dec_l4)

    return np.array(cp_alg), np.array(ce_alg)
