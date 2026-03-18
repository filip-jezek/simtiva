import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# --- 1. Parameters for Marsh Model (70kg patient) ---
weight = 70  # kg
vc = 0.228 * weight
k10 = 0.119
k12 = 0.112
k13 = 0.0419
k21 = 0.055
k31 = 0.0033
ke0 = 1.21

# Simulation time (minutes)
sim_time_min = 120
t_eval = np.linspace(0, sim_time_min, 1200)
dt = sim_time_min / 1200

# --- 2. Infusion Profile ---
bolus_dose_mg = 2 * weight
maint_rate_mg_h = 10 * weight
maint_rate_mg_min = maint_rate_mg_h / 60

def infusion_rate(t):
    """Returns infusion rate in mg/min at time t."""
    if t < 1.0:
        return bolus_dose_mg / 1.0
    elif t < 60.0:
        return maint_rate_mg_min
    else:
        return 0.0

# --- 3. ODE Model (Method A) ---
def pk_model_ode(t, y):
    x1, x2, x3, ce = y
    u = infusion_rate(t)
    
    # Differential Equations
    dx1_dt = -(k10 + k12 + k13) * x1 + k21 * x2 + k31 * x3 + u
    dx2_dt = k12 * x1 - k21 * x2
    dx3_dt = k13 * x1 - k31 * x3
    
    # Effect site
    cp = x1 / vc
    dce_dt = ke0 * (cp - ce)
    
    return [dx1_dt, dx2_dt, dx3_dt, dce_dt]

# --- 4. Run Numerical ODE Solver ---
print("Running ODE solver...")
y0 = [0.0, 0.0, 0.0, 0.0]
sol = solve_ivp(pk_model_ode, [0, sim_time_min], y0, t_eval=t_eval, method='Radau', rtol=1e-6, atol=1e-8)

cp_ode = sol.y[0] / vc
ce_ode = sol.y[3]

# --- 5. Analytical Algebraic Solver (Method B - SimTIVA equivalent) ---
# Roots of the characteristic polynomial for Marsh
print("Running algebraic exact solver...")
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

# Sort roots ascending (most negative first) to strictly match JS logic
roots = np.sort(roots)

# lambdas are the positive magnitudes of the roots (since cube() negates them)
lam1 = roots[0]
lam2 = roots[1]
lam3 = roots[2]
lam4 = ke0 

# Coefficients array building exactly like JS
p_coef = np.zeros(3)
p_coef[0] = (k21 - lam1) * (k31 - lam1) / ((lam1 - lam2) * (lam1 - lam3) * vc * lam1)
p_coef[1] = (k21 - lam2) * (k31 - lam2) / ((lam2 - lam1) * (lam2 - lam3) * vc * lam2)
p_coef[2] = (k21 - lam3) * (k31 - lam3) / ((lam3 - lam2) * (lam3 - lam1) * vc * lam3)

e_coef = np.zeros(4)
e_coef[0] = p_coef[0] / (ke0 - lam1) * ke0
e_coef[1] = p_coef[1] / (ke0 - lam2) * ke0
e_coef[2] = p_coef[2] / (ke0 - lam3) * ke0
e_coef[3] = (ke0 - k21) * (ke0 - k31) / ((lam1 - ke0) * (lam2 - ke0) * (lam3 - ke0) * vc)

# Simulate step-by-step
cp_alg = []
ce_alg = []
p_state = np.zeros(3)
e_state = np.zeros(4)

dec_L1 = np.exp(-lam1 * dt)
dec_L2 = np.exp(-lam2 * dt)
dec_L3 = np.exp(-lam3 * dt)
dec_L4 = np.exp(-lam4 * dt)

for t in t_eval:
    u = infusion_rate(t)
    
    # Store current snapshot (before adding new infusion for this dt, to align with continuous ODE state)
    cp_alg.append(np.sum(p_state))
    ce_alg.append(np.sum(e_state))
    
    # Exponential Update Loop exactly as in JS
    p_state[0] = p_state[0] * dec_L1 + p_coef[0] * u * (1 - dec_L1)
    p_state[1] = p_state[1] * dec_L2 + p_coef[1] * u * (1 - dec_L2)
    p_state[2] = p_state[2] * dec_L3 + p_coef[2] * u * (1 - dec_L3)
    
    e_state[0] = e_state[0] * dec_L1 + e_coef[0] * u * (1 - dec_L1)
    e_state[1] = e_state[1] * dec_L2 + e_coef[1] * u * (1 - dec_L2)
    e_state[2] = e_state[2] * dec_L3 + e_coef[2] * u * (1 - dec_L3)
    e_state[3] = e_state[3] * dec_L4 + e_coef[3] * u * (1 - dec_L4)

cp_alg = np.array(cp_alg)
ce_alg = np.array(ce_alg)

# --- 6. Verification & Plotting ---
rmse_cp = np.sqrt(np.mean((cp_ode - cp_alg)**2))
rmse_ce = np.sqrt(np.mean((ce_ode - ce_alg)**2))

print(f"\n--- Results Verification ---")
print(f"RMSE (Plasma Concentration): {rmse_cp:.8f} mcg/ml")
print(f"RMSE (Effect Concentration): {rmse_ce:.8f} mcg/ml")

if rmse_cp < 1e-2:
    print("SUCCESS: The ODE Model perfectly aligns with the exact analytical model.")
else:
    print("WARNING: Divergence detected.")

# Plot the graph overlay
plt.figure(figsize=(10, 6))
plt.plot(t_eval, cp_ode, label="Cp (ODE Solver)", linewidth=3, color="blue", alpha=0.6)
plt.plot(t_eval, cp_alg, label="Cp (Exact Algebraic JS Method)", linestyle='dashed', color="cyan")

plt.plot(t_eval, ce_ode, label="Ce (ODE Solver)", linewidth=3, color="red", alpha=0.6)
plt.plot(t_eval, ce_alg, label="Ce (Exact Algebraic JS Method)", linestyle='dashed', color="orange")

# Limit the axes to the ODE values so diverging JS values don't squish the chart
plt.ylim([0, max(max(cp_ode), max(ce_ode)) * 1.2])

plt.title("Propofol Marsh Model Simulation (70kg)")
plt.xlabel("Time (minutes)")
plt.ylabel("Concentration (mcg/ml)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("pk_verification.png")
print("Plot saved to pk_verification.png")

print(f"\n--- Diagnostics ---")
print(f"Lam:  lam1={lam1:.5f}, lam2={lam2:.5f}, lam3={lam3:.5f}, lam4={lam4:.5f}")
print(f"Dec:  dec1={dec_L1:.5e}, dec2={dec_L2:.5e}, dec3={dec_L3:.5e}, dec4={dec_L4:.5e}")
print(f"P_coef: {p_coef}")
print(f"E_coef: {e_coef}")

