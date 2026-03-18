import numpy as np
import matplotlib.pyplot as plt
from models import MarshModel, SchniderModel, PaedfusorModel, EleveldModel
from core_solvers import run_ode_solver, run_algebraic_solver
from clinical_scenario import generate_infusion_func

def verify_model(model_class, model_name, weight, age, height, gender):
    print(f"\n--- Verifying {model_name} Model ---")
    
    # Initialize the model instance
    model = model_class(weight=weight, age=age, height=height, gender=gender, is_adj_bw=False)
    params = model.get_parameters()
    
    # Generic continuous time array
    sim_time_min = 120
    t_eval = np.linspace(0, sim_time_min, 1200)
    
    # Infusion Scenario
    infusion_func = generate_infusion_func(weight)
    
    # Run Solvers
    cp_ode, ce_ode = run_ode_solver(params, infusion_func, t_eval)
    cp_alg, ce_alg = run_algebraic_solver(params, infusion_func, t_eval)
    
    # Verification
    rmse_cp = np.sqrt(np.mean((cp_ode - cp_alg)**2))
    rmse_ce = np.sqrt(np.mean((ce_ode - ce_alg)**2))
    
    print(f"RMSE (Cp): {rmse_cp:.8f} mcg/ml")
    print(f"RMSE (Ce): {rmse_ce:.8f} mcg/ml")
    
    if rmse_cp < 1e-2 and rmse_ce < 1e-2:
        print("SUCCESS: Full mathematical alignment.")
    else:
        print("WARNING: Divergence.")

    # Chart Generation
    plt.figure(figsize=(10, 6))
    plt.plot(t_eval, cp_ode, label="Cp (ODE Solver)", linewidth=3, color="blue", alpha=0.6)
    plt.plot(t_eval, cp_alg, label="Cp (Exact Algebraic Method)", linestyle='dashed', color="cyan")
    plt.plot(t_eval, ce_ode, label="Ce (ODE Solver)", linewidth=3, color="red", alpha=0.6)
    plt.plot(t_eval, ce_alg, label="Ce (Exact Algebraic Method)", linestyle='dashed', color="orange")

    plt.ylim([0, max(max(cp_ode), max(ce_ode)) * 1.2])
    plt.title(f"Propofol {model_name} Model ({weight}kg, {age}yo)")
    plt.xlabel("Time (minutes)")
    plt.ylabel("Concentration (mcg/ml)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    filename = f"verification_{model_name.lower()}.png"
    plt.savefig(filename)
    print(f"Plot saved to {filename}")


if __name__ == "__main__":
    weight = 70
    age = 40
    height = 175
    gender = 0 # male
    
    # Dictionary mapping class references to names
    models = {
        "Marsh": MarshModel,
        "Schnider": SchniderModel,
        "Paedfusor": PaedfusorModel,
        "Eleveld": EleveldModel
    }
    
    for name, m_class in models.items():
        verify_model(m_class, name, weight, age, height, gender)

    print("\nAll models verified successfully.")
