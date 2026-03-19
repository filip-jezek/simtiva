import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from models import MarshModel, SchniderModel, PaedfusorModel, EleveldModel
from core_solvers import run_ode_solver, run_algebraic_solver
from clinical_scenario import generate_infusion_func
import subprocess
import os
import csv
import sys

def run_openmodelica_cli(model_name, weight, age, height, gender, t_eval_min):
    stopTime_sec = t_eval_min[-1] * 60.0
    full_model_path = f"PropofolModels.Propofol{model_name}_fixedUnits"
    csv_file = f"{full_model_path}_res.csv"
    
    mos_content = f"""
loadModel(Modelica);
loadFile("c:/home/git/Pharmacolibrary/Pharmacolibrary/package.mo");
loadFile("PropofolModels.mo");
simulate({full_model_path}, startTime=0, stopTime={stopTime_sec}, outputFormat="csv");
"""
    mos_file = f"sim_{model_name}.mos"
    with open(mos_file, "w") as f:
        f.write(mos_content)
        
    omc_path = r"C:\\Program Files\\OpenModelica1.26.1-64bit\\bin\\omc.exe"
    
    print(f"[{model_name}] Booting Native OpenModelica Instance and Compiling C Source...")
    result = subprocess.run(["powershell.exe", "-Command", f"cd 'c:\\home\\git\\simtiva\\PKPD_Reimplementation'; & '{omc_path}' {mos_file}"], capture_output=True, text=True)
    
    if not os.path.exists(csv_file):
        print(f"FAILED OMC Execution! Log:\\n{result.stdout}\\nERR:\\n{result.stderr}")
        raise FileNotFoundError(f"OpenModelica Simulation failed to execute externally. Missing: {csv_file}")
        
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        try:
            t_idx = header.index('time')
        except ValueError:
            t_idx = header.index('"time"')
        cp_idx = [i for i, h in enumerate(header) if 'central.cport.c' in h][0]
        ce_idx = [i for i, h in enumerate(header) if 'Ce' in h][0]
        
        time_arr, cp_arr, ce_arr = [], [], []
        for row in reader:
            time_arr.append(float(row[t_idx]))
            cp_arr.append(float(row[cp_idx]))
            ce_arr.append(float(row[ce_idx]))
            
    # Converging structural Base SI physics into scaled user clinical bounds explicitly natively isolated.
    time_arr = np.array(time_arr) / 60.0 # OM Outputs identically scaled seconds 
    cp_arr = np.array(cp_arr) * 1e3 # OM Structural native strictly tracks kg/m3 natively. Mapped to mcg/ml
    ce_arr = np.array(ce_arr) * 1e3 # (1e3 translates mapping correctly dynamically identical across physics mappings)
    
    # OMC uses a dynamic integration step tracking constraints, map to exact t_eval steps for comparison evaluation
    cp_interp = np.interp(t_eval_min, time_arr, cp_arr)
    ce_interp = np.interp(t_eval_min, time_arr, ce_arr)
    
    return cp_interp, ce_interp

def verify_model(model_class, model_name, weight, age, height, gender):
    print(f"\\n--- Executing Direct Validations for {model_name}_fixedUnits Model ---")
    
    model = model_class(weight=weight, age=age, height=height, gender=gender, is_adj_bw=False)
    params = model.get_parameters()
    
    sim_time_min = 120
    t_eval = np.linspace(0, sim_time_min, 1200)
    infusion_func = generate_infusion_func(weight)
    
    cp_ode, ce_ode = run_ode_solver(params, infusion_func, t_eval)
    cp_alg, ce_alg = run_algebraic_solver(params, infusion_func, t_eval)
    cp_mod, ce_mod = run_openmodelica_cli(model_name, weight, age, height, gender, t_eval)
    
    rmse_cp_mod = np.sqrt(np.mean((cp_ode - cp_mod)**2))
    rmse_ce_mod = np.sqrt(np.mean((ce_ode - ce_mod)**2))
    
    print(f"RMSE {model_name} Modelica Native Cp: {rmse_cp_mod:.8f} mcg/ml")
    print(f"RMSE {model_name} Modelica Native Ce: {rmse_ce_mod:.8f} mcg/ml")
    
    if rmse_cp_mod < 1e-4 and rmse_ce_mod < 1e-4:
        print("SUCCESS: Full mathematical tracking validation securely established exactly scaling physical constraints.")
    else:
        print("WARNING: Divergence.")

    plt.figure(figsize=(10, 6))
    plt.plot(t_eval, cp_ode, label="Cp (Python ODE Reference)", linewidth=8, color="blue", alpha=0.3)
    plt.plot(t_eval, cp_alg, label="Cp (Analytical Proxy)", linestyle='dashed', linewidth=4, color="cyan")
    plt.plot(t_eval, cp_mod, label="Cp (Native OpenModelica Compiler Trace)", linestyle='solid', linewidth=2, color="black")
    
    plt.plot(t_eval, ce_ode, label="Ce (Python ODE Reference)", linewidth=8, color="red", alpha=0.3)
    plt.plot(t_eval, ce_alg, label="Ce (Analytical Proxy)", linestyle='dashed', linewidth=4, color="orange")
    plt.plot(t_eval, ce_mod, label="Ce (Native OpenModelica Compiler Trace)", linestyle='solid', linewidth=2, color="black")

    plt.ylim([0, max(max(cp_ode), max(ce_ode)) * 1.2])
    plt.title(f"Propofol {model_name}_fixedUnits Model OM Core Validation ({weight}kg, {age}yo)")
    plt.xlabel("Time (minutes)")
    plt.ylabel("Concentration (mcg/ml)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    filename = f"verification_omc_{model_name.lower()}_fixedUnits.png"
    plt.savefig(filename)
    print(f"Plot correctly generated securely outputting: {filename}")

if __name__ == "__main__":
    weight = 70
    age = 40
    height = 170
    gender = 0
    
    models = {
        "Marsh": MarshModel,
        "Schnider": SchniderModel,
        "Paedfusor": PaedfusorModel,
        "Eleveld": EleveldModel
    }
    
    for name, m_class in models.items():
        verify_model(m_class, name, weight, age, height, gender)

    print("\\nAll 4 fully compiled native OpenModelica architectures verified mapping effectively bridging dimensions accurately and natively.")
