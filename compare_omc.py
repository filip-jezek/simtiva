import csv
import numpy as np
import matplotlib.pyplot as plt

print("Loading OpenModelica CSV Outputs via csv module...")

def load_csv_data(filepath):
    with open(filepath, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        
        # Find indices
        t_idx = header.index('"time"')
        try:
            cp_idx = header.index('"central.cport.c"')
        except ValueError:
            cp_idx = [i for i, h in enumerate(header) if 'central.cport.c' in h][0]
            
        try:
            ce_idx = header.index('"Ce"')
        except ValueError:
            ce_idx = [i for i, h in enumerate(header) if 'Ce' in h][0]
            
        time_arr, cp_arr, ce_arr = [], [], []
        for row in reader:
            time_arr.append(float(row[t_idx]))
            cp_arr.append(float(row[cp_idx]))
            ce_arr.append(float(row[ce_idx]))
            
    return np.array(time_arr), np.array(cp_arr), np.array(ce_arr)

try:
    time_flawed, cp_flawed, ce_flawed = load_csv_data("PropofolMarsh_res.csv")
    time_fixed, cp_fixed, ce_fixed = load_csv_data("PropofolMarsh_fixedUnits_res.csv")
    
    # Original model used k-values without dividing by 60
    # Thus, everything processes 60 times faster mathematically (time is highly compressed).
    # And massFlow infused 60x the mass (per minute instead of per second)
    # The fix scaled time by assigning k / 60, meaning the timeline stretches by * 60.
    
    # We rescale the original flawed time vector explicitly:
    rescaled_time_flawed = time_flawed * 60.0
    
    # We rescale the concentrations. The mass flow pumped 60x too much mg per second. 
    # Therefore the total mass bounds should be 60x higher than reality.
    rescaled_cp_flawed = cp_flawed / 60.0
    rescaled_ce_flawed = ce_flawed / 60.0

    print("Interpolating data to evaluate strict numerical mapping RMSE...")
    # Because Modelica DASSL outputs dynamic time steps, we must interpolate to compare bounds
    cp_flawed_interp = np.interp(time_fixed, rescaled_time_flawed, rescaled_cp_flawed)
    ce_flawed_interp = np.interp(time_fixed, rescaled_time_flawed, rescaled_ce_flawed)

    rmse_cp = np.sqrt(np.mean((cp_fixed - cp_flawed_interp)**2))
    rmse_ce = np.sqrt(np.mean((ce_fixed - ce_flawed_interp)**2))
    
    print(f"Rescaled RMSE (Cp): {rmse_cp:.8f} mcg/ml")
    print(f"Rescaled RMSE (Ce): {rmse_ce:.8f} mcg/ml")

    if rmse_cp < 1e-2:
        print("SUCCESS: The flawed model linearly transforms perfectly into the fixed unit model.")

    # Visualization
    plt.figure(figsize=(10, 6))
    plt.plot(time_fixed, cp_fixed, label="Cp (Fixed SI Units)", color="blue", linewidth=3)
    plt.plot(rescaled_time_flawed, rescaled_cp_flawed, label="Cp (Flawed Extrapolated)", color="cyan", linestyle="dashed")
    
    plt.plot(time_fixed, ce_fixed, label="Ce (Fixed SI Units)", color="red", linewidth=3)
    plt.plot(rescaled_time_flawed, rescaled_ce_flawed, label="Ce (Flawed Extrapolated)", color="orange", linestyle="dashed")
    
    plt.title("Propofol Marsh OpenModelica Analysis: Scaling Proof")
    plt.xlabel("Time (seconds)")
    plt.ylabel("Concentration (mcg/ml)")
    plt.legend()
    plt.grid(True)
    plt.savefig("omc_scaling_proof.png")
    print("Saved plot to omc_scaling_proof.png")
    
except Exception as e:
    print(f"Error occurred analyzing files: {e}")
