import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from models.marsh import MarshModel
from models.schnider import SchniderModel
from models.paedfusor import PaedfusorModel
from models.eleveld import EleveldModel
from models.eleveld_nonmem import EleveldNonMem
from core_solvers import run_algebraic_solver

def get_infusion_func(dosing_rows):
    """
    Pre-computes tuples bounded into high-speed native Python lists for rapid 300k+ step integration.
    """
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


def plot_summary(all_rmse, patient_ids):
    """Generates a box plot of RMSE distribution per model — scales to any N patients."""
    model_names = list(all_rmse.keys())
    colors_map = {'Marsh': 'red', 'Schnider': 'blue', 'Paedfusor': 'green',
                  'Eleveld': 'purple', 'ENONMEM': 'orange', 'ELEVELD_PD': 'cyan'}

    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # --- Left panel: Box plot (N-agnostic) ---
    ax = axes[0]
    data = [all_rmse[m] for m in model_names]
    bp = ax.boxplot(data, patch_artist=True, notch=False, vert=True)
    for patch, name in zip(bp['boxes'], model_names):
        patch.set_facecolor(colors_map.get(name, 'gray'))
        patch.set_alpha(0.7)
    ax.set_xticklabels(model_names, rotation=15, ha='right')
    ax.set_ylabel('RMSE (mg/L)')
    ax.set_title(f'RMSE Distribution by Model (N={len(patient_ids)} patients)')
    ax.grid(True, axis='y', alpha=0.3)

    # Annotate median
    for i, name in enumerate(model_names):
        median = np.median(all_rmse[name])
        ax.text(i + 1, median + 0.05, f'{median:.2f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    # --- Right panel: Mean ± SD bar chart ---
    ax2 = axes[1]
    means = [np.mean(all_rmse[m]) for m in model_names]
    stds  = [np.std(all_rmse[m])  for m in model_names]
    x = np.arange(len(model_names))
    bars = ax2.bar(x, means, yerr=stds, capsize=5, alpha=0.8,
                   color=[colors_map.get(m, 'gray') for m in model_names])
    ax2.set_xticks(x)
    ax2.set_xticklabels(model_names, rotation=15, ha='right')
    ax2.set_ylabel('Mean RMSE ± SD (mg/L)')
    ax2.set_title('Mean RMSE ± SD by Model')
    ax2.grid(True, axis='y', alpha=0.3)
    for bar, m, s in zip(bars, means, stds):
        ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + s + 0.05,
                 f'{m:.2f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.tight_layout()
    plt.savefig('rmse_summary.png', dpi=300)
    plt.close()
    print("Summary plot saved to rmse_summary.png")

    print(f"\n{'Model':<15} {'Mean':>8} {'Median':>8} {'SD':>8} {'Min':>8} {'Max':>8}")
    print("-" * 55)
    for name in model_names:
        v = all_rmse[name]
        print(f"{name:<15} {np.mean(v):>8.3f} {np.median(v):>8.3f} {np.std(v):>8.3f} {np.min(v):>8.3f} {np.max(v):>8.3f}")


def main():
    parser = argparse.ArgumentParser(description='Simulate PK models for a range of patients.')
    parser.add_argument('--pid-start', type=int, default=1, help='First patient ID (inclusive)')
    parser.add_argument('--pid-end',   type=int, default=10, help='Last patient ID (inclusive)')
    parser.add_argument('--no-plots',  action='store_true', help='Skip per-patient plots (for large runs)')
    args = parser.parse_args()

    print("Loading NONMEM PK and PD Datasets...")
    dataset_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'supplementary_digital_content_1.txt')
    df = pd.read_csv(dataset_path, sep=r'\s+')
    
    pd_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'supplementary_digital_content_3.txt')
    df_pd = pd.read_csv(pd_path, sep=r'\s+')

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'simulated_patients')
    os.makedirs(out_dir, exist_ok=True)

    all_rmse = {}   # {model_name: [rmse per patient]}
    patient_ids = []

    for pid in range(args.pid_start, args.pid_end + 1):
        p1 = df[df['ID'] == pid].sort_values('TIME').copy()
        if len(p1) == 0:
            continue

        pd_rows = df_pd[df_pd['ID'] == pid]
        ev1 = ev2 = ev3 = ecl = eq2 = eq3 = None
        if len(pd_rows) > 0:
            first_pd = pd_rows.iloc[0]
            ev1 = first_pd.get('EV1', None)
            ev2 = first_pd.get('EV2', None)
            ev3 = first_pd.get('EV3', None)
            ecl = first_pd.get('ECL', None)
            eq2 = first_pd.get('EQ2', None)
            eq3 = first_pd.get('EQ3', None)
        
        first_row = p1.iloc[0]
        weight = first_row['WGT']
        age    = first_row['AGE']
        height = first_row['HGT']
        gender = 0 if first_row['M1F2'] == 1 else 1 
        opioid = first_row['TECH'] == 2
        
        dosing_rows = p1[p1['EVID'] == 1]
        obs_rows    = p1[(p1['EVID'] == 0) & (p1['DV'] > 0)]
        if len(obs_rows) == 0:
            continue

        infusion_func = get_infusion_func(dosing_rows)
        
        max_time  = p1['TIME'].max() + 10 
        time_eval = np.linspace(0, max_time, int(max_time * 60))
        if len(time_eval) == 0:
            continue

        # --- ENONMEM_MAPPED derivation (uncomment together with the entry below) ---
        # _p_elv = EleveldModel(weight=weight, age=age, height=height, gender=gender, opioid=opioid).get_parameters()
        # evL1   = _p_elv['vc']
        # eL_cl  = _p_elv['k10'] * evL1
        # eL_q2  = _p_elv['k12'] * evL1;  eL_v2 = eL_q2 / _p_elv['k21']
        # eL_q3  = _p_elv['k13'] * evL1;  eL_v3 = eL_q3 / _p_elv['k31']

        models = {
            'Marsh':    MarshModel(weight=weight, age=age, height=height, gender=gender),
            'Schnider': SchniderModel(weight=weight, age=age, height=height, gender=gender),
            'Paedfusor':PaedfusorModel(weight=weight, age=age, height=height, gender=gender),
            'Eleveld':  EleveldModel(weight=weight, age=age, height=height, gender=gender, opioid=opioid),
            'ENONMEM':  EleveldNonMem(weight=weight, age=age, height=height, gender=gender, 
                                      pma=first_row['PMA'], tech=first_row['TECH'], a1v2=first_row['A1V2']),
            # 'ENONMEM_MAPPED': EleveldNonMem(  # re-enable to verify JS vs NONMEM alignment
            #     weight=weight, age=age, height=height, gender=gender,
            #     pma=first_row['PMA'], tech=first_row['TECH'], a1v2=first_row['A1V2'],
            #     ev1=evL1, ev2=eL_v2, ev3=eL_v3, ecl=eL_cl, eq2=eL_q2, eq3=eL_q3),
            # 'ELEVELD_PD': EleveldNonMem(  # Bayesian individual params from Dataset 3 — re-enable to compare
            #     weight=weight, age=age, height=height, gender=gender,
            #     pma=first_row['PMA'], tech=first_row['TECH'], a1v2=first_row['A1V2'],
            #     ev1=ev1, ev2=ev2, ev3=ev3, ecl=ecl, eq2=eq2, eq3=eq3),
        }

        colors = {'Marsh': 'red', 'Schnider': 'blue', 'Paedfusor': 'green',
                  'Eleveld': 'purple', 'ENONMEM': 'orange', 'ELEVELD_PD': 'cyan',
                  'ENONMEM_MAPPED': 'black'}

        overall_rmse = {}
        cp_curves = {}
        for name, model in models.items():
            params = model.get_parameters()
            cp, _  = run_algebraic_solver(params, infusion_func, time_eval)
            cp_curves[name] = cp

            sq_errors = []
            for _, obs in obs_rows.iterrows():
                t_idx = np.searchsorted(time_eval, obs['TIME'])
                if t_idx < len(time_eval):
                    sq_errors.append((cp[t_idx] - np.exp(obs['DV'])) ** 2)
            if sq_errors:
                rmse = np.sqrt(np.mean(sq_errors))
                overall_rmse[name] = rmse
                all_rmse.setdefault(name, []).append(rmse)

        patient_ids.append(pid)
        print(f"P{pid:04d}  " + "  ".join(f"{n}={v:.2f}" for n, v in overall_rmse.items()))

        if not args.no_plots:
            fig, ax1 = plt.subplots(figsize=(14, 8))

            for name, cp in cp_curves.items():
                is_pd = name in ['ENONMEM', 'ELEVELD_PD']
                ax1.plot(time_eval, cp, label=f'{name}', color=colors[name],
                         linewidth=is_pd and 3.0 or 2.0, alpha=is_pd and 0.9 or 0.5)

            ax1.scatter(obs_rows['TIME'], np.exp(obs_rows['DV']),
                        color='black', marker='X', s=120, label='Observed DV', zorder=10)

            rate_array = np.array([infusion_func(t) for t in time_eval])
            ax2 = ax1.twinx()
            ax2.fill_between(time_eval, rate_array, alpha=0.12, color='gray', step='post', label='Infusion Rate')
            ax2.set_ylabel('Infusion Rate (mg/min)', color='gray')
            ax2.tick_params(axis='y', labelcolor='gray')

            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax1.legend(lines + lines2, labels + labels2, loc='upper left', fontsize=9)

            rmse_text = "\n".join([f"{k}: {v:.2f}" for k, v in overall_rmse.items()])
            ax1.text(0.98, 0.98, rmse_text, transform=ax1.transAxes, fontsize=10,
                     va='top', ha='right', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

            plt.title(f"Patient {pid:04d}  ({weight}kg, {age:.0f}yo, {'F' if gender==1 else 'M'})")
            ax1.set_xlabel("Time (Minutes)")
            ax1.set_ylabel("Propofol Plasma Cp (mg/L)")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(os.path.join(out_dir, f'patient_{pid:04d}_overlay.png'), dpi=150)
            plt.close()

    if patient_ids:
        plot_summary(all_rmse, patient_ids)


if __name__ == "__main__":
    main()
