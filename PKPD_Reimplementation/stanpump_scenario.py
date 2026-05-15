"""
Scenario: propofol TCI using STANPUMP CET algorithm.

Protocol (translated from Czech):
    Scenario — default:
      Ce target = 4 mcg/mL
      Procedure duration = 90 min
      Syringe: 50 mL of 1 % propofol  (500 mg, 10 mg/mL)
      Syringe changes (cyclic alternation):
        1st change: 180 s  (prolonged, e.g. nurse unavailable)
        2nd change:  20 s  (quick)
        3rd change: 180 s  … and so on
      Rate step: 0.1 mL/h  (pump resolution)

Patients simulated:
    1  Male   5 yr  18 kg 110 cm
    2  Female 10 yr 30 kg 140 cm
    3  Male   20 yr 80 kg 180 cm
    4  Female 20 yr 50 kg 150 cm
    5  Male   40 yr 80 kg 180 cm
    6  Female 40 yr 50 kg 150 cm
    7  Male   40 yr 120 kg 180 cm

Models: Marsh, Schnider, Paedfusor, Eleveld   (all without opioid for Eleveld)
"""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

from stanpump_tci import simulate_stanpump_cet
from models.marsh      import MarshModel
from models.schnider   import SchniderModel
from models.paedfusor  import PaedfusorModel
from models.eleveld    import EleveldModel


# ---------------------------------------------------------------------------
# Protocol constants
# ---------------------------------------------------------------------------
CE_TARGET          = 4.0        # mcg/mL
DURATION_MIN       = 90
DURATION_S         = DURATION_MIN * 60
SYRINGE_CAP_MG     = 500.0      # 50 mL × 10 mg/mL
INFUSATE_MG_ML     = 10.0       # 1 % propofol
CHANGE_DURATIONS   = [180, 20]  # cyclic: 1st=180 s, 2nd=20 s, 3rd=180 s, …
DELTA_SECONDS      = 1          # matching SimTIVA JS
WAKEUP_CE          = 1.5        # mcg/mL — propofol awakening threshold
RATE_STEP_ML_H     = 0.1        # pump resolution (mL/h)


# ---------------------------------------------------------------------------
# Patient definitions
# ---------------------------------------------------------------------------
PATIENTS = [
    {'id': 1, 'label': 'M 5yr 18kg 110cm',  'gender': 0, 'age':  5, 'weight': 18,  'height': 110},
    {'id': 2, 'label': 'F 10yr 30kg 140cm', 'gender': 1, 'age': 10, 'weight': 30,  'height': 140},
    {'id': 3, 'label': 'M 20yr 80kg 180cm', 'gender': 0, 'age': 20, 'weight': 80,  'height': 180},
    {'id': 4, 'label': 'F 20yr 50kg 150cm', 'gender': 1, 'age': 20, 'weight': 50,  'height': 150},
    {'id': 5, 'label': 'M 40yr 80kg 180cm', 'gender': 0, 'age': 40, 'weight': 80,  'height': 180},
    {'id': 6, 'label': 'F 40yr 50kg 150cm', 'gender': 1, 'age': 40, 'weight': 50,  'height': 150},
    {'id': 7, 'label': 'M 40yr 120kg 180cm','gender': 0, 'age': 40, 'weight': 120, 'height': 180},
]

MODEL_NAMES = ['Marsh', 'Schnider', 'Paedfusor', 'Eleveld']

# Colour cycle — one per patient
PATIENT_COLORS = [
    '#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
    '#ff7f00', '#a65628', '#f781bf',
]
# Line styles — one per model
MODEL_STYLES = ['-', '--', '-.', ':']

# Short labels for axes
MODEL_LABELS = {
    'Marsh':    'Marsh 1991',
    'Schnider': 'Schnider 1998',
    'Paedfusor':'Paedfusor 2003',
    'Eleveld':  'Eleveld 2018',
}


# ---------------------------------------------------------------------------
# Helper: build PK params for a given model & patient
# ---------------------------------------------------------------------------

def get_pk_params(model_name: str, patient: dict) -> dict:
    w, a, h, g = patient['weight'], patient['age'], patient['height'], patient['gender']

    if model_name == 'Marsh':
        return MarshModel(w, a, h, g).get_parameters()
    elif model_name == 'Schnider':
        return SchniderModel(w, a, h, g).get_parameters()
    elif model_name == 'Paedfusor':
        return PaedfusorModel(w, a, h, g).get_parameters()
    elif model_name == 'Eleveld':
        return EleveldModel(w, a, h, g, opioid=False).get_parameters()
    else:
        raise ValueError(f'Unknown model: {model_name}')


# ---------------------------------------------------------------------------
# Run all simulations
# ---------------------------------------------------------------------------

def run_all() -> dict:
    """Return nested dict results[model][patient_id] = simulation dict."""
    results = {m: {} for m in MODEL_NAMES}

    for model_name in MODEL_NAMES:
        print(f'\n=== Model: {model_name} ===')
        for pat in PATIENTS:
            pid = pat['id']
            try:
                pk = get_pk_params(model_name, pat)
                res = simulate_stanpump_cet(
                    pk_params             = pk,
                    ce_target             = CE_TARGET,
                    duration_s            = DURATION_S,
                    syringe_capacity_mg   = SYRINGE_CAP_MG,
                    infusate_conc_mg_ml   = INFUSATE_MG_ML,
                    syringe_change_durations = CHANGE_DURATIONS,
                    delta_seconds         = DELTA_SECONDS,
                    max_rate_ml_h         = 1800.0,
                    extend_until_ce       = WAKEUP_CE,
                    rate_step_ml_h        = RATE_STEP_ML_H,
                )
                results[model_name][pid] = res
                n_ch = res['n_changes']
                bols = res['initial_bolus_mg']
                peak = res['peak_time_s']
                dose = res['total_dose_mg']
                print(f"  Patient {pid} ({pat['label']}): "
                      f"bolus={bols} mg  peak={peak} s  "
                      f"total={dose:.1f} mg  syringe_changes={n_ch}")
            except Exception as exc:
                print(f"  Patient {pid} ({pat['label']}): ERROR — {exc}")
                results[model_name][pid] = None

    return results


# ---------------------------------------------------------------------------
# Plotting utilities
# ---------------------------------------------------------------------------

def _minutes(time_s: np.ndarray) -> np.ndarray:
    return time_s / 60.0


def plot_ce_comparison(results: dict, out_path: str) -> None:
    """4-panel plot: Ce vs time, one panel per model, all 7 patients."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=True, sharey=True)
    fig.suptitle(
        f'Effect-site concentration (Ce) — STANPUMP CET  |  '
        f'Target Ce = {CE_TARGET} mcg/mL  |  Duration = {DURATION_MIN} min',
        fontsize=13, fontweight='bold'
    )

    for ax, model_name in zip(axes.flat, MODEL_NAMES):
        ax.axhline(CE_TARGET, color='k', lw=1.2, ls='--', alpha=0.5, label='Target')
        ax.set_title(MODEL_LABELS[model_name], fontsize=11)
        ax.set_ylim(0, 7)
        ax.set_xlim(0, DURATION_MIN)
        ax.set_ylabel('Ce (mcg/mL)')
        ax.set_xlabel('Time (min)')
        ax.grid(True, alpha=0.3)

        for pat, color in zip(PATIENTS, PATIENT_COLORS):
            pid = pat['id']
            res = results[model_name].get(pid)
            if res is None:
                continue
            t_min = res['time_min']
            ce    = res['ce']
            # Downsample to every 10 s for plotting speed
            idx = np.arange(0, len(t_min), 10)
            ax.plot(t_min[idx], ce[idx], color=color, lw=1.3,
                    label=pat['label'])

        # Mark syringe change times for a representative patient (patient 5)
        rep = results[model_name].get(5)
        if rep:
            for (ts, dur) in rep['syringe_changes']:
                ax.axvspan(ts / 60, (ts + dur) / 60, color='gray',
                           alpha=0.15, zorder=0)

    # Shared legend (patients)
    legend_lines = [Line2D([0], [0], color=c, lw=2) for c in PATIENT_COLORS]
    legend_labels = [p['label'] for p in PATIENTS]
    legend_lines  += [Line2D([0], [0], color='k', lw=1.2, ls='--')]
    legend_labels += [f'Target {CE_TARGET} mcg/mL']
    fig.legend(legend_lines, legend_labels, loc='lower center', ncol=4,
               fontsize=8, bbox_to_anchor=(0.5, 0.0))

    plt.tight_layout(rect=[0, 0.07, 1, 1])
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


def plot_rate_comparison(results: dict, out_path: str) -> None:
    """4-panel plot: infusion rate (ml/h) vs time, one panel per model."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=True)
    fig.suptitle(
        'Infusion rate (ml/h) — STANPUMP CET  |  '
        '1% Propofol (10 mg/mL)  |  50 mL syringe',
        fontsize=13, fontweight='bold'
    )

    for ax, model_name in zip(axes.flat, MODEL_NAMES):
        ax.set_title(MODEL_LABELS[model_name], fontsize=11)
        ax.set_xlim(0, DURATION_MIN)
        ax.set_ylabel('Rate (ml/h)')
        ax.set_xlabel('Time (min)')
        ax.grid(True, alpha=0.3)

        for pat, color in zip(PATIENTS, PATIENT_COLORS):
            pid = pat['id']
            res = results[model_name].get(pid)
            if res is None:
                continue
            t_min    = res['time_min']
            rate_mlh = res['rate_mlh']
            idx = np.arange(0, len(t_min), 10)
            ax.plot(t_min[idx], rate_mlh[idx], color=color, lw=1.2,
                    label=pat['label'])

        # Shade syringe-change pauses (patient 5)
        rep = results[model_name].get(5)
        if rep:
            for (ts, dur) in rep['syringe_changes']:
                ax.axvspan(ts / 60, (ts + dur) / 60, color='gray',
                           alpha=0.2, zorder=0)

    legend_lines  = [Line2D([0], [0], color=c, lw=2) for c in PATIENT_COLORS]
    legend_labels = [p['label'] for p in PATIENTS]
    fig.legend(legend_lines, legend_labels, loc='lower center', ncol=4,
               fontsize=8, bbox_to_anchor=(0.5, 0.0))

    plt.tight_layout(rect=[0, 0.07, 1, 1])
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


def plot_dose_comparison(results: dict, out_path: str) -> None:
    """Bar chart: total propofol dose (mg) for each patient × model."""
    n_patients = len(PATIENTS)
    n_models   = len(MODEL_NAMES)
    x = np.arange(n_patients)
    width = 0.18

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.set_title('Total propofol dose over 90 min — STANPUMP CET', fontsize=13, fontweight='bold')

    hatch_patterns = ['', '//', 'xx', '..']
    bar_colors = ['#4c72b0', '#dd8452', '#55a868', '#c44e52']

    for mi, (model_name, hatch, color) in enumerate(zip(MODEL_NAMES, hatch_patterns, bar_colors)):
        doses = []
        for pat in PATIENTS:
            res = results[model_name].get(pat['id'])
            doses.append(res['total_dose_mg'] if res else 0.0)
        offset = (mi - (n_models - 1) / 2) * width
        bars = ax.bar(x + offset, doses, width, label=MODEL_LABELS[model_name],
                      color=color, hatch=hatch, edgecolor='white', alpha=0.85)
        for bar, d in zip(bars, doses):
            if d > 0:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 2,
                        f'{d:.0f}', ha='center', va='bottom', fontsize=7, rotation=45)

    ax.set_xlabel('Patient', fontsize=11)
    ax.set_ylabel('Total dose (mg)', fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels([p['label'] for p in PATIENTS], rotation=20, ha='right', fontsize=9)
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


def plot_patient_model_grid(results: dict, out_path: str) -> None:
    """
    7×4 grid: each row = patient, each column = model.
    Shows Ce (blue) and Cp (orange) together.
    """
    n_rows = len(PATIENTS)
    n_cols = len(MODEL_NAMES)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 22),
                             sharex=True, sharey=True)
    fig.suptitle(
        f'STANPUMP CET — Ce (blue) and Cp (orange)\n'
        f'Target Ce = {CE_TARGET} mcg/mL  |  90 min  |  50 mL × 1% propofol',
        fontsize=13, fontweight='bold'
    )

    for ri, pat in enumerate(PATIENTS):
        for ci, model_name in enumerate(MODEL_NAMES):
            ax = axes[ri][ci]
            pid = pat['id']
            res = results[model_name].get(pid)

            if ri == 0:
                ax.set_title(MODEL_LABELS[model_name], fontsize=9, fontweight='bold')
            if ci == 0:
                ax.set_ylabel(pat['label'], fontsize=8)

            ax.axhline(CE_TARGET, color='k', lw=0.8, ls='--', alpha=0.5)
            ax.set_xlim(0, DURATION_MIN)
            ax.set_ylim(0, 8)
            ax.grid(True, alpha=0.25)
            ax.tick_params(labelsize=7)

            if ri == n_rows - 1:
                ax.set_xlabel('min', fontsize=7)

            if res is None:
                ax.text(0.5, 0.5, 'N/A', transform=ax.transAxes,
                        ha='center', va='center', color='red', fontsize=10)
                continue

            t_min = res['time_min']
            idx = np.arange(0, len(t_min), 5)

            ax.plot(t_min[idx], res['ce'][idx], '#2166ac', lw=1.2, label='Ce')
            ax.plot(t_min[idx], res['cp'][idx], '#d6604d', lw=0.9, ls='--', alpha=0.7, label='Cp')

            # Gray bands for syringe changes
            for (ts, dur) in res['syringe_changes']:
                ax.axvspan(ts / 60, (ts + dur) / 60, color='gray', alpha=0.2, zorder=0)

            # Annotate bolus & n_changes
            ax.text(0.98, 0.97,
                    f'Bolus {res["initial_bolus_mg"]} mg\n{res["n_changes"]} changes\n{res["total_dose_mg"]:.0f} mg total',
                    transform=ax.transAxes, ha='right', va='top',
                    fontsize=6.5, bbox=dict(boxstyle='round,pad=0.2', fc='wheat', alpha=0.7))

    # Legend
    handles = [
        Line2D([0], [0], color='#2166ac', lw=1.5, label='Ce'),
        Line2D([0], [0], color='#d6604d', lw=1.0, ls='--', label='Cp'),
        Line2D([0], [0], color='k',       lw=0.8, ls='--', label='Target'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=3,
               fontsize=9, bbox_to_anchor=(0.5, 0.0))

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(out_path, dpi=130, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


def plot_model_overlay_by_patient(results: dict, out_path: str) -> None:
    """
    7-panel plot: one panel per patient, all 4 models overlaid (Ce vs time).
    Shows how model choice affects predicted Ce.
    """
    fig, axes = plt.subplots(4, 2, figsize=(14, 18))
    axes = axes.flat
    fig.suptitle(
        'Ce vs time — 4 models overlaid per patient\n'
        f'Target Ce = {CE_TARGET} mcg/mL  |  STANPUMP CET',
        fontsize=13, fontweight='bold'
    )

    for ax, pat in zip(axes, PATIENTS):
        ax.axhline(CE_TARGET, color='k', lw=1.0, ls='--', alpha=0.5)
        ax.set_title(pat['label'], fontsize=10, fontweight='bold')
        ax.set_xlim(0, DURATION_MIN)
        ax.set_ylim(0, 8)
        ax.set_xlabel('Time (min)', fontsize=9)
        ax.set_ylabel('Ce (mcg/mL)', fontsize=9)
        ax.grid(True, alpha=0.3)

        pid = pat['id']
        for model_name, ls in zip(MODEL_NAMES, MODEL_STYLES):
            res = results[model_name].get(pid)
            if res is None:
                continue
            t_min = res['time_min']
            idx   = np.arange(0, len(t_min), 5)
            ax.plot(t_min[idx], res['ce'][idx], ls=ls, lw=1.5,
                    label=MODEL_LABELS[model_name])

        ax.legend(fontsize=7, loc='upper right')

    # Hide the spare 8th panel
    try:
        list(axes)[7].set_visible(False)
    except (IndexError, StopIteration):
        pass

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


def plot_syringe_summary(results: dict, out_path: str) -> None:
    """Bar chart: number of syringe changes per patient × model."""
    n_patients = len(PATIENTS)
    n_models   = len(MODEL_NAMES)
    x = np.arange(n_patients)
    width = 0.18
    bar_colors = ['#4c72b0', '#dd8452', '#55a868', '#c44e52']

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))

    # --- Left: syringe changes ---
    ax1.set_title('Number of syringe changes', fontsize=12, fontweight='bold')
    for mi, (model_name, color) in enumerate(zip(MODEL_NAMES, bar_colors)):
        counts = []
        for pat in PATIENTS:
            res = results[model_name].get(pat['id'])
            counts.append(res['n_changes'] if res else 0)
        offset = (mi - (n_models - 1) / 2) * width
        ax1.bar(x + offset, counts, width, label=MODEL_LABELS[model_name],
                color=color, edgecolor='white', alpha=0.85)

    ax1.set_xticks(x)
    ax1.set_xticklabels([p['label'] for p in PATIENTS], rotation=20, ha='right', fontsize=9)
    ax1.set_ylabel('N syringe changes')
    ax1.legend(fontsize=8)
    ax1.grid(axis='y', alpha=0.3)
    ax1.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # --- Right: mg/kg total dose ---
    ax2.set_title('Total dose per kg body weight', fontsize=12, fontweight='bold')
    for mi, (model_name, color) in enumerate(zip(MODEL_NAMES, bar_colors)):
        dose_per_kg = []
        for pat in PATIENTS:
            res = results[model_name].get(pat['id'])
            if res:
                dose_per_kg.append(res['total_dose_mg'] / pat['weight'])
            else:
                dose_per_kg.append(0.0)
        offset = (mi - (n_models - 1) / 2) * width
        bars = ax2.bar(x + offset, dose_per_kg, width, label=MODEL_LABELS[model_name],
                       color=color, edgecolor='white', alpha=0.85)
        for bar, d in zip(bars, dose_per_kg):
            if d > 0:
                ax2.text(bar.get_x() + bar.get_width() / 2,
                         bar.get_height() + 0.2,
                         f'{d:.0f}', ha='center', va='bottom', fontsize=6, rotation=45)

    ax2.set_xticks(x)
    ax2.set_xticklabels([p['label'] for p in PATIENTS], rotation=20, ha='right', fontsize=9)
    ax2.set_ylabel('Dose (mg/kg)')
    ax2.legend(fontsize=8)
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary_table(results: dict) -> None:
    print('\n' + '=' * 120)
    print(f"{'STANPUMP CET  —  SUMMARY TABLE':^120}")
    print(f"{'Target Ce = 4 mcg/mL  |  90 min  |  50 mL 1% propofol (500 mg)':^120}")
    print('=' * 120)
    hdr = f"{'Patient':<22}"
    for m in MODEL_NAMES:
        hdr += f'  {m:>22}'
    print(hdr)
    print('-' * 120)

    metrics = [
        ('Bolus (mg)',       lambda r: f"{r['initial_bolus_mg']}"),
        ('Peak time (s)',    lambda r: f"{r['peak_time_s']}"),
        ('Total dose (mg)',  lambda r: f"{r['total_dose_mg']:.1f}"),
        ('Dose (mg/kg)',     lambda r: f"{r['total_dose_mg']/r.get('_weight',1):.1f}"),
        ('# changes',       lambda r: f"{r['n_changes']}"),
    ]

    for pat in PATIENTS:
        pid = pat['id']
        # inject weight for mg/kg
        for m in MODEL_NAMES:
            if results[m].get(pid):
                results[m][pid]['_weight'] = pat['weight']

        print(f"\n  {pat['label']}")
        for metric_name, metric_fn in metrics:
            row = f"    {metric_name:<18}"
            for m in MODEL_NAMES:
                res = results[m].get(pid)
                val = metric_fn(res) if res else 'N/A'
                row += f'  {val:>22}'
            print(row)

    print('=' * 120)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    outdir = os.path.join(os.path.dirname(__file__), 'figures')
    os.makedirs(outdir, exist_ok=True)

    print('Running STANPUMP CET simulations …')
    results = run_all()

    print('\nGenerating figures …')
    plot_ce_comparison(
        results,
        os.path.join(outdir, 'stanpump_ce_by_model.png'))
    plot_rate_comparison(
        results,
        os.path.join(outdir, 'stanpump_rate_by_model.png'))
    plot_dose_comparison(
        results,
        os.path.join(outdir, 'stanpump_dose_comparison.png'))
    plot_patient_model_grid(
        results,
        os.path.join(outdir, 'stanpump_grid_ce_cp.png'))
    plot_model_overlay_by_patient(
        results,
        os.path.join(outdir, 'stanpump_model_overlay.png'))
    plot_syringe_summary(
        results,
        os.path.join(outdir, 'stanpump_syringe_summary.png'))

    print_summary_table(results)
    print('\nDone.')
