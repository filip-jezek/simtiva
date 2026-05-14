#!/usr/bin/env python3
"""
STANPUMP CET clinical scenario for all 1033 Eleveld-dataset patients.

Outputs:
  patient_scenario/figures/patient_<ID:04d>.png  — per-patient time-course
  patient_scenario/patient_stats.csv             — one row per patient × model
  patient_scenario/population_both.png           — covariate scatter (all)
  patient_scenario/population_male.png           — covariate scatter (male)
  patient_scenario/population_female.png         — covariate scatter (female)
"""

import os, sys, time, traceback
from multiprocessing import Pool, cpu_count

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from scipy import stats as sp_stats

from stanpump_scenario import (
    CE_TARGET, DURATION_MIN, DURATION_S, WAKEUP_CE,
    SYRINGE_CAP_MG, INFUSATE_MG_ML, CHANGE_DURATIONS, DELTA_SECONDS,
    MODEL_NAMES, MODEL_LABELS, get_pk_params,
)
from stanpump_tci import simulate_stanpump_cet

# ── constants ─────────────────────────────────────────────────────────────────
MODEL_COLORS = {
    'Marsh':     '#2166ac',
    'Schnider':  '#d6604d',
    'Paedfusor': '#1a9641',
    'Eleveld':   '#8e44ad',
}
ADULT_MODELS = ['Marsh', 'Schnider', 'Eleveld']

DATA_FILE = os.path.normpath(
    os.path.join(os.path.dirname(__file__), '..', 'data',
                 'supplementary_digital_content_1.txt'))
OUT_DIR   = os.path.join(os.path.dirname(__file__), 'patient_scenario')
FIG_DIR   = os.path.join(OUT_DIR, 'figures')
STATS_CSV = os.path.join(OUT_DIR, 'patient_stats.csv')

STATS_COLS = [
    'patient_id', 'age', 'weight', 'height', 'bmi', 'sex',
    'model', 'bolus_mg', 'peak_time_s', 'total_dose_mg',
    'dose_per_kg', 'n_swaps', 'wakeup_time_min',
]

plt.rcParams.update({
    'font.family':    'sans-serif',
    'font.size':      8,
    'axes.linewidth': 0.7,
})


# ─────────────────────────────────────────────────────────────────────────────
# Simulation
# ─────────────────────────────────────────────────────────────────────────────

def _build_pat(row) -> dict:
    gender = 0 if int(row['M1F2']) == 1 else 1
    return {
        'id':     int(row['ID']),
        'age':    float(row['AGE']),
        'weight': float(row['WGT']),
        'height': float(row['HGT']),
        'gender': gender,
        'label':  (f'{"M" if gender == 0 else "F"} '
                   f'{float(row["AGE"]):.1f}yr '
                   f'{float(row["WGT"]):.0f}kg '
                   f'{float(row["HGT"]):.0f}cm'),
    }


def run_all_models(pat: dict) -> dict:
    results = {}
    for mn in MODEL_NAMES:
        try:
            pk = get_pk_params(mn, pat)
            results[mn] = simulate_stanpump_cet(
                pk_params                = pk,
                ce_target                = CE_TARGET,
                duration_s               = DURATION_S,
                syringe_capacity_mg      = SYRINGE_CAP_MG,
                infusate_conc_mg_ml      = INFUSATE_MG_ML,
                syringe_change_durations = CHANGE_DURATIONS,
                delta_seconds            = DELTA_SECONDS,
                max_rate_ml_h            = 1800.0,
                extend_until_ce          = WAKEUP_CE,
            )
        except Exception:
            results[mn] = None
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Stats extraction
# ─────────────────────────────────────────────────────────────────────────────

def extract_stats(pat: dict, results: dict) -> list:
    bmi = pat['weight'] / (pat['height'] / 100) ** 2
    sex = 'M' if pat['gender'] == 0 else 'F'
    rows = []
    for mn in MODEL_NAMES:
        res = results[mn]
        rows.append({
            'patient_id':      pat['id'],
            'age':             pat['age'],
            'weight':          pat['weight'],
            'height':          pat['height'],
            'bmi':             round(bmi, 2),
            'sex':             sex,
            'model':           mn,
            'bolus_mg':        res['initial_bolus_mg'] if res else None,
            'peak_time_s':     res['peak_time_s']      if res else None,
            'total_dose_mg':   round(res['total_dose_mg'], 1) if res else None,
            'dose_per_kg':     round(res['total_dose_mg'] / pat['weight'], 2) if res else None,
            'n_swaps':         res['n_changes']         if res else None,
            'wakeup_time_min': (round(res['wakeup_time_min'], 1)
                                if res and res['wakeup_time_min'] is not None
                                else None),
        })
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# Per-patient figure
# ─────────────────────────────────────────────────────────────────────────────

def _maint_rate_max(res: dict) -> float:
    peak  = res['peak_time_s']
    maint = res['rate_mlh'][peak:]
    if not len(maint) or not (maint > 0).any():
        return 10.0
    return float(np.percentile(maint[maint > 0], 99))


def plot_patient_fast(pat: dict, results: dict, out_path: str) -> None:
    pid = pat['id']

    x_max_min = DURATION_MIN + 5
    conc_max  = CE_TARGET
    rate_max  = 5.0
    for mn in ADULT_MODELS:
        res = results[mn]
        if res is None:
            continue
        x_max_min = max(x_max_min, float(res['time_min'][-1]) + 2)
        conc_max  = max(conc_max, float(res['ce'].max()),
                        float(np.percentile(res['cp'], 99.5)))
        rate_max  = max(rate_max, _maint_rate_max(res))

    conc_ymax = conc_max * 1.2
    rate_ymax = rate_max * 1.6

    fig = plt.figure(figsize=(12, 5.8))
    gs  = fig.add_gridspec(2, 1, height_ratios=[5, 1], hspace=0.06)
    ax  = fig.add_subplot(gs[0])
    ax_s = fig.add_subplot(gs[1])
    ax_s.axis('off')

    fig.suptitle(
        f'ID {pid:04d}: {pat["label"]}  ·  STANPUMP CET  ·  '
        f'Ce = {CE_TARGET} μg/mL  ·  {DURATION_MIN} min  ·  ≤1800 mL/h',
        fontsize=8.5, fontweight='bold',
    )

    ax2 = ax.twinx()
    ax.set_xlim(0, x_max_min)
    ax.set_ylim(0, conc_ymax)
    ax.set_xlabel('Time (min)', fontsize=7)
    ax.set_ylabel('Concentration (μg mL⁻¹)', fontsize=7)
    ax.tick_params(labelsize=7)
    ax.grid(True, alpha=0.18, zorder=0)
    ax.set_zorder(ax2.get_zorder() + 1)
    ax.patch.set_visible(False)

    ax2.set_ylim(0, rate_ymax)
    ax2.set_ylabel('Rate (ml h⁻¹)', fontsize=7, color='#777777')
    ax2.tick_params(labelsize=6, colors='#777777')
    ax2.yaxis.label.set_color('#777777')
    ax2.spines['right'].set_color('#cccccc')

    ax.axhline(CE_TARGET, color='#555555', lw=0.9, ls='--', zorder=4, alpha=0.6)
    ax.axhline(WAKEUP_CE, color='#aaaaaa', lw=0.7, ls=':',  zorder=3, alpha=0.7)
    ax.axvline(DURATION_MIN, color='#bbbbbb', lw=0.9, ls='--', zorder=3)
    ax.text(DURATION_MIN + 0.3, conc_ymax * 0.08,
            'off', ha='left', va='bottom', fontsize=6, color='#999999')

    step = 10
    stats_lines = []

    for mn in ADULT_MODELS:
        res   = results[mn]
        color = MODEL_COLORS[mn]
        if res is None:
            stats_lines.append(f'{MODEL_LABELS[mn].split()[0]:9s}  N/A')
            continue

        t_m  = res['time_min']
        t_ds = t_m[::step]
        r_ds = np.clip(res['rate_mlh'][::step], 0, rate_ymax)

        ax2.fill_between(t_ds, 0, r_ds, step='post', color=color, alpha=0.07, zorder=1)
        ax2.step(t_ds, r_ds, where='post', color=color, lw=1.5, alpha=0.50, zorder=2)
        ax.plot(t_ds, res['cp'][::step], color=color, lw=0.8, ls='--',
                alpha=0.40, zorder=5)
        ax.plot(t_ds, res['ce'][::step], color=color, lw=2.0, zorder=6,
                label=MODEL_LABELS[mn])

        for ts, _dur in res['syringe_changes']:
            ax.axvline(ts / 60.0, color=color, lw=1.0,
                       ls=(0, (5, 3)), alpha=0.70, zorder=3)

        wu = res.get('wakeup_time_min')
        wu_str = f'{wu:.0f} min' if wu is not None else '>2 h'
        if wu is not None and wu <= x_max_min:
            wi = min(int(np.searchsorted(t_m, wu)), len(t_m) - 1)
            ax.scatter(t_m[wi], res['ce'][wi], marker='D', s=28,
                       color=color, edgecolors='white', lw=0.6, zorder=8)

        stats_lines.append(
            f'{MODEL_LABELS[mn].split()[0]:9s}'
            f'  bolus {res["initial_bolus_mg"]:3d} mg'
            f'  peak {res["peak_time_s"]:3d} s'
            f'  {res["total_dose_mg"]:5.0f} mg'
            f' ({res["total_dose_mg"] / pat["weight"]:.1f} mg/kg)'
            f'  {res["n_changes"]} swaps'
            f'  wakeup {wu_str}'
        )

    ax.legend(loc='upper right', fontsize=7.5, framealpha=0.92, edgecolor='#cccccc')
    ax_s.text(0.01, 0.95, '\n'.join(stats_lines),
              transform=ax_s.transAxes,
              fontsize=7.0, family='monospace', va='top',
              bbox=dict(boxstyle='round,pad=0.3', fc='#f8f8f8',
                        ec='#cccccc', alpha=0.95))

    plt.savefig(out_path, dpi=100, bbox_inches='tight')
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Worker (called in a pool process)
# ─────────────────────────────────────────────────────────────────────────────

def _process_one(row_dict: dict) -> list:
    try:
        pat     = _build_pat(row_dict)
        results = run_all_models(pat)
        fig_path = os.path.join(FIG_DIR, f'patient_{pat["id"]:04d}.png')
        plot_patient_fast(pat, results, fig_path)
        return extract_stats(pat, results)
    except Exception:
        pid = row_dict.get('ID', '?')
        print(f'  ERROR patient {pid}: {traceback.format_exc(limit=2)}')
        return []


# ─────────────────────────────────────────────────────────────────────────────
# Population covariate plots (large N)
# ─────────────────────────────────────────────────────────────────────────────

def _regress_ci(x: np.ndarray, y: np.ndarray, x_line: np.ndarray):
    n = len(x)
    if n < 3:
        return None
    sl, ic, r, p, se = sp_stats.linregress(x.astype(float), y.astype(float))
    y_fit   = sl * x_line + ic
    t_crit  = sp_stats.t.ppf(0.975, df=n - 2)
    x_mean  = x.mean()
    SSx     = float(np.sum((x - x_mean) ** 2))
    ci      = (t_crit * se * np.sqrt(1.0 / n + (x_line - x_mean) ** 2 / SSx)
               if SSx > 0 else np.zeros_like(x_line))
    return y_fit, y_fit - ci, y_fit + ci, {'slope': sl, 'r': r, 'p': p, 'n': n}


def plot_population_covariate(df_stats: pd.DataFrame,
                               out_path: str, title_sex: str) -> None:
    covariates = [
        ('Age',    'Age (years)',   'age'),
        ('Weight', 'Weight (kg)',   'weight'),
        ('Height', 'Height (cm)',   'height'),
        ('BMI',    'BMI (kg m⁻²)', 'bmi'),
    ]
    legend_handles = [
        Line2D([0], [0], color=MODEL_COLORS[mn], lw=2.0,
               marker='o', markersize=5,
               markeredgecolor='#333333', markeredgewidth=0.4,
               label=MODEL_LABELS[mn])
        for mn in MODEL_NAMES
    ]

    n_pats = df_stats['patient_id'].nunique()
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    fig.suptitle(
        f'Population dosing covariate analysis  ·  STANPUMP CET  ·  '
        f'Ce = {CE_TARGET} μg/mL  ·  {DURATION_MIN} min  ·  '
        f'{title_sex}  (n = {n_pats} patients)\n'
        'Dots = individual patients  ·  Line = linear regression  ·  '
        'Band = 95% CI  ·  Inset = model n / mean±SD / r / p',
        fontsize=10, fontweight='bold',
    )

    for ax, (cov_title, xlabel, col) in zip(axes.flat, covariates):
        ax.set_title(cov_title, fontsize=9, fontweight='bold')
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel('Total dose per kg (mg kg⁻¹)', fontsize=8)
        ax.grid(True, alpha=0.18)
        ax.tick_params(labelsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        stat_lines = []
        for mn in MODEL_NAMES:
            color = MODEL_COLORS[mn]
            sub   = df_stats[(df_stats['model'] == mn)].dropna(
                subset=[col, 'dose_per_kg'])
            if sub.empty:
                continue
            x_arr = sub[col].values
            y_arr = sub['dose_per_kg'].values

            ax.scatter(x_arr, y_arr, color=color, s=5,
                       alpha=0.18, edgecolors='none', zorder=3)

            if len(x_arr) >= 3:
                x_line = np.linspace(x_arr.min(), x_arr.max(), 120)
                fit = _regress_ci(x_arr, y_arr, x_line)
                if fit:
                    y_fit, y_lo, y_hi, st = fit
                    ax.plot(x_line, y_fit, color=color, lw=2.0, alpha=0.85, zorder=4)
                    ax.fill_between(x_line, y_lo, y_hi,
                                    color=color, alpha=0.14, zorder=2)
                    p_str = 'p<0.001' if st['p'] < 0.001 else f'p={st["p"]:.3f}'
                    stat_lines.append(
                        f'{MODEL_LABELS[mn].split()[0]:9s}'
                        f' n={st["n"]}'
                        f' μ={y_arr.mean():.1f}±{y_arr.std():.1f}'
                        f' r={st["r"]:.2f} {p_str}'
                    )

        if stat_lines:
            ax.text(0.02, 0.97, '\n'.join(stat_lines),
                    transform=ax.transAxes, fontsize=6.2,
                    family='monospace', va='top',
                    bbox=dict(boxstyle='round,pad=0.3', fc='white',
                              ec='#cccccc', alpha=0.90))

    fig.legend(handles=legend_handles, loc='lower center', ncol=4,
               fontsize=9, frameon=True, bbox_to_anchor=(0.5, -0.02))
    plt.savefig(out_path, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    os.makedirs(FIG_DIR, exist_ok=True)

    # Load demographics (one row per patient)
    df_full = pd.read_csv(DATA_FILE, sep=' ')
    patients_df = (df_full.drop_duplicates('ID')
                   [['ID', 'AGE', 'WGT', 'HGT', 'M1F2']]
                   .reset_index(drop=True))
    n_pats = len(patients_df)
    print(f'Loaded {n_pats} patients  '
          f'(M={( patients_df["M1F2"]==1).sum()}  '
          f'F={(patients_df["M1F2"]==2).sum()})')
    print(f'Figures → {FIG_DIR}')
    print(f'Stats   → {STATS_CSV}')

    # ── Simulate + plot in parallel ───────────────────────────────────────
    rows_list = patients_df.to_dict('records')
    n_workers = max(1, cpu_count() - 1)
    print(f'\nRunning with {n_workers} worker processes …\n')
    t0 = time.time()

    all_stats = []
    with Pool(n_workers) as pool:
        for i, stat_rows in enumerate(
                pool.imap(_process_one, rows_list, chunksize=4), 1):
            all_stats.extend(stat_rows)
            if i % 100 == 0 or i == n_pats:
                elapsed = time.time() - t0
                eta     = elapsed / i * (n_pats - i)
                print(f'  {i:4d}/{n_pats}  '
                      f'{elapsed:5.0f}s elapsed  '
                      f'ETA {eta:.0f}s')

    # ── Write CSV ─────────────────────────────────────────────────────────
    df_stats = pd.DataFrame(all_stats, columns=STATS_COLS)
    df_stats.to_csv(STATS_CSV, index=False, float_format='%.3f')
    print(f'\nStats CSV: {len(df_stats)} rows  ({df_stats["patient_id"].nunique()} patients)')

    # Quick text summary
    summary_path = os.path.join(OUT_DIR, 'population_summary.txt')
    with open(summary_path, 'w') as f:
        f.write('STANPUMP CET — Population statistics\n')
        f.write(f'Protocol: Ce = {CE_TARGET} μg/mL, {DURATION_MIN} min, '
                f'≤1800 mL/h, syringe 500 mg\n')
        f.write(f'Patients: n = {df_stats["patient_id"].nunique()}\n\n')
        for mn in MODEL_NAMES:
            sub = df_stats[df_stats['model'] == mn].dropna(subset=['total_dose_mg'])
            if sub.empty:
                continue
            f.write(f'{"─"*60}\n{mn}\n{"─"*60}\n')
            for col, label in [
                ('bolus_mg',        'Bolus (mg)'),
                ('total_dose_mg',   'Total dose (mg)'),
                ('dose_per_kg',     'Dose/kg (mg/kg)'),
                ('n_swaps',         'Syringe swaps'),
                ('wakeup_time_min', 'Wake-up time (min)'),
            ]:
                v = sub[col].dropna()
                if v.empty:
                    continue
                f.write(f'  {label:<24}  '
                        f'median {v.median():7.1f}  '
                        f'mean {v.mean():7.1f} ± {v.std():5.1f}  '
                        f'[{v.min():.1f} – {v.max():.1f}]\n')
            f.write('\n')
    print(f'Summary  → {summary_path}')

    # ── Population covariate plots ─────────────────────────────────────────
    print('\nGenerating population covariate plots …')
    plot_population_covariate(
        df_stats,
        os.path.join(OUT_DIR, 'population_both.png'),
        'Male & Female combined',
    )
    plot_population_covariate(
        df_stats[df_stats['sex'] == 'M'],
        os.path.join(OUT_DIR, 'population_male.png'),
        'Male',
    )
    plot_population_covariate(
        df_stats[df_stats['sex'] == 'F'],
        os.path.join(OUT_DIR, 'population_female.png'),
        'Female',
    )

    total = time.time() - t0
    print(f'\nAll done in {total:.0f}s  ({total/60:.1f} min)')
