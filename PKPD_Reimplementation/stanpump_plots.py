"""
Publication-quality plots for the STANPUMP CET scenario.

Figures produced:
  stanpump_summary.png              — 3-panel grouped bar chart
  stanpump_patient_<N>.png  (×7)   — adult models overlaid, extended to wake-up
  stanpump_population_both.png      — covariate scatter, all patients
  stanpump_population_male.png      — covariate scatter, male only
  stanpump_population_female.png    — covariate scatter, female only

Run as:  python3 stanpump_plots.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy import stats as sp_stats

from stanpump_scenario import (
    run_all, PATIENTS, MODEL_NAMES, MODEL_LABELS,
    CE_TARGET, DURATION_MIN, DURATION_S, WAKEUP_CE,
    SYRINGE_CAP_MG, INFUSATE_MG_ML, RATE_STEP_ML_H,
)

# ── permanent colour per model ────────────────────────────────────────────────
MODEL_COLORS = {
    'Marsh':     '#2166ac',
    'Schnider':  '#d6604d',
    'Paedfusor': '#1a9641',
    'Eleveld':   '#8e44ad',
}
ADULT_MODELS = ['Marsh', 'Schnider', 'Eleveld']
PAT_SHORT    = [p['label'] for p in PATIENTS]

plt.rcParams.update({
    'font.family':    'sans-serif',
    'font.size':      9,
    'axes.linewidth': 0.8,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
})


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _ds(arr: np.ndarray, step: int = 5) -> np.ndarray:
    return arr[::step]


def _maintenance_rate_max(res: dict) -> float:
    peak  = res['peak_time_s']
    maint = res['rate_mlh'][peak:]
    if len(maint) == 0 or not (maint > 0).any():
        return 10.0
    return float(np.percentile(maint[maint > 0], 99))


def _unified_scales(pid: int, results: dict, model_list=None) -> tuple:
    if model_list is None:
        model_list = MODEL_NAMES
    conc_max = CE_TARGET
    rate_max = 5.0
    for mn in model_list:
        res = results[mn].get(pid)
        if res is None:
            continue
        conc_max = max(conc_max,
                       float(res['ce'].max()),
                       float(np.percentile(res['cp'], 99.5)))
        rate_max = max(rate_max, _maintenance_rate_max(res))
    return conc_max * 1.2, rate_max * 1.6


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — Summary bar chart
# ─────────────────────────────────────────────────────────────────────────────

def plot_summary(results: dict, out_path: str) -> None:
    n_pat   = len(PATIENTS)
    n_mod   = len(ADULT_MODELS)
    x       = np.arange(n_pat)
    w       = 0.18
    offsets = np.linspace(-(n_mod - 1) / 2, (n_mod - 1) / 2, n_mod) * w

    fig, axes = plt.subplots(3, 1, figsize=(13, 11), constrained_layout=True)
    fig.suptitle(
        'STANPUMP CET — Protocol summary\n'
        f'Ce target = {CE_TARGET} μg/mL  ·  Duration = {DURATION_MIN} min  ·  '
        f'Syringe = 50 mL × 1% propofol  ·  Max rate = 1800 mL/h  ·  '
        f'Rate step = {RATE_STEP_ML_H} mL/h\n'
        'Syringe-change pauses: 1st = 180 s  ·  2nd = 20 s  ·  alternating (cyclic)',
        fontsize=11, fontweight='bold',
    )

    panels = [
        ('Number of syringe changes',    lambda r, p: r['n_changes'],                    True),
        ('Total propofol dose  (mg)',     lambda r, p: r['total_dose_mg'],                False),
        ('Total dose per kg  (mg kg⁻¹)', lambda r, p: r['total_dose_mg'] / p['weight'],  False),
    ]

    for ax, (title, getter, integer_y) in zip(axes, panels):
        ax.set_title(title, fontsize=10, fontweight='bold', loc='left', pad=4)
        ax.set_xticks(x)
        ax.set_xticklabels(PAT_SHORT, rotation=25, ha='right', fontsize=8)
        ax.grid(axis='y', color='#cccccc', linewidth=0.6, zorder=0)
        ax.set_axisbelow(True)

        for mi, mn in enumerate(ADULT_MODELS):
            vals = [getter(results[mn].get(p['id']), p)
                    if results[mn].get(p['id']) else 0.0
                    for p in PATIENTS]
            bars = ax.bar(
                x + offsets[mi], vals, w,
                label=MODEL_LABELS[mn], color=MODEL_COLORS[mn],
                edgecolor='white', linewidth=0.5, alpha=0.88, zorder=3,
            )
            for bar, v in zip(bars, vals):
                if v <= 0:
                    continue
                txt = f'{int(v)}' if integer_y else f'{v:.0f}'
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + ax.get_ylim()[1] * 0.005,
                    txt, ha='center', va='bottom',
                    fontsize=6.5, rotation=55, color='#333333',
                )

        if integer_y:
            ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax.set_ylim(0, ax.get_ylim()[1] * 1.18)

    handles = [Patch(facecolor=MODEL_COLORS[mn], label=MODEL_LABELS[mn],
                     edgecolor='white', alpha=0.88)
               for mn in ADULT_MODELS]
    fig.legend(handles=handles, loc='lower center', ncol=4,
               fontsize=9, frameon=True, bbox_to_anchor=(0.5, -0.02))

    plt.savefig(out_path, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 (×7) — Per-patient: adult models overlaid, extended to wake-up
# ─────────────────────────────────────────────────────────────────────────────

def plot_patient(pat: dict, results: dict, out_path: str) -> None:
    pid = pat['id']

    # x-axis: infusion duration + washout to latest wake-up + 5 min buffer
    x_max_min = DURATION_MIN + 5
    for mn in ADULT_MODELS:
        res = results[mn].get(pid)
        if res is None:
            continue
        t_end = float(res['time_min'][-1])
        x_max_min = max(x_max_min, t_end + 3)

    conc_ymax, rate_ymax = _unified_scales(pid, results, ADULT_MODELS)

    # Two-row layout: main plot (tall) + stats row (short, no axes frame)
    fig = plt.figure(figsize=(14, 7.0))
    gs  = fig.add_gridspec(2, 1, height_ratios=[5, 1], hspace=0.08)
    ax  = fig.add_subplot(gs[0])
    ax_stats = fig.add_subplot(gs[1])
    ax_stats.axis('off')

    fig.suptitle(
        f'Patient {pid}: {pat["label"]}  ·  STANPUMP CET  ·  '
        f'Ce = {CE_TARGET} μg/mL  ·  {DURATION_MIN} min  ·  '
        f'50 mL syringe  ·  ≤1800 mL/h  ·  '
        f'washout to Ce < {WAKEUP_CE} μg/mL',
        fontsize=10, fontweight='bold',
    )

    # secondary rate axis
    ax2 = ax.twinx()
    ax.set_xlim(0, x_max_min)
    ax.set_ylim(0, conc_ymax)
    ax.set_xlabel('Time (min)', fontsize=8)
    ax.set_ylabel('Concentration (μg mL⁻¹)', fontsize=8)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.2, zorder=0)
    ax.set_zorder(ax2.get_zorder() + 1)
    ax.patch.set_visible(False)

    ax2.set_ylim(0, rate_ymax)
    ax2.set_ylabel('Rate (ml h⁻¹)', fontsize=8, color='#777777')
    ax2.tick_params(labelsize=7, colors='#777777')
    ax2.yaxis.label.set_color('#777777')
    ax2.spines['right'].set_color('#bbbbbb')

    # Target line and wake-up threshold
    ax.axhline(CE_TARGET, color='#444444', lw=1.0, ls='--', zorder=4, alpha=0.7)
    ax.axhline(WAKEUP_CE, color='#999999', lw=0.8, ls=':', zorder=3, alpha=0.8)
    ax.text(x_max_min * 0.99, WAKEUP_CE + conc_ymax * 0.01,
            f'wake-up ≈ {WAKEUP_CE} μg/mL',
            ha='right', va='bottom', fontsize=7, color='#888888')

    # Infusion end marker
    ax.axvline(DURATION_MIN, color='#999999', lw=1.0, ls='--', zorder=3, alpha=0.6)
    ax.text(DURATION_MIN + 0.4, conc_ymax * 0.12,
            'infusion\noff', ha='left', va='bottom',
            fontsize=6.5, color='#888888', linespacing=1.3)

    step = 5
    stats_lines = []

    for mn in ADULT_MODELS:
        res = results[mn].get(pid)
        color = MODEL_COLORS[mn]

        if res is None:
            stats_lines.append(f'{MODEL_LABELS[mn].split()[0]:9s}  N/A')
            continue

        t_m     = res['time_min']
        t_ds    = _ds(t_m, step)
        cp_ds   = _ds(res['cp'], step)
        ce_ds   = _ds(res['ce'], step)
        rate_ds = np.clip(_ds(res['rate_mlh'], step), 0, rate_ymax)

        # Rate: very light fill + thick step line
        ax2.fill_between(t_ds, 0, rate_ds, step='post',
                         color=color, alpha=0.07, zorder=1)
        ax2.step(t_ds, rate_ds, where='post',
                 color=color, lw=1.8, alpha=0.55, zorder=2)

        # Cp (thin dashed) and Ce (thick solid)
        ax.plot(t_ds, cp_ds, color=color, lw=0.9, ls='--', alpha=0.45, zorder=5)
        ax.plot(t_ds, ce_ds, color=color, lw=2.2, zorder=6,
                label=MODEL_LABELS[mn])

        # Syringe changes: dashed verticals, no shading
        for i, (ts, dur) in enumerate(res['syringe_changes']):
            xv = ts / 60.0
            ax.axvline(xv, color=color, lw=1.1, ls=(0, (5, 3)), alpha=0.75, zorder=3)
            ax.text(xv + 0.3, conc_ymax * (0.85 - i * 0.12),
                    f'↕{dur}s', ha='left', va='top',
                    fontsize=6.0, color=color, alpha=0.85)

        # Wake-up diamond marker on Ce curve
        wu = res.get('wakeup_time_min')
        wu_str = f'{wu:.0f} min' if wu is not None else '>2 h'
        if wu is not None and wu <= x_max_min:
            wu_idx = np.searchsorted(t_m, wu)
            wu_idx = min(wu_idx, len(t_m) - 1)
            ax.scatter(t_m[wu_idx], res['ce'][wu_idx],
                       marker='D', s=40, color=color,
                       edgecolors='white', lw=0.8, zorder=8)

        # Stats block: 2 lines per model
        stab = res.get('stabilization_time_min')
        stab_str = f'{stab:.1f} min' if stab is not None else '>90 min'
        stats_lines.append(
            f'{MODEL_LABELS[mn].split()[0]:9s}'
            f'  bolus {res["initial_bolus_mg"]:3d} mg'
            f'  peak {res["peak_time_s"]:3d} s'
            f'  dose {res["total_dose_mg"]:5.0f} mg ({res["total_dose_mg"] / pat["weight"]:.1f} mg/kg)'
            f'  peak-rate {res["peak_rate_mlh"]:6.1f} mL/h'
            f'  AUC {res["auc_rate_mlh_min"]:6.0f} mL·min'
        )
        stats_lines.append(
            f'{"":9s}'
            f'  swaps {res["n_changes"]}'
            f'  zero {res["zero_infusion_dur_s"]:4d} s'
            f'  CV {res["infusion_cv"]:.2f}'
            f'  stab {stab_str}'
            f'  wake-up {wu_str}'
        )

    # Legend inside axes — upper right (data will be well below there during washout)
    ax.legend(loc='upper right', fontsize=8.5, framealpha=0.92,
              edgecolor='#cccccc')

    # Stats block in the dedicated lower row
    ax_stats.text(
        0.01, 0.95,
        '\n'.join(stats_lines),
        transform=ax_stats.transAxes,
        fontsize=8.0, family='monospace', va='top',
        bbox=dict(boxstyle='round,pad=0.4', fc='#f8f8f8',
                  ec='#cccccc', alpha=0.95),
    )

    plt.savefig(out_path, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — Population covariate scatter (3 sex variants)
# ─────────────────────────────────────────────────────────────────────────────

def _regress_ci(x: np.ndarray, y: np.ndarray, x_line: np.ndarray):
    """Linear regression + 95% CI band for the mean response.  Returns
    (y_fit, y_lo, y_hi, stats_dict) or None if < 3 points."""
    n = len(x)
    if n < 2:
        return None
    slope, intercept, r, p, se = sp_stats.linregress(x.astype(float),
                                                      y.astype(float))
    y_fit  = slope * x_line + intercept
    if n < 3:
        return y_fit, y_fit, y_fit, dict(slope=slope, r=r, p=p, n=n)
    t_crit = sp_stats.t.ppf(0.975, df=n - 2)
    x_mean = x.mean()
    SSx    = float(np.sum((x - x_mean) ** 2))
    if SSx == 0:
        ci = np.zeros_like(x_line)
    else:
        ci = t_crit * se * np.sqrt(1.0 / n + (x_line - x_mean) ** 2 / SSx)
    return y_fit, y_fit - ci, y_fit + ci, dict(slope=slope, r=r, p=p, n=n)


def _pop_covariate_panel(ax, patients, results_dict, cov_fn, xlabel) -> None:
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel('Total dose per kg (mg kg⁻¹)', fontsize=8)
    ax.grid(True, alpha=0.20)
    ax.tick_params(labelsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    stat_lines = []
    for mn in ADULT_MODELS:
        color  = MODEL_COLORS[mn]
        x_vals, y_vals = [], []
        for p in patients:
            res = results_dict[mn].get(p['id'])
            if res is None:
                continue
            x_vals.append(cov_fn(p))
            y_vals.append(res['total_dose_mg'] / p['weight'])
        if not x_vals:
            continue
        x_arr, y_arr = np.array(x_vals), np.array(y_vals)

        ax.scatter(x_arr, y_arr, color=color, s=55,
                   edgecolors='#333333', linewidths=0.6,
                   alpha=0.85, zorder=4, label=MODEL_LABELS[mn])

        if len(x_arr) >= 2:
            x_line = np.linspace(x_arr.min(), x_arr.max(), 80)
            fit = _regress_ci(x_arr, y_arr, x_line)
            if fit is not None:
                y_fit, y_lo, y_hi, st = fit
                ax.plot(x_line, y_fit, color=color, lw=2.0, alpha=0.80, zorder=3)
                ax.fill_between(x_line, y_lo, y_hi,
                                color=color, alpha=0.12, zorder=2)
                p_str = f'p={st["p"]:.2f}' if st['p'] >= 0.005 else f'p={st["p"]:.3f}'
                stat_lines.append(
                    f'{MODEL_LABELS[mn].split()[0]:9s}'
                    f'  n={st["n"]}'
                    f'  μ={y_arr.mean():.1f}±{y_arr.std():.1f} mg/kg'
                    f'  r={st["r"]:.2f}  {p_str}'
                )

    if stat_lines:
        ax.text(0.02, 0.97, '\n'.join(stat_lines),
                transform=ax.transAxes,
                fontsize=6.2, family='monospace', va='top',
                bbox=dict(boxstyle='round,pad=0.3', fc='white',
                          ec='#cccccc', alpha=0.88))


def plot_population_covariate(results: dict, patients: list,
                               outdir: str) -> None:
    covariates = [
        ('Age (years)',    lambda p: p['age']),
        ('Weight (kg)',    lambda p: p['weight']),
        ('Height (cm)',    lambda p: p['height']),
        ('BMI (kg m⁻²)',  lambda p: p['weight'] / (p['height'] / 100) ** 2),
    ]

    sex_variants = [
        ('both',   'Male & Female combined',  None),
        ('male',   'Male',                    0),
        ('female', 'Female',                  1),
    ]

    # Figure-level legend handles (shared across all panels)
    legend_handles = [
        Line2D([0], [0], color=MODEL_COLORS[mn], lw=2.2,
               marker='o', markersize=6, markeredgecolor='#333333',
               markeredgewidth=0.5, label=MODEL_LABELS[mn])
        for mn in ADULT_MODELS
    ]

    for sex_key, sex_label, sex_val in sex_variants:
        pats = [p for p in patients
                if sex_val is None or p['gender'] == sex_val]

        fig, axes = plt.subplots(2, 2, figsize=(12, 9),
                                 constrained_layout=True)
        n_str = f'n = {len(pats)} patients'
        fig.suptitle(
            f'Population dosing covariate analysis  ·  STANPUMP CET  ·  '
            f'Ce = {CE_TARGET} μg/mL  ·  {DURATION_MIN} min  ·  '
            f'{sex_label}  ({n_str})\n'
            'Line = linear regression  ·  Band = 95% CI  ·  '
            'Stats: model mean ± SD dose/kg  ·  r = Pearson  ·  ◆ = model colour',
            fontsize=10, fontweight='bold',
        )

        for ax, (xlabel, cov_fn) in zip(axes.flat, covariates):
            cov_title = xlabel.split('(')[0].strip()
            ax.set_title(cov_title, fontsize=9, fontweight='bold')
            _pop_covariate_panel(ax, pats, results, cov_fn, xlabel)

        fig.legend(handles=legend_handles, loc='lower center', ncol=4,
                   fontsize=9, frameon=True, bbox_to_anchor=(0.5, -0.02))

        out_path = os.path.join(outdir, f'stanpump_population_{sex_key}.png')
        plt.savefig(out_path, dpi=160, bbox_inches='tight')
        plt.close()
        print(f'Saved: {out_path}')


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 — Primary + secondary endpoint comparison (adult patients only)
# ─────────────────────────────────────────────────────────────────────────────

def plot_endpoint_comparison(results: dict, out_path: str) -> None:
    """2×4 bar chart: 8 endpoints × adult patients, one colour per model."""
    adult_pats = [p for p in PATIENTS if p['age'] >= 18]
    n_pat  = len(adult_pats)
    x      = np.arange(n_pat)
    w      = 0.22
    n_mod  = len(ADULT_MODELS)
    offsets = np.linspace(-(n_mod - 1) / 2, (n_mod - 1) / 2, n_mod) * w

    endpoints = [
        ('dose_5min_mg',          'Dose @ 5 min (mg)'),
        ('dose_30min_mg',         'Dose @ 30 min (mg)'),
        ('dose_60min_mg',         'Dose @ 60 min (mg)'),
        ('dose_end_mg',           'Total dose (mg)'),
        ('peak_rate_mlh',         'Peak rate (mL/h)'),
        ('auc_rate_mlh_min',      'AUC rate (mL·min)'),
        ('stabilization_time_min','Stab. time (min)'),
        ('infusion_cv',           'Infusion CV (maintenance)'),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(18, 8), constrained_layout=True)
    fig.suptitle(
        'Primary & secondary endpoint comparison — adult patients (≥ 18 yr)  ·  STANPUMP CET  ·  '
        f'Ce = {CE_TARGET} μg/mL  ·  {DURATION_MIN} min  ·  rate step 0.1 mL/h  ·  '
        'cyclic syringe changes [180 s, 20 s, …]',
        fontsize=10, fontweight='bold',
    )

    for ax, (col, label) in zip(axes.flat, endpoints):
        ax.set_title(label, fontsize=9, fontweight='bold', loc='left')
        for mi, mn in enumerate(ADULT_MODELS):
            vals = []
            for p in adult_pats:
                res = results[mn].get(p['id'])
                v = res.get(col) if res else None
                vals.append(float(v) if v is not None else 0.0)
            bars = ax.bar(x + offsets[mi], vals, w,
                          label=MODEL_LABELS[mn], color=MODEL_COLORS[mn],
                          edgecolor='white', linewidth=0.4, alpha=0.88, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels([p['label'] for p in adult_pats],
                           rotation=25, ha='right', fontsize=7.5)
        ax.grid(axis='y', color='#cccccc', linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=7.5)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.12)

    handles = [Patch(facecolor=MODEL_COLORS[mn], label=MODEL_LABELS[mn],
                     edgecolor='white', alpha=0.88)
               for mn in ADULT_MODELS]
    fig.legend(handles=handles, loc='lower center', ncol=3,
               fontsize=9, frameon=True, bbox_to_anchor=(0.5, -0.02))

    plt.savefig(out_path, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


def plot_secondary_endpoints(results: dict, out_path: str) -> None:
    """1×3 bar chart: secondary endpoint comparison (swaps, zero-infusion, rate-change freq)."""
    adult_pats = [p for p in PATIENTS if p['age'] >= 18]
    n_pat  = len(adult_pats)
    x      = np.arange(n_pat)
    w      = 0.22
    n_mod  = len(ADULT_MODELS)
    offsets = np.linspace(-(n_mod - 1) / 2, (n_mod - 1) / 2, n_mod) * w

    endpoints = [
        ('n_changes',          'Syringe swaps (#)'),
        ('zero_infusion_dur_s','Zero-infusion (s)'),
        ('rate_change_freq',   'Rate-change freq (/min)'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True)
    fig.suptitle(
        'Secondary endpoints — adult patients (≥ 18 yr)  ·  STANPUMP CET',
        fontsize=11, fontweight='bold',
    )

    for ax, (col, label) in zip(axes, endpoints):
        ax.set_title(label, fontsize=10, fontweight='bold', loc='left')
        for mi, mn in enumerate(ADULT_MODELS):
            vals = []
            for p in adult_pats:
                res = results[mn].get(p['id'])
                v = res.get(col) if res else None
                vals.append(float(v) if v is not None else 0.0)
            ax.bar(x + offsets[mi], vals, w,
                   label=MODEL_LABELS[mn], color=MODEL_COLORS[mn],
                   edgecolor='white', linewidth=0.4, alpha=0.88, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels([p['label'] for p in adult_pats],
                           rotation=25, ha='right', fontsize=8)
        ax.grid(axis='y', color='#cccccc', linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if col == 'n_changes':
            ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    handles = [Patch(facecolor=MODEL_COLORS[mn], label=MODEL_LABELS[mn],
                     edgecolor='white', alpha=0.88)
               for mn in ADULT_MODELS]
    fig.legend(handles=handles, loc='lower center', ncol=3,
               fontsize=9, frameon=True, bbox_to_anchor=(0.5, -0.02))

    plt.savefig(out_path, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


# ─────────────────────────────────────────────────────────────────────────────
# Figure 5 — Ce time courses + inter-model divergence (adult patients)
# ─────────────────────────────────────────────────────────────────────────────

def plot_ce_divergence(results: dict, out_path: str) -> None:
    """One panel per adult patient: Ce curves for all 3 adult models + gray divergence band."""
    adult_pats = [p for p in PATIENTS if p['age'] >= 18]
    n = len(adult_pats)

    fig, axes = plt.subplots(1, n, figsize=(3.6 * n, 4.8), constrained_layout=True)
    fig.suptitle(
        'Effect-site concentration over time — adult patients  ·  3 adult models  ·  '
        'gray band = inter-model range',
        fontsize=10, fontweight='bold',
    )

    step = 5
    for ax, pat in zip(axes, adult_pats):
        pid = pat['id']
        ax.set_title(pat['label'], fontsize=8, fontweight='bold')
        ax.axhline(CE_TARGET, color='#444', lw=0.9, ls='--', alpha=0.6, zorder=4)
        ax.axhline(WAKEUP_CE, color='#aaa', lw=0.7, ls=':', alpha=0.8, zorder=3)
        ax.axvline(DURATION_MIN, color='#bbb', lw=0.8, ls='--', zorder=3)
        ax.set_xlabel('Time (min)', fontsize=7)
        ax.set_ylabel('Ce (μg mL⁻¹)', fontsize=7)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.18)
        ax.set_xlim(0, DURATION_MIN + 5)

        ces, t_ref = [], None
        for mn in ADULT_MODELS:
            res = results[mn].get(pid)
            if res is None:
                continue
            t_full = res['time_min'][:DURATION_S + 1]
            ce_full = res['ce'][:DURATION_S + 1]
            ax.plot(t_full[::step], ce_full[::step],
                    color=MODEL_COLORS[mn], lw=1.8, alpha=0.90, zorder=5,
                    label=MODEL_LABELS[mn].split()[0])
            ces.append(ce_full)
            t_ref = t_full

        if len(ces) >= 2 and t_ref is not None:
            min_len = min(len(c) for c in ces)
            ce_arr  = np.array([c[:min_len] for c in ces])
            t_plot  = t_ref[:min_len:step]
            ax.fill_between(t_plot,
                            ce_arr[:, ::step].min(axis=0),
                            ce_arr[:, ::step].max(axis=0),
                            color='#888888', alpha=0.18, zorder=2, label='Model range')

        # Syringe changes from Marsh
        ref = results['Marsh'].get(pid)
        if ref:
            for ts, dur in ref['syringe_changes']:
                ax.axvspan(ts / 60, (ts + dur) / 60,
                           color='#cccccc', alpha=0.35, zorder=1)

        ymax = CE_TARGET * 1.8
        ax.set_ylim(0, ymax)
        ax.text(DURATION_MIN + 0.5, ymax * 0.06,
                'off', ha='left', va='bottom', fontsize=6, color='#999')
        ax.legend(fontsize=6.5, loc='upper right', framealpha=0.85)

    plt.savefig(out_path, dpi=160, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_path}')


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    outdir = os.path.join(os.path.dirname(__file__), 'figures')
    os.makedirs(outdir, exist_ok=True)

    print('Running simulations …')
    results = run_all()

    print('\nPlotting …')
    plot_summary(results, os.path.join(outdir, 'stanpump_summary.png'))

    for pat in PATIENTS:
        plot_patient(pat, results,
                     os.path.join(outdir, f'stanpump_patient_{pat["id"]}.png'))

    plot_population_covariate(results, PATIENTS, outdir)

    plot_endpoint_comparison(results, os.path.join(outdir, 'stanpump_endpoints_primary.png'))
    plot_secondary_endpoints(results, os.path.join(outdir, 'stanpump_endpoints_secondary.png'))
    plot_ce_divergence(results, os.path.join(outdir, 'stanpump_ce_divergence.png'))

    print('\nAll done.')
