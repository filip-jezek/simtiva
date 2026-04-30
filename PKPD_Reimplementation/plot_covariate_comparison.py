"""
plot_covariate_comparison.py

Runs all four models (Marsh, Paedfusor, Eleveld, EleveldUpdated) on the full
1033-patient dataset using the algebraic solver, then generates comprehensive
covariate analysis plots for each model and a combined comparison figure.

If optimized_theta.npy is found, EleveldUpdated uses the optimized parameters.
Otherwise it uses the default (Eleveld-equivalent) sigmoid initialization.

Outputs (all in figures/):
    covariate_Marsh.png
    covariate_Paedfusor.png
    covariate_Eleveld.png
    covariate_EleveldUpdated.png
    covariate_comparison_all.png      <- combined panel
"""

import os
import sys
import time
import warnings
import multiprocessing as mp
from functools import partial

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH  = os.path.join(SCRIPT_DIR, '..', 'data', 'supplementary_digital_content_1.txt')
FIG_DIR    = os.path.join(SCRIPT_DIR, 'figures')
THETA_NPY  = os.path.join(SCRIPT_DIR, 'optimized_theta.npy')
os.makedirs(FIG_DIR, exist_ok=True)

sys.path.insert(0, SCRIPT_DIR)
from core_solvers import run_algebraic_solver
from models.marsh import MarshModel
from models.paedfusor import PaedfusorModel
from models.eleveld import EleveldModel
from models.eleveld_updated import EleveldUpdatedModel, THETA_DEFAULT

# ── load optimized theta ──────────────────────────────────────────────────────
theta_opt = None
opt_label = 'default'
if os.path.isfile(THETA_NPY):
    theta_opt = np.load(THETA_NPY)
    opt_label = 'optimized'
    print(f"Using optimized theta from {THETA_NPY}")
else:
    print("No optimized_theta.npy — EleveldUpdated will use default sigmoid init.")

MODEL_CFG = {
    'Marsh':         {'color': '#e74c3c', 'marker': 'o'},
    'Paedfusor':     {'color': '#2ecc71', 'marker': 's'},
    'Eleveld':       {'color': '#9b59b6', 'marker': 'D'},
    'EleveldUpdated':{'color': '#f39c12', 'marker': '^'},
}

# ──────────────────────────────────────────────────────────────────────────────
# DATA LOADING
# ──────────────────────────────────────────────────────────────────────────────

def load_data(max_patients=None):
    df = pd.read_csv(DATA_PATH, sep=r'\s+')
    pids = df['ID'].unique()
    if max_patients:
        pids = pids[:max_patients]

    patients = []
    for pid in pids:
        p = df[df['ID'] == pid].sort_values('TIME')
        first = p.iloc[0]

        weight = first['WGT']; age = first['AGE']; height = first['HGT']
        m1f2   = int(first['M1F2'])
        gender = 0 if m1f2 == 1 else 1
        opioid = (first['TECH'] == 2)
        stdy   = first.get('STDY', np.nan)
        pma    = first.get('PMA', age * 52.1429 + 40)

        dosing = p[p['EVID'] == 1]
        infusion_bounds = []
        for _, row in dosing.iterrows():
            t = row['TIME']
            if row['RATE'] > 0:
                dur = row['AMT'] / row['RATE']
                infusion_bounds.append((t, t + dur, row['RATE']))
            elif row['AMT'] > 0:
                infusion_bounds.append((t, t + 1.0, row['AMT']))

        obs = p[(p['EVID'] == 0) & (p['DV'] > 0)]
        if len(obs) == 0:
            continue

        t_obs = obs['TIME'].values.astype(float)
        c_obs = np.exp(obs['DV'].values.astype(float))
        max_t = t_obs.max()
        t_eval = np.linspace(0, max_t + 10, max(int(max_t * 10), 100))
        obs_idx = np.clip(np.searchsorted(t_eval, t_obs), 0, len(t_eval) - 1)

        patients.append({
            'pid': int(pid), 'weight': weight, 'age': age, 'height': height,
            'gender': gender, 'opioid': opioid, 'stdy': stdy, 'pma': pma,
            'm1f2': m1f2,
            'infusion_bounds': infusion_bounds,
            't_obs': t_obs, 'c_obs': c_obs, 't_eval': t_eval, 'obs_idx': obs_idx,
        })
    return patients


def make_infusion_func(bounds):
    def f(t):
        r = 0.0
        for s, e, rate in bounds:
            if s <= t < e:
                r += rate
        return r
    return f


# ──────────────────────────────────────────────────────────────────────────────
# SINGLE-PATIENT EVALUATION
# ──────────────────────────────────────────────────────────────────────────────

def evaluate_patient(pat, theta_updated):
    results = {}
    infusion_func = make_infusion_func(pat['infusion_bounds'])
    w, a, h, g, o = pat['weight'], pat['age'], pat['height'], pat['gender'], pat['opioid']

    model_instances = {
        'Marsh':     MarshModel(weight=w, age=a, height=h, gender=g),
        'Paedfusor': PaedfusorModel(weight=w, age=a, height=h, gender=g),
        'Eleveld':   EleveldModel(weight=w, age=a, height=h, gender=g, opioid=o),
        'EleveldUpdated': EleveldUpdatedModel(weight=w, age=a, height=h, gender=g,
                                               opioid=o, theta=theta_updated),
    }

    for name, model in model_instances.items():
        try:
            params = model.get_parameters()
            if params.get('vc', 0) <= 0:
                results[name] = np.nan
                continue
            cp, _ = run_algebraic_solver(params, infusion_func, pat['t_eval'])
            cp_obs = cp[pat['obs_idx']]
            valid = cp_obs > 0
            if valid.sum() == 0:
                results[name] = np.nan
                continue
            rmse = np.sqrt(np.mean((pat['c_obs'][valid] - cp_obs[valid]) ** 2))
            results[name] = rmse
        except Exception:
            results[name] = np.nan

    return {
        'pid': pat['pid'], 'age': pat['age'], 'weight': pat['weight'],
        'height': pat['height'], 'gender': pat['gender'], 'opioid': pat['opioid'],
        'stdy': pat['stdy'],
        **results
    }


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE HELPERS
# ──────────────────────────────────────────────────────────────────────────────

DARK_BG   = '#ffffff'
PANEL_BG  = '#ffffff'
GRID_CLR  = '#e0e0e0'

def style_ax(ax):
    ax.set_facecolor(PANEL_BG)
    for spine in ax.spines.values():
        spine.set_edgecolor('#cccccc')
    ax.tick_params(colors='black', labelsize=8)
    ax.xaxis.label.set_color('black')
    ax.yaxis.label.set_color('black')
    ax.title.set_color('black')
    ax.grid(True, alpha=0.3, color='gray', linewidth=0.5)


def make_single_model_figure(df_m, model_name, color, out_path):
    """Per-model covariate figure: 3×3 grid identical to original analyze_covariates.py."""
    fig = plt.figure(figsize=(20, 14))
    fig.patch.set_facecolor(DARK_BG)
    plt.rcParams.update({'text.color': 'black'})

    fig.suptitle(
        f'{model_name} — RMSE vs Patient Covariates  (N={len(df_m)})',
        fontsize=15, fontweight='bold', color='black', y=0.99
    )
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.48, wspace=0.36)

    df_m = df_m.copy()
    df_m['Sex']      = df_m['gender'].map({0: 'Male', 1: 'Female'})
    df_m['BMI']      = df_m['weight'] / (df_m['height'] / 100) ** 2
    df_m['AgeGroup'] = pd.cut(df_m['age'], bins=[0, 35, 60, 200],
                               labels=['Young\n(<35)', 'Middle\n(35-60)', 'Senior\n(>60)'])

    # ── scatter helper ────────────────────────────────────────────────────────
    def scatter_lowess(ax, xcol, xlabel):
        style_ax(ax)
        sex_style = [('Male', 'o', '#5b9bd5'), ('Female', '^', '#e36b6b')]
        for sex, mk, clr in sex_style:
            sub = df_m[df_m['Sex'] == sex].dropna(subset=[xcol, 'RMSE'])
            ax.scatter(sub[xcol], sub['RMSE'], alpha=0.20, s=12, color=clr, marker=mk, label=sex)
            if len(sub) > 10:
                sm = lowess(sub['RMSE'], sub[xcol], frac=0.4, return_sorted=True)
                ax.plot(sm[:, 0], sm[:, 1], color=clr, linewidth=2.5)
        valid = df_m.dropna(subset=[xcol, 'RMSE'])
        rho, pval = stats.spearmanr(valid[xcol], valid['RMSE'])
        ax.set_xlabel(xlabel); ax.set_ylabel('RMSE (mg/L)')
        ax.set_title(f'vs {xlabel}  (ρ={rho:+.2f}, p={pval:.3f})')
        ax.legend(fontsize=8, labelcolor='black')
        ax.set_ylim(0, 15)

    scatter_lowess(fig.add_subplot(gs[0, 0]), 'age', 'Age (years)')
    scatter_lowess(fig.add_subplot(gs[0, 1]), 'weight', 'Weight (kg)')
    scatter_lowess(fig.add_subplot(gs[0, 2]), 'BMI', 'BMI (kg/m²)')

    # ── heatmap Age × Weight ──────────────────────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 0]); style_ax(ax4)
    df_m['AgeBin'] = pd.cut(df_m['age'],    np.arange(0, 101, 10))
    df_m['WgtBin'] = pd.cut(df_m['weight'], np.arange(30, 151, 10))
    heat = df_m.groupby(['AgeBin', 'WgtBin'], observed=True)['RMSE'].median().unstack('WgtBin')
    im = ax4.imshow(heat.values, aspect='auto', origin='lower', cmap='RdYlGn_r', vmin=0)
    plt.colorbar(im, ax=ax4, label='Median RMSE (mg/L)')
    ax4.set_xticks(range(len(heat.columns)))
    ax4.set_xticklabels([f'{c.mid:.0f}' for c in heat.columns], rotation=45, fontsize=7)
    ax4.set_yticks(range(len(heat.index)))
    ax4.set_yticklabels([f'{int(c.mid)}' for c in heat.index], fontsize=7)
    ax4.set_xlabel('Weight (kg)'); ax4.set_ylabel('Age (years)')
    ax4.set_title('Median RMSE: Age × Weight')

    # ── box: age group × sex ──────────────────────────────────────────────────
    ax5 = fig.add_subplot(gs[1, 1]); style_ax(ax5)
    age_grps = ['Young\n(<35)', 'Middle\n(35-60)', 'Senior\n(>60)']
    sex_styles = [('Male', '#5b9bd5'), ('Female', '#e36b6b')]
    pos, labs, data_boxes, clrs = [], [], [], []
    base = 0
    for ag in age_grps:
        for sex, clr in sex_styles:
            sub = df_m[(df_m['AgeGroup'] == ag) & (df_m['Sex'] == sex)]['RMSE'].dropna()
            pos.append(base); labs.append(f'{sex[0]}\n{ag}')
            data_boxes.append(sub.values); clrs.append(clr); base += 0.9
        base += 1.0
    bp = ax5.boxplot(data_boxes, positions=pos, widths=0.75, patch_artist=True, notch=False,
                     medianprops=dict(color='black', linewidth=2),
                     whiskerprops=dict(color='#888'), capprops=dict(color='#888'),
                     flierprops=dict(marker='.', alpha=0.3, markersize=4))
    for patch, clr in zip(bp['boxes'], clrs):
        patch.set_facecolor(clr); patch.set_alpha(0.7)
    ax5.set_xticks(pos); ax5.set_xticklabels(labs, fontsize=6.5)
    ax5.set_ylabel('RMSE (mg/L)'); ax5.set_title('RMSE by Sex × Age Group')
    ax5.set_ylim(0, 7.5)

    # ── box: opioid ───────────────────────────────────────────────────────────
    ax6 = fig.add_subplot(gs[1, 2]); style_ax(ax6)
    opi_no  = df_m[df_m['opioid'] == False]['RMSE'].dropna()
    opi_yes = df_m[df_m['opioid'] == True]['RMSE'].dropna()
    bp6 = ax6.boxplot([opi_no, opi_yes], patch_artist=True, labels=['No opioid', 'With opioid'],
                      medianprops=dict(color='black', linewidth=2))
    bp6['boxes'][0].set_facecolor('#5b9bd5'); bp6['boxes'][0].set_alpha(0.8)
    bp6['boxes'][1].set_facecolor('#e36b6b'); bp6['boxes'][1].set_alpha(0.8)
    _, p6 = stats.mannwhitneyu(opi_no, opi_yes)
    ax6.set_ylabel('RMSE (mg/L)'); ax6.set_title(f'Opioid co-admin  (MW p={p6:.3f})')
    ax6.set_ylim(0, 7.5)

    # ── correlation bars ──────────────────────────────────────────────────────
    ax7 = fig.add_subplot(gs[2, :2]); style_ax(ax7)
    cov_cols = [('age', 'Age'), ('weight', 'Weight'), ('height', 'Height'),
                ('BMI', 'BMI'), ('opioid', 'Opioid'), ('stdy', 'Study Site')]
    rhos, pvs, labs7 = [], [], []
    for col, label in cov_cols:
        sub7 = df_m[[col, 'RMSE']].dropna().astype(float)
        if len(sub7) < 5:
            continue
        r, p = stats.spearmanr(sub7[col], sub7['RMSE'])
        rhos.append(r); pvs.append(p); labs7.append(label)
    order = np.argsort(np.abs(rhos))[::-1]
    rhos_s = [rhos[i] for i in order]
    pvs_s  = [pvs[i]  for i in order]
    labs_s = [labs7[i] for i in order]
    bar_clrs = ['#e74c3c' if r > 0 else '#3498db' for r in rhos_s]
    ax7.barh(range(len(rhos_s)), rhos_s, color=bar_clrs, alpha=0.85)
    ax7.set_yticks(range(len(labs_s))); ax7.set_yticklabels(labs_s)
    ax7.axvline(0, color='black', linewidth=0.8)
    for i, (r, p) in enumerate(zip(rhos_s, pvs_s)):
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        ax7.text(r + (0.005 if r >= 0 else -0.005), i, f'{r:+.3f}{sig}',
                 va='center', ha='left' if r >= 0 else 'right', fontsize=9, color='black')
    ax7.set_xlabel('Spearman ρ (RMSE)'); ax7.set_title('Covariate Correlation with RMSE (ranked)')

    # ── study site box ────────────────────────────────────────────────────────
    ax8 = fig.add_subplot(gs[2, 2]); style_ax(ax8)
    sites = sorted(df_m['stdy'].dropna().unique())
    site_data = [df_m[df_m['stdy'] == s]['RMSE'].dropna().values for s in sites]
    if site_data:
        bp8 = ax8.boxplot(site_data, patch_artist=True, labels=[f'S{int(s)}' for s in sites],
                          medianprops=dict(color='black', linewidth=1.5))
        for patch in bp8['boxes']:
            patch.set_facecolor(color); patch.set_alpha(0.65)
    ax8.set_ylabel('RMSE (mg/L)'); ax8.set_title('RMSE by Study Site')
    ax8.tick_params(axis='x', labelsize=7)

    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor=DARK_BG)
    plt.close()
    print(f"  Saved: {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# COMBINED COMPARISON FIGURE
# ──────────────────────────────────────────────────────────────────────────────

def make_comparison_figure(results_df, out_path):
    """
    4-panel summary: RMSE boxplot, age scatter, weight scatter, opioid stratification.
    All 4 models overlaid for direct comparison.
    """
    models = list(MODEL_CFG.keys())
    fig = plt.figure(figsize=(20, 14))
    fig.patch.set_facecolor(DARK_BG)
    plt.rcParams.update({'text.color': 'black'})
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.35)
    gs.update(left=0.06, right=0.98, top=0.91, bottom=0.08)

    # ── Panel A: RMSE box plot ────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0, :])
    style_ax(ax1)
    data = [results_df[m].dropna().values for m in models]
    medians = [np.median(d) for d in data]
    bp = ax1.boxplot(
        data, positions=range(len(models)), widths=0.55,
        patch_artist=True, notch=False, vert=True,
        medianprops=dict(color='black', linewidth=2.5),
        whiskerprops=dict(color='#888', linewidth=1.5),
        capprops=dict(color='#888'), flierprops=dict(marker='.', alpha=0.2, markersize=3)
    )
    for patch, m in zip(bp['boxes'], models):
        patch.set_facecolor(MODEL_CFG[m]['color']); patch.set_alpha(0.75)
    ax1.set_xticks(range(len(models)))
    ax1.set_xticklabels(models, fontsize=11, color='black')
    ax1.set_ylabel('RMSE (mg/L)', fontsize=11)
    ax1.set_title('Population-Level RMSE Distribution — All Models', fontsize=13, pad=8)
    for i, (med, d) in enumerate(zip(medians, data)):
        ax1.text(i, med + 0.03, f'{med:.3f}', ha='center', va='bottom',
                 fontsize=10, fontweight='bold', color='black')
    ax1.set_ylim(0, 15)

    # ── Panel B: RMSE vs Age (LOWESS, all models) ─────────────────────────────
    ax2 = fig.add_subplot(gs[1, 0]); style_ax(ax2)
    for m in models:
        sub = results_df[['age', m]].dropna()
        sm = lowess(sub[m], sub['age'], frac=0.35, return_sorted=True)
        ax2.plot(sm[:, 0], sm[:, 1], color=MODEL_CFG[m]['color'], lw=2.5, label=m)
    ax2.set_xlabel('Age (years)'); ax2.set_ylabel('RMSE (mg/L)')
    ax2.set_title('LOWESS RMSE vs Age'); ax2.legend(fontsize=8, labelcolor='black')
    ax2.set_ylim(0, 15)

    # ── Panel C: RMSE vs Weight ───────────────────────────────────────────────
    ax3 = fig.add_subplot(gs[1, 1]); style_ax(ax3)
    for m in models:
        sub = results_df[['weight', m]].dropna()
        sm = lowess(sub[m], sub['weight'], frac=0.45, return_sorted=True)
        ax3.plot(sm[:, 0], sm[:, 1], color=MODEL_CFG[m]['color'], lw=2.5, label=m)
    ax3.set_xlabel('Weight (kg)'); ax3.set_ylabel('RMSE (mg/L)')
    ax3.set_title('LOWESS RMSE vs Weight'); ax3.legend(fontsize=8, labelcolor='black')
    ax3.set_ylim(0, 15)

    # ── Panel D: Opioid stratification ───────────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 2]); style_ax(ax4)
    x_pos = np.arange(len(models))
    for offset, opioid, ls in [(-0.2, False, '--'), (0.2, True, '-')]:
        sub = results_df[results_df['opioid'] == opioid]
        means = [sub[m].mean() for m in models]
        stds  = [sub[m].std()  for m in models]
        label = 'With opioid' if opioid else 'No opioid'
        ax4.bar(x_pos + offset, means, 0.35, yerr=stds, capsize=4,
                color=[MODEL_CFG[m]['color'] for m in models], alpha=0.75 if opioid else 0.45,
                label=label, error_kw=dict(ecolor='black', elinewidth=1))
    ax4.set_xticks(x_pos); ax4.set_xticklabels(models, rotation=12, fontsize=8)
    ax4.set_ylabel('Mean RMSE ± SD (mg/L)'); ax4.set_title('RMSE by Opioid Co-administration')
    ax4.legend(fontsize=9, labelcolor='black')
    ax4.set_ylim(0, 7.5)

    fig.suptitle(
        'Model Comparison: Marsh | Paedfusor | Eleveld | EleveldUpdated\n'
        f'N={"~1029"} patients, 11,204 observations — EleveldUpdated theta: {opt_label}',
        fontsize=13, fontweight='bold', color='black', y=0.97
    )
    plt.savefig(out_path, dpi=160, bbox_inches='tight', facecolor=DARK_BG)
    plt.close()
    print(f"  Saved: {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    print("Loading patient data...")
    t0 = time.time()
    patients = load_data()
    print(f"  {len(patients)} patients loaded in {time.time()-t0:.1f}s")

    theta_updated = list(theta_opt) if theta_opt is not None else list(THETA_DEFAULT)

    print("Evaluating all models (multiprocessing)...")
    t0 = time.time()
    n_cores = min(mp.cpu_count(), 8)
    eval_fn = partial(evaluate_patient, theta_updated=theta_updated)
    with mp.Pool(processes=n_cores) as pool:
        records = pool.map(eval_fn, patients)
    print(f"  Done in {time.time()-t0:.1f}s")

    results_df = pd.DataFrame(records)
    results_df['RMSE'] = results_df['EleveldUpdated']  # default RMSE column for single-model figs

    # Print summary table
    print("\n{:<20} {:>10} {:>10} {:>10} {:>10}".format(
        'Model', 'Mean', 'Median', 'P90', 'Max'))
    print('-' * 62)
    for m in MODEL_CFG.keys():
        vals = results_df[m].dropna()
        print("{:<20} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f}".format(
            m, vals.mean(), vals.median(), vals.quantile(0.90), vals.max()))

    # Per-model figures
    for model_name, cfg in MODEL_CFG.items():
        df_m = results_df.copy()
        df_m['RMSE'] = df_m[model_name]
        out = os.path.join(FIG_DIR, f'covariate_{model_name}.png')
        print(f"\nGenerating {model_name} covariate figure...")
        make_single_model_figure(df_m, model_name, cfg['color'], out)

    # Combined comparison
    print("\nGenerating combined comparison figure...")
    make_comparison_figure(results_df, os.path.join(FIG_DIR, 'covariate_comparison_all.png'))

    print("\nAll figures written to:", FIG_DIR)
    return results_df


if __name__ == '__main__':
    results = main()
