"""
plot_pk_metrics.py

Evaluates Marsh, Paedfusor, Eleveld, and EleveldUpdated using standard PK metrics:
MDPE (Median Directional Performance Error) - measures BIAS (overshoot/undershoot)
MDAPE (Median Absolute Performance Error) - measures INACCURACY

Outputs standard PK visualizations:
    figures/mdape_comparison_all.png
    figures/mdpe_bias_analysis.png
    figures/mdape_Marsh.png ...
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
from models.schnider import SchniderModel
from models.marsh import MarshModel
from models.paedfusor import PaedfusorModel
from models.eleveld import EleveldModel
from models.eleveld_updated import EleveldUpdatedModel, THETA_DEFAULT

theta_opt = None
opt_label = 'default'
if os.path.isfile(THETA_NPY):
    theta_opt = np.load(THETA_NPY)
    opt_label = 'optimized'

MODEL_CFG = {
    'Schnider':      {'color': '#3498db', 'marker': 'v'},
    'Marsh':         {'color': '#e74c3c', 'marker': 'o'},
    'Paedfusor':     {'color': '#2ecc71', 'marker': 's'},
    'Eleveld':       {'color': '#9b59b6', 'marker': 'D'},
    'EleveldUpdated':{'color': '#f39c12', 'marker': '^'},
}

# ──────────────────────────────────────────────────────────────────────────────
# DATA LOADING (Identical to before)
# ──────────────────────────────────────────────────────────────────────────────

def load_data(max_patients=None):
    df = pd.read_csv(DATA_PATH, sep=r'\s+')
    pids = df['ID'].unique()
    if max_patients: pids = pids[:max_patients]

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
                infusion_bounds.append((t, t + row['AMT']/row['RATE'], row['RATE']))
            elif row['AMT'] > 0:
                infusion_bounds.append((t, t + 1.0, row['AMT']))

        obs = p[(p['EVID'] == 0) & (p['DV'] > 0)]
        if len(obs) == 0: continue

        t_obs = obs['TIME'].values.astype(float)
        c_obs = np.exp(obs['DV'].values.astype(float))
        max_t = t_obs.max()
        t_eval = np.linspace(0, max_t + 10, max(int(max_t * 10), 100))
        obs_idx = np.clip(np.searchsorted(t_eval, t_obs), 0, len(t_eval) - 1)

        patients.append({
            'pid': int(pid), 'weight': weight, 'age': age, 'height': height,
            'gender': gender, 'opioid': opioid, 'stdy': stdy, 'pma': pma,
            'm1f2': m1f2, 'infusion_bounds': infusion_bounds,
            't_obs': t_obs, 'c_obs': c_obs, 't_eval': t_eval, 'obs_idx': obs_idx,
        })
    return patients

def make_infusion_func(bounds):
    def f(t):
        r = 0.0
        for s, e, rate in bounds:
            if s <= t < e: r += rate
        return r
    return f

# ──────────────────────────────────────────────────────────────────────────────
# EVALUATION (Now computing MDAPE and MDPE)
# ──────────────────────────────────────────────────────────────────────────────

def evaluate_patient(pat, theta_updated):
    results = {}
    inf_f = make_infusion_func(pat['infusion_bounds'])
    w, a, h, g, o = pat['weight'], pat['age'], pat['height'], pat['gender'], pat['opioid']

    model_instances = {
        'Schnider':  SchniderModel(weight=w, age=a, height=h, gender=g),
        'Marsh':     MarshModel(weight=w, age=a, height=h, gender=g),
        'Paedfusor': PaedfusorModel(weight=w, age=a, height=h, gender=g),
        'Eleveld':   EleveldModel(weight=w, age=a, height=h, gender=g, opioid=o),
        'EleveldUpdated': EleveldUpdatedModel(weight=w, age=a, height=h, gender=g, opioid=o, theta=theta_updated),
    }

    for name, model in model_instances.items():
        try:
            params = model.get_parameters()
            if params.get('vc', 0) <= 0:
                raise ValueError
            cp, _ = run_algebraic_solver(params, inf_f, pat['t_eval'])
            cp_obs = cp[pat['obs_idx']]
            valid = cp_obs > 0
            if valid.sum() == 0: raise ValueError
            
            # Classical PE formula: (Measured - Predicted) / Predicted * 100
            # Positive PE -> Model Underpredicts (Measured > Predicted)
            # Negative PE -> Model Overpredicts (Measured < Predicted)
            pe = (pat['c_obs'][valid] - cp_obs[valid]) / cp_obs[valid] * 100.0
            
            results[f'{name}_MDAPE'] = np.median(np.abs(pe))
            results[f'{name}_MDPE']  = np.median(pe)
        except Exception:
            results[f'{name}_MDAPE'] = np.nan
            results[f'{name}_MDPE']  = np.nan

    return {
        'pid': pat['pid'], 'age': pat['age'], 'weight': pat['weight'],
        'height': pat['height'], 'gender': pat['gender'], 'opioid': pat['opioid'],
        'stdy': pat['stdy'],
        **results
    }

# ──────────────────────────────────────────────────────────────────────────────
# STYLING
# ──────────────────────────────────────────────────────────────────────────────

def style_ax(ax):
    ax.set_facecolor('#ffffff')
    for spine in ax.spines.values(): spine.set_edgecolor('#cccccc')
    ax.tick_params(colors='black', labelsize=8)
    ax.xaxis.label.set_color('black')
    ax.yaxis.label.set_color('black')
    ax.title.set_color('black')
    ax.grid(True, alpha=0.3, color='gray', linewidth=0.5)

# ──────────────────────────────────────────────────────────────────────────────
# FIGURE: MDAPE (Inaccuracy) Single Model
# ──────────────────────────────────────────────────────────────────────────────

def make_single_model_figure(df, model_name, color, out_path):
    fig = plt.figure(figsize=(20, 14)); fig.patch.set_facecolor('#ffffff')
    plt.rcParams.update({'text.color': 'black'})
    fig.suptitle(f'{model_name} — MDAPE (Inaccuracy) vs Patient Covariates',
                 fontsize=15, fontweight='bold', color='black', y=0.99)
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.48, wspace=0.36)

    mdape_col = f'{model_name}_MDAPE'
    df_m = df.copy()
    df_m['Sex']      = df_m['gender'].map({0: 'Male', 1: 'Female'})
    df_m['BMI']      = df_m['weight'] / (df_m['height'] / 100) ** 2
    df_m['AgeGroup'] = pd.cut(df_m['age'], bins=[0, 35, 60, 200],
                               labels=['Young\n(<35)', 'Middle\n(35-60)', 'Senior\n(>60)'])

    def scatter_lowess(ax, xcol, xlabel):
        style_ax(ax)
        for sex, mk, clr in [('Male', 'o', '#5b9bd5'), ('Female', '^', '#e36b6b')]:
            sub = df_m[df_m['Sex'] == sex].dropna(subset=[xcol, mdape_col])
            ax.scatter(sub[xcol], sub[mdape_col], alpha=0.25, s=12, color=clr, marker=mk, label=sex)
            if len(sub) > 10:
                sm = lowess(sub[mdape_col], sub[xcol], frac=0.4, return_sorted=True)
                ax.plot(sm[:, 0], sm[:, 1], color=clr, linewidth=2.5)
        valid = df_m.dropna(subset=[xcol, mdape_col])
        if len(valid) > 5:
            rho, pval = stats.spearmanr(valid[xcol], valid[mdape_col])
            ax.set_title(f'vs {xlabel}  (ρ={rho:+.2f}, p={pval:.3f})')
        ax.set_xlabel(xlabel); ax.set_ylabel('MDAPE (%)')
        ax.legend(fontsize=8, labelcolor='black'); ax.set_ylim(0, 100)

    scatter_lowess(fig.add_subplot(gs[0, 0]), 'age', 'Age (years)')
    scatter_lowess(fig.add_subplot(gs[0, 1]), 'weight', 'Weight (kg)')
    scatter_lowess(fig.add_subplot(gs[0, 2]), 'BMI', 'BMI (kg/m²)')

    # Heatmap
    ax4 = fig.add_subplot(gs[1, 0]); style_ax(ax4)
    df_m['AgeBin'] = pd.cut(df_m['age'], np.arange(0, 101, 10))
    df_m['WgtBin'] = pd.cut(df_m['weight'], np.arange(30, 151, 10))
    heat = df_m.groupby(['AgeBin', 'WgtBin'], observed=True)[mdape_col].median().unstack('WgtBin')
    im = ax4.imshow(heat.values, aspect='auto', origin='lower', cmap='RdYlGn_r', vmin=0, vmax=60)
    plt.colorbar(im, ax=ax4, label='Median MDAPE (%)')
    ax4.set_xticks(range(len(heat.columns))); ax4.set_xticklabels([f'{c.mid:.0f}' for c in heat.columns], rotation=45, fontsize=7)
    ax4.set_yticks(range(len(heat.index))); ax4.set_yticklabels([f'{int(c.mid)}' for c in heat.index], fontsize=7)
    ax4.set_xlabel('Weight (kg)'); ax4.set_ylabel('Age (years)'); ax4.set_title('Median MDAPE: Age × Weight')

    # Box: age group × sex
    ax5 = fig.add_subplot(gs[1, 1]); style_ax(ax5)
    pos, labs, data_boxes, clrs = [], [], [], []
    base = 0
    for ag in ['Young\n(<35)', 'Middle\n(35-60)', 'Senior\n(>60)']:
        for sex, clr in [('Male', '#5b9bd5'), ('Female', '#e36b6b')]:
            sub = df_m[(df_m['AgeGroup'] == ag) & (df_m['Sex'] == sex)][mdape_col].dropna()
            pos.append(base); labs.append(f'{sex[0]}\n{ag}')
            data_boxes.append(sub.values); clrs.append(clr); base += 0.9
        base += 1.0
    bp = ax5.boxplot(data_boxes, positions=pos, widths=0.75, patch_artist=True, notch=False,
                     medianprops=dict(color='black', linewidth=2),
                     flierprops=dict(marker='.', alpha=0.3, markersize=4))
    for patch, clr in zip(bp['boxes'], clrs): patch.set_facecolor(clr); patch.set_alpha(0.7)
    ax5.set_xticks(pos); ax5.set_xticklabels(labs, fontsize=6.5)
    ax5.set_ylabel('MDAPE (%)'); ax5.set_title('MDAPE by Sex × Age Group'); ax5.set_ylim(0, 80)

    # Box: opioid
    ax6 = fig.add_subplot(gs[1, 2]); style_ax(ax6)
    o0 = df_m[df_m['opioid'] == False][mdape_col].dropna()
    o1 = df_m[df_m['opioid'] == True][mdape_col].dropna()
    bp6 = ax6.boxplot([o0, o1], patch_artist=True, labels=['No opioid', 'With opioid'], medianprops=dict(color='black', linewidth=2))
    bp6['boxes'][0].set_facecolor('#5b9bd5'); bp6['boxes'][0].set_alpha(0.8)
    bp6['boxes'][1].set_facecolor('#e36b6b'); bp6['boxes'][1].set_alpha(0.8)
    _, p6 = stats.mannwhitneyu(o0, o1) if len(o0) and len(o1) else (0, 1)
    ax6.set_ylabel('MDAPE (%)'); ax6.set_title(f'Opioid co-admin (MW p={p6:.3f})'); ax6.set_ylim(0, 80)

    plt.savefig(out_path, dpi=150, bbox_inches='tight', facecolor='#ffffff')
    plt.close()

# ──────────────────────────────────────────────────────────────────────────────
# FIGURE: MDAPE (Inaccuracy) All Models Comparison
# ──────────────────────────────────────────────────────────────────────────────

def make_mdape_comparison(df, out_path):
    models = list(MODEL_CFG.keys())
    fig = plt.figure(figsize=(20, 14)); fig.patch.set_facecolor('#ffffff')
    plt.rcParams.update({'text.color': 'black'})
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.35)
    gs.update(left=0.06, right=0.98, top=0.91, bottom=0.08)

    ax1 = fig.add_subplot(gs[0, :]); style_ax(ax1)
    data = [df[f'{m}_MDAPE'].dropna().values for m in models]
    medians = [np.median(d) for d in data]
    bp = ax1.boxplot(data, positions=range(len(models)), widths=0.55, patch_artist=True,
                     medianprops=dict(color='black', linewidth=2.5), flierprops=dict(marker='.', alpha=0.2, markersize=3))
    for patch, m in zip(bp['boxes'], models): patch.set_facecolor(MODEL_CFG[m]['color']); patch.set_alpha(0.75)
    ax1.set_xticks(range(len(models))); ax1.set_xticklabels(models, fontsize=11, color='black')
    ax1.set_ylabel('MDAPE (%)', fontsize=11); ax1.set_title('Global Inaccuracy (MDAPE) — All Models', fontsize=13)
    for i, med in enumerate(medians):
        ax1.text(i, med + 1.0, f'{med:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold', color='black')
    ax1.set_ylim(0, 100)

    for i, (col, label) in enumerate([('age', 'Age (years)'), ('weight', 'Weight (kg)')]):
        ax = fig.add_subplot(gs[1, i]); style_ax(ax)
        for m in models:
            sub = df[[col, f'{m}_MDAPE']].dropna()
            sm = lowess(sub[f'{m}_MDAPE'], sub[col], frac=0.4, return_sorted=True)
            ax.plot(sm[:, 0], sm[:, 1], color=MODEL_CFG[m]['color'], lw=2.5, label=m)
        ax.set_xlabel(label); ax.set_ylabel('MDAPE (%)'); ax.set_title(f'LOWESS MDAPE vs {label}')
        ax.legend(fontsize=8, labelcolor='black'); ax.set_ylim(0, 80)

    ax4 = fig.add_subplot(gs[1, 2]); style_ax(ax4)
    x_pos = np.arange(len(models))
    for offset, opioid, ls in [(-0.2, False, '--'), (0.2, True, '-')]:
        sub = df[df['opioid'] == opioid]
        means = [sub[f'{m}_MDAPE'].mean() for m in models]
        stds  = [sub[f'{m}_MDAPE'].std() for m in models]
        ax4.bar(x_pos + offset, means, 0.35, yerr=stds, capsize=4,
                color=[MODEL_CFG[m]['color'] for m in models], alpha=0.75 if opioid else 0.45,
                label='With opioid' if opioid else 'No opioid', error_kw=dict(ecolor='black', elinewidth=1))
    ax4.set_xticks(x_pos); ax4.set_xticklabels(models, rotation=12, fontsize=8)
    ax4.set_ylabel('Mean MDAPE ± SD (%)'); ax4.set_title('MDAPE by Opioid Co-administration'); ax4.legend(fontsize=9, labelcolor='black')
    ax4.set_ylim(0, 80)

    fig.suptitle('Model Inaccuracy Comparison (MDAPE)\n'
                 'Lower MDAPE = Better Accuracy. Proportional Error evaluated across 1,033 patients.',
                 fontsize=13, fontweight='bold', color='black', y=0.97)
    plt.savefig(out_path, dpi=160, bbox_inches='tight', facecolor='#ffffff')
    plt.close()


# ──────────────────────────────────────────────────────────────────────────────
# FIGURE: MDPE (Bias / Directionality)
# ──────────────────────────────────────────────────────────────────────────────

def make_mdpe_bias_figure(df, out_path):
    models = list(MODEL_CFG.keys())
    fig = plt.figure(figsize=(22, 14)); fig.patch.set_facecolor('#ffffff')
    plt.rcParams.update({'text.color': 'black'})
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.25)
    gs.update(left=0.05, right=0.97, top=0.91, bottom=0.08)

    # MDPE By Global
    ax1 = fig.add_subplot(gs[0, 0]); style_ax(ax1)
    data = [df[f'{m}_MDPE'].dropna().values for m in models]
    medians = [np.median(d) for d in data]
    bp = ax1.boxplot(data, positions=range(len(models)), widths=0.55, patch_artist=True,
                     medianprops=dict(color='black', linewidth=2.5), flierprops=dict(marker='.', alpha=0.2, markersize=3))
    for patch, m in zip(bp['boxes'], models): patch.set_facecolor(MODEL_CFG[m]['color']); patch.set_alpha(0.75)
    ax1.axhline(0, color='red', linewidth=2, linestyle='--', zorder=0, label='Zero Bias (Perfect)')
    ax1.set_xticks(range(len(models))); ax1.set_xticklabels(models, rotation=30, fontsize=10, color='black')
    ax1.set_ylabel('MDPE (%)', fontsize=12); ax1.set_title('Global Bias (MDPE)', fontsize=14)
    ax1.set_ylim(-100, 150)
    for i, med in enumerate(medians):
        ax1.text(i, med + (2 if med > 0 else -3), f'{med:+.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    textstr = "MDPE > 0%: Under-predicts\nMDPE < 0%: Over-predicts"
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='#ccc')
    ax1.text(0.015, 0.95, textstr, transform=ax1.transAxes, fontsize=10, verticalalignment='top', bbox=props)

    # MDPE By Sex
    ax_sex = fig.add_subplot(gs[0, 1]); style_ax(ax_sex)
    df_sex = df.copy()
    df_sex['Sex'] = df_sex['gender'].map({0: 'Male', 1: 'Female'})
    pos, labs, data_boxes, clrs = [], [], [], []
    base = 0
    for sex, clr in [('Male', '#5b9bd5'), ('Female', '#e36b6b')]:
        for m in models:
            sub = df_sex[df_sex['Sex'] == sex][f'{m}_MDPE'].dropna()
            pos.append(base); labs.append(m)
            data_boxes.append(sub.values); clrs.append(clr); base += 0.9
        base += 1.0
    bp_sex = ax_sex.boxplot(data_boxes, positions=pos, widths=0.75, patch_artist=True, notch=False,
                            medianprops=dict(color='black', linewidth=1.5), flierprops=dict(marker='.', alpha=0.2, markersize=3))
    for patch, clr in zip(bp_sex['boxes'], clrs): patch.set_facecolor(clr); patch.set_alpha(0.7)
    ax_sex.axhline(0, color='red', linewidth=1.5, linestyle='--', zorder=0)
    ax_sex.set_xticks(pos); ax_sex.set_xticklabels(labs, rotation=45, fontsize=8)
    ax_sex.set_ylabel('MDPE (%)'); ax_sex.set_title('Bias By Sex')
    ax_sex.set_ylim(-100, 150)
    ax_sex.text(0, -85, "Blue=Male, Red=Female", color='black', fontsize=9, fontweight='bold')

    # MDPE By Opioid
    ax_op = fig.add_subplot(gs[0, 2]); style_ax(ax_op)
    df_op = df.copy()
    df_op['OpioidLabel'] = df_op['opioid'].map({False: 'No Opioid', True: 'With Opioid'})
    pos, labs, data_boxes, clrs = [], [], [], []
    base = 0
    for op_lbl, clr in [('No Opioid', '#95a5a6'), ('With Opioid', '#e67e22')]:
        for m in models:
            sub = df_op[df_op['OpioidLabel'] == op_lbl][f'{m}_MDPE'].dropna()
            pos.append(base); labs.append(m)
            data_boxes.append(sub.values); clrs.append(clr); base += 0.9
        base += 1.0
    bp_op = ax_op.boxplot(data_boxes, positions=pos, widths=0.75, patch_artist=True, notch=False,
                            medianprops=dict(color='black', linewidth=1.5), flierprops=dict(marker='.', alpha=0.2, markersize=3))
    for patch, clr in zip(bp_op['boxes'], clrs): patch.set_facecolor(clr); patch.set_alpha(0.7)
    ax_op.axhline(0, color='red', linewidth=1.5, linestyle='--', zorder=0)
    ax_op.set_xticks(pos); ax_op.set_xticklabels(labs, rotation=45, fontsize=8)
    ax_op.set_ylabel('MDPE (%)'); ax_op.set_title('Bias By Opioid Co-administration')
    ax_op.set_ylim(-100, 150)
    ax_op.text(0, -85, "Gray=No Opioid, Orange=With Opioid", color='black', fontsize=9, fontweight='bold')


    def plot_lowess_covariate(ax_obj, col_name, xlabel, text_under_x, text_under_y, text_over_y):
        ax_obj.axhline(0, color='black', linewidth=1.5, linestyle=':', alpha=0.7)
        for m in models:
            sub = df[[col_name, f'{m}_MDPE']].dropna()
            if len(sub) > 10:
                sm = lowess(sub[f'{m}_MDPE'], sub[col_name], frac=0.45, return_sorted=True)
                ax_obj.plot(sm[:, 0], sm[:, 1], color=MODEL_CFG[m]['color'], lw=3.0, label=m)
        ax_obj.set_xlabel(xlabel); ax_obj.set_ylabel('MDPE (%) / Bias')
        ax_obj.set_title(f'Bias vs {xlabel.split(" ")[0]} (LOWESS)'); ax_obj.legend(fontsize=9, labelcolor='black')
        ax_obj.set_ylim(-40, 40)
        ax_obj.text(text_under_x, text_under_y, "Under-predicting", color='gray', fontsize=10, fontweight='bold')
        ax_obj.text(text_under_x, text_over_y, "Over-predicting", color='gray', fontsize=10, fontweight='bold')

    ax2 = fig.add_subplot(gs[1, 0]); style_ax(ax2)
    plot_lowess_covariate(ax2, 'age', 'Age (years)', 10, 32, -35)

    ax3 = fig.add_subplot(gs[1, 1]); style_ax(ax3)
    plot_lowess_covariate(ax3, 'weight', 'Weight (kg)', 120, 32, -35)

    ax_ht = fig.add_subplot(gs[1, 2]); style_ax(ax_ht)
    plot_lowess_covariate(ax_ht, 'height', 'Height (cm)', 190, 32, -35)

    fig.suptitle('Model Bias Comparison (MDPE): Directionality of Prediction Errors\n'
                 'Evaluated on 1,033 patients. PE = (Measured - Predicted) / Predicted',
                 fontsize=14, fontweight='bold', color='black', y=0.98)
    plt.savefig(out_path, dpi=160, bbox_inches='tight', facecolor='#ffffff')
    plt.close()


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    print("Loading patient data...")
    t0 = time.time()
    patients = load_data()
    print(f"  {len(patients)} patients loaded in {time.time()-t0:.1f}s")

    theta_updated = list(theta_opt) if theta_opt is not None else list(THETA_DEFAULT)

    print("Evaluating all models (MDAPE, MDPE) with multiprocessing...")
    t0 = time.time()
    n_cores = min(mp.cpu_count(), 8)
    eval_fn = partial(evaluate_patient, theta_updated=theta_updated)
    with mp.Pool(processes=n_cores) as pool:
        records = pool.map(eval_fn, patients)
    print(f"  Done in {time.time()-t0:.1f}s")

    df = pd.DataFrame(records)

    print("\n{:<20} {:>12} {:>12} | {:>12} {:>12}".format('Model', 'MDAPE Mean', 'MDAPE Med', 'MDPE Mean', 'MDPE Med'))
    print('-' * 74)
    for m in MODEL_CFG.keys():
        mdape = df[f'{m}_MDAPE'].dropna()
        mdpe  = df[f'{m}_MDPE'].dropna()
        print("{:<20} {:>11.2f}% {:>11.2f}% | {:>11.2f}% {:>11.2f}%".format(
            m, mdape.mean(), mdape.median(), mdpe.mean(), mdpe.median()))

    print("\nGenerating MDAPE per-model figures...")
    for m, cfg in MODEL_CFG.items():
        make_single_model_figure(df, m, cfg['color'], os.path.join(FIG_DIR, f'mdape_{m}.png'))

    print("Generating MDAPE comparison figure...")
    make_mdape_comparison(df, os.path.join(FIG_DIR, 'mdape_comparison_all.png'))

    print("Generating MDPE (Bias) analysis figure...")
    make_mdpe_bias_figure(df, os.path.join(FIG_DIR, 'mdpe_bias_analysis.png'))

    print("\nAll figures written to:", FIG_DIR)

if __name__ == '__main__':
    main()
