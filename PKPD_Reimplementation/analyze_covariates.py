"""
analyze_covariates.py
Covariate-stratified RMSE analysis: how well does a PK model fit
each patient, as a function of their clinical characteristics?

Usage:
  python3 analyze_covariates.py --model Marsh
  python3 analyze_covariates.py --model ALL
"""
import os
import re
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from statsmodels.nonparametric.smoothers_lowess import lowess

warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR   = os.path.join(SCRIPT_DIR, '..', 'data')
SIM_DIR    = os.path.join(SCRIPT_DIR, 'simulated_patients')

MODEL_COLORS = {
    'Marsh':    'firebrick',
    'Schnider': 'steelblue',
    'Paedfusor':'forestgreen',
    'Eleveld':  'purple',
    'ENONMEM':  'darkorange',
    'ELEVELD_PD':'teal',
}

# ── 1. Parse summary.txt ──────────────────────────────────────────────────────

def parse_summary(summary_path):
    """Returns DataFrame with columns: ID, model_name, RMSE"""
    rows = []
    with open(summary_path) as f:
        for line in f:
            line = line.strip()
            if not line.startswith('P'):
                continue
            parts = line.split()
            if not parts or not re.match(r'^P\d{4}$', parts[0]):
                continue
            pid = int(parts[0][1:])
            for tok in parts[1:]:
                m = re.match(r'(\w+)=([\d.]+)', tok)
                if m:
                    rows.append({'ID': pid, 'model': m.group(1), 'RMSE': float(m.group(2))})
    return pd.DataFrame(rows)


# ── 2. Load demographics ──────────────────────────────────────────────────────

def load_demographics():
    df = pd.read_csv(os.path.join(DATA_DIR, 'supplementary_digital_content_1.txt'), sep=r'\s+')
    demo = df.groupby('ID').first().reset_index()
    demo['BMI']    = demo['WGT'] / (demo['HGT'] / 100) ** 2
    demo['Sex']    = demo['M1F2'].map({1: 'Male', 2: 'Female'})
    demo['Opioid'] = demo['TECH'].map({1: False, 2: True})
    demo['AgeGroup'] = pd.cut(demo['AGE'], bins=[0, 35, 60, 200],
                              labels=['Young (<35)', 'Adult (35-60)', 'Senior (>60)'])
    demo['WgtGroup'] = pd.cut(demo['WGT'], bins=[0, 60, 80, 200],
                              labels=['Light (<60)', 'Normal (60-80)', 'Heavy (>80)'])
    return demo[['ID', 'AGE', 'WGT', 'HGT', 'BMI', 'Sex', 'Opioid',
                 'STDY', 'AgeGroup', 'WgtGroup']].copy()


# ── 3. Single-model analysis figure ──────────────────────────────────────────

def make_figure(df_m, model_name, out_path):
    """df_m: merged DataFrame with RMSE + demographics for one model."""
    color = MODEL_COLORS.get(model_name, 'gray')
    fig = plt.figure(figsize=(20, 14))
    fig.suptitle(f'{model_name} model — RMSE vs Patient Covariates  (N={len(df_m)})',
                 fontsize=16, fontweight='bold', y=0.98)

    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.35)

    # ── panel helper: scatter + LOWESS ───────────────────────────────────────
    def scatter_lowess(ax, xcol, xlabel):
        for sex, marker, clr in [('Male', 'o', 'steelblue'), ('Female', '^', 'salmon')]:
            sub = df_m[df_m['Sex'] == sex]
            ax.scatter(sub[xcol], sub['RMSE'], alpha=0.25, s=14,
                       color=clr, marker=marker, label=sex)
            if len(sub) > 10:
                sm = lowess(sub['RMSE'], sub[xcol], frac=0.4, return_sorted=True)
                ax.plot(sm[:, 0], sm[:, 1], color=clr, linewidth=2.5)
        rho, pval = stats.spearmanr(df_m[xcol], df_m['RMSE'])
        ax.set_xlabel(xlabel)
        ax.set_ylabel('RMSE (mg/L)')
        ax.set_title(f'vs {xlabel}  (ρ={rho:+.2f}, p={pval:.3f})')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Panel 1 — Age
    scatter_lowess(fig.add_subplot(gs[0, 0]), 'AGE', 'Age (years)')
    # Panel 2 — Weight
    scatter_lowess(fig.add_subplot(gs[0, 1]), 'WGT', 'Weight (kg)')
    # Panel 3 — BMI
    scatter_lowess(fig.add_subplot(gs[0, 2]), 'BMI', 'BMI (kg/m²)')

    # ── Panel 4: Age × Weight 2-D heatmap ───────────────────────────────────
    ax4 = fig.add_subplot(gs[1, 0])
    age_bins = np.arange(0, 101, 10)
    wgt_bins = np.arange(30, 151, 10)
    df_m['AgeBin'] = pd.cut(df_m['AGE'], age_bins)
    df_m['WgtBin'] = pd.cut(df_m['WGT'], wgt_bins)
    heat = df_m.groupby(['AgeBin', 'WgtBin'], observed=True)['RMSE'].median().unstack('WgtBin')
    im = ax4.imshow(heat.values, aspect='auto', origin='lower', cmap='RdYlGn_r')
    plt.colorbar(im, ax=ax4, label='Median RMSE')
    ax4.set_xticks(range(len(heat.columns)))
    ax4.set_xticklabels([str(c.mid) for c in heat.columns], rotation=45, fontsize=7)
    ax4.set_yticks(range(len(heat.index)))
    ax4.set_yticklabels([str(int(c.mid)) for c in heat.index], fontsize=7)
    ax4.set_xlabel('Weight (kg)')
    ax4.set_ylabel('Age (years)')
    ax4.set_title('Median RMSE: Age × Weight')

    # ── Panel 5: Sex × AgeGroup box plot ─────────────────────────────────────
    ax5 = fig.add_subplot(gs[1, 1])
    groups = [('Male', 'steelblue'), ('Female', 'salmon')]
    age_grps = df_m['AgeGroup'].cat.categories.tolist()
    positions, tick_labels, all_data, colors_list = [], [], [], []
    base = 0
    for age_g in age_grps:
        for i, (sex, clr) in enumerate(groups):
            sub = df_m[(df_m['AgeGroup'] == age_g) & (df_m['Sex'] == sex)]['RMSE'].dropna()
            pos = base + i * 0.9
            positions.append(pos)
            tick_labels.append(f'{sex[0]} {age_g.split()[0]}')
            all_data.append(sub.values)
            colors_list.append(clr)
        base += 2.5
    bp = ax5.boxplot(all_data, positions=positions, widths=0.75, patch_artist=True, notch=False)
    for patch, clr in zip(bp['boxes'], colors_list):
        patch.set_facecolor(clr); patch.set_alpha(0.7)
    ax5.set_xticks(positions)
    ax5.set_xticklabels(tick_labels, rotation=45, fontsize=7)
    ax5.set_ylabel('RMSE (mg/L)')
    ax5.set_title('RMSE by Sex × Age Group')
    ax5.grid(True, axis='y', alpha=0.3)
    # Mann-Whitney p for sex difference
    mw_u, mw_p = stats.mannwhitneyu(
        df_m[df_m['Sex'] == 'Male']['RMSE'].dropna(),
        df_m[df_m['Sex'] == 'Female']['RMSE'].dropna())
    ax5.text(0.97, 0.97, f'Sex MW p={mw_p:.3f}', transform=ax5.transAxes,
             ha='right', va='top', fontsize=8)

    # ── Panel 6: Opioid & Study site ─────────────────────────────────────────
    ax6 = fig.add_subplot(gs[1, 2])
    opi_yes = df_m[df_m['Opioid'] == True]['RMSE'].dropna()
    opi_no  = df_m[df_m['Opioid'] == False]['RMSE'].dropna()
    bp6 = ax6.boxplot([opi_no.values, opi_yes.values], patch_artist=True,
                      labels=['No opioid', 'With opioid'])
    bp6['boxes'][0].set_facecolor('lightblue');  bp6['boxes'][0].set_alpha(0.8)
    bp6['boxes'][1].set_facecolor('lightsalmon'); bp6['boxes'][1].set_alpha(0.8)
    _, p6 = stats.mannwhitneyu(opi_no, opi_yes)
    ax6.set_ylabel('RMSE (mg/L)')
    ax6.set_title(f'RMSE by Opioid Co-admin  (MW p={p6:.3f})')
    ax6.grid(True, axis='y', alpha=0.3)

    # ── Panel 7: Spearman correlation bar chart ───────────────────────────────
    ax7 = fig.add_subplot(gs[2, :2])
    cov_cols = [('AGE', 'Age'), ('WGT', 'Weight'), ('HGT', 'Height'), ('BMI', 'BMI'),
                ('Opioid', 'Opioid (bool)'), ('STDY', 'Study Site')]
    rhos, pvals, labels_c = [], [], []
    for col, label in cov_cols:
        x = df_m[col].astype(float)
        mask = x.notna() & df_m['RMSE'].notna()
        rho, pval = stats.spearmanr(x[mask], df_m['RMSE'][mask])
        rhos.append(rho); pvals.append(pval); labels_c.append(label)
    # Sort by |rho|
    order = np.argsort(np.abs(rhos))[::-1]
    rhos_s = [rhos[i] for i in order]
    pvs_s  = [pvals[i] for i in order]
    labs_s = [labels_c[i] for i in order]
    bar_colors = ['firebrick' if r > 0 else 'steelblue' for r in rhos_s]
    bars = ax7.barh(range(len(rhos_s)), rhos_s, color=bar_colors, alpha=0.8)
    ax7.set_yticks(range(len(labs_s)))
    ax7.set_yticklabels(labs_s)
    ax7.axvline(0, color='black', linewidth=0.8)
    for i, (r, p) in enumerate(zip(rhos_s, pvs_s)):
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        ax7.text(r + (0.005 if r >= 0 else -0.005), i, f'{r:+.3f}{sig}',
                 va='center', ha='left' if r >= 0 else 'right', fontsize=9)
    ax7.set_xlabel('Spearman ρ with RMSE')
    ax7.set_title('Covariate Correlation with RMSE (ranked)')
    ax7.grid(True, axis='x', alpha=0.3)

    # ── Panel 8: Study site box ───────────────────────────────────────────────
    ax8 = fig.add_subplot(gs[2, 2])
    sites = sorted(df_m['STDY'].dropna().unique())
    site_data = [df_m[df_m['STDY'] == s]['RMSE'].dropna().values for s in sites]
    bp8 = ax8.boxplot(site_data, patch_artist=True, labels=[f'S{int(s)}' for s in sites])
    for patch in bp8['boxes']:
        patch.set_facecolor(color); patch.set_alpha(0.6)
    ax8.set_ylabel('RMSE (mg/L)')
    ax8.set_title('RMSE by Study Site')
    ax8.tick_params(axis='x', labelsize=8)
    ax8.grid(True, axis='y', alpha=0.3)

    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_path}")


# ── 4. Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', default='Marsh',
                        help='Model name or ALL. E.g. Marsh, Eleveld, ALL')
    args = parser.parse_args()

    summary_path = os.path.join(SIM_DIR, 'summary.txt')
    rmse_df  = parse_summary(summary_path)
    demo_df  = load_demographics()

    available_models = rmse_df['model'].unique().tolist()
    models_to_run = available_models if args.model == 'ALL' else [args.model]

    for model in models_to_run:
        sub = rmse_df[rmse_df['model'] == model].copy()
        if sub.empty:
            print(f"No data for model '{model}'. Available: {available_models}")
            continue

        merged = sub.merge(demo_df, on='ID', how='inner')
        out_path = os.path.join(SIM_DIR, f'covariate_analysis_{model}.png')
        print(f"Processing {model} (N={len(merged)})...")
        make_figure(merged, model, out_path)


if __name__ == '__main__':
    main()
