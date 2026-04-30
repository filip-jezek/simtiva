"""
analyze_population_ranges.py

Population covariate statistics + visualization of Eleveld equations
with Hill sigmoid fits and extrapolation over the data ranges.
"""

import math
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

DATA_PATH = '/mnt/c/home/git/simtiva/data/supplementary_digital_content_1.txt'
OUT_PATH  = '/mnt/c/home/git/simtiva/PKPD_Reimplementation/eleveld_hill_analysis.png'

# ── Eleveld building blocks ──────────────────────────────────────────────────

def fsigmoid(x, y, z):
    return x**z / (x**z + y**z)

def fcentral(weight):
    return fsigmoid(weight, 33.6, 1)

def fclmaturation(pma):
    return fsigmoid(pma, 42.3, 9.06)

def fq3maturation(age_weeks):
    return fsigmoid(age_weeks + 40, 68.3, 1)

def fffm(weight, height, age, gender):
    bmi = weight / (height / 100.0) ** 2
    if gender == 1:  # male
        factor = 0.88 + (1-0.88) / (1 + (age/13.4)**(-12.7))
        return factor * (9270*weight) / (6680 + 216*bmi)
    else:
        factor = 1.11 + (1-1.11) / (1 + (age/7.1)**(-1.1))
        return factor * (9270*weight) / (8780 + 244*bmi)

# ── Hill function for fitting ────────────────────────────────────────────────

def hill_rising(x, base, emax, x50, gamma):
    """Rising Hill: starts at base, approaches base+emax."""
    return base + emax * x**gamma / (x50**gamma + x**gamma)

def hill_falling(x, top, drop, x50, gamma):
    """Falling Hill: starts at top, drops by 'drop' as x grows."""
    return top - drop * x**gamma / (x50**gamma + x**gamma)

# ── Load data ────────────────────────────────────────────────────────────────

df = pd.read_csv(DATA_PATH, sep=r'\s+')
demo = df.groupby('ID').first().reset_index()
demo['BMI'] = demo['WGT'] / (demo['HGT'] / 100)**2
demo['FFM'] = demo.apply(lambda r: fffm(r['WGT'], r['HGT'], r['AGE'], r['M1F2']), axis=1)
demo['PMA_weeks'] = demo['AGE'] * 52.1429 + 40

ffmref = (0.88 + (1-0.88)/(1+(35/13.4)**(-12.7))) * ((9270*70)/(6680+216*24.22145))

# ── Compute Eleveld outputs per patient ──────────────────────────────────────

records = []
for _, r in demo.iterrows():
    w = r['WGT']; a = r['AGE']; h = r['HGT']; g = int(r['M1F2']); pma = r['PMA_weeks']
    opioid = (r['TECH'] == 2)
    ffm = fffm(w, h, a, g)
    
    vc  = 6.28 * fcentral(w) / fcentral(70)
    v2  = 25.5 * (w/70) * math.exp(-0.0156*(a-35))
    v3  = 273 * ffm / ffmref * (math.exp(-0.0138*a) if opioid else 1.0)
    
    mat = fclmaturation(pma) / fclmaturation(35*52.1429+40)
    age_decay = math.exp(-0.00286*a)
    base_cl = 1.79 if g == 1 else 2.1
    cl1 = base_cl * (w/70)**0.75 * mat * (age_decay if opioid else 1.0)
    cl2 = 1.75 * (v2/25.5)**0.75 * (1 + 1.3*(1 - fq3maturation(a*52.1429)))
    cl3 = 1.11 * (v3/273)**0.75 * (fq3maturation(a*52.1429) / fq3maturation(35*52.1429))
    
    # Isolated sub-functions for plotting
    v2_age_factor = math.exp(-0.0156*(a - 35))
    v3_age_factor = math.exp(-0.0138*a) if opioid else 1.0
    cl1_age_factor = age_decay if opioid else 1.0
    cl_mat_val = fclmaturation(pma)
    q3_mat_val = fq3maturation(a * 52.1429)
    
    records.append({
        'age': a, 'weight': w, 'ffm': ffm, 'pma': pma, 'gender': g, 'opioid': opioid,
        'vc': vc, 'v2': v2, 'v3': v3, 'cl1': cl1, 'cl2': cl2, 'cl3': cl3,
        'v2_age_factor': v2_age_factor, 'v3_age_factor': v3_age_factor,
        'cl1_age_factor': cl1_age_factor, 'cl_mat': cl_mat_val, 'q3_mat': q3_mat_val,
        'vc_weight_factor': fcentral(w) / fcentral(70),
    })

pop = pd.DataFrame(records)

# ── Fit Hill functions to Eleveld outputs ────────────────────────────────────

# 1. Vc vs Weight — sigmoid(W, 33.6, 1) / sigmoid(70, 33.6, 1) * 6.28
wgts_fine = np.linspace(0.5, 180, 500)
vc_eleveld_fine = np.array([6.28 * fcentral(w) / fcentral(70) for w in wgts_fine])
popt_vc, _ = curve_fit(hill_rising, pop['weight'], pop['vc'], p0=[0, 6.28, 33.6, 1.0], maxfev=10000)

# 2. V2 age factor: exp(-0.0156*(Age-35))
ages_fine = np.linspace(0, 100, 500)
v2_age_eleveld = np.exp(-0.0156*(ages_fine - 35))
# Fit falling hill to the factor
popt_v2age, _ = curve_fit(hill_falling, pop['age'], pop['v2_age_factor'],
                           p0=[1.7, 1.3, 60, 2], maxfev=10000,
                           bounds=([0.5, 0.1, 10, 0.1], [5, 4, 150, 20]))

# 3. V3 age factor (opioid only): exp(-0.0138*age)
opi = pop[pop['opioid'] == True].copy()
v3_age_eleveld = np.exp(-0.0138*ages_fine)
popt_v3age, _ = curve_fit(hill_falling, opi['age'], opi['v3_age_factor'],
                           p0=[1.0, 0.8, 50, 2], maxfev=10000,
                           bounds=([0.3, 0.1, 5, 0.1], [2, 2, 150, 20]))

# 4. CL maturation vs PMA
pmas_fine = np.linspace(30, 5000, 500)
cl_mat_eleveld = np.array([fclmaturation(p) for p in pmas_fine])
popt_clmat, _ = curve_fit(hill_rising, pop['pma'], pop['cl_mat'], p0=[0, 1, 42.3, 9.06], maxfev=10000,
                           bounds=([0, 0.5, 1, 0.1], [0.01, 2, 200, 20]))

# 5. Q3 maturation vs Age (weeks)
q3_input_fine = ages_fine * 52.1429
q3_mat_eleveld = np.array([fq3maturation(aw) for aw in q3_input_fine])
popt_q3, _ = curve_fit(hill_rising, pop['age']*52.1429 + 40, pop['q3_mat'], 
                        p0=[0, 1, 68.3, 1.0], maxfev=10000,
                        bounds=([0, 0.5, 1, 0.1], [0.01, 2, 200, 20]))

# 6. CL1 age factor (opioid only): exp(-0.00286*age)
cl1_age_eleveld = np.exp(-0.00286*ages_fine)
opi_cl = pop[pop['opioid'] == True].copy()
popt_cl1age, _ = curve_fit(hill_falling, opi_cl['age'], opi_cl['cl1_age_factor'],
                            p0=[1.0, 0.3, 60, 2], maxfev=10000,
                            bounds=([0.5, 0.01, 10, 0.1], [2, 1, 200, 20]))

# ── PLOTTING ─────────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(22, 16))
fig.suptitle('Eleveld Equations: Population Data + Hill Sigmoid Fits & Extrapolation',
             fontsize=16, fontweight='bold', y=0.98)
gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.3)

def style_ax(ax, title, xlabel, ylabel):
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.legend(fontsize=9, loc='best')
    ax.grid(True, alpha=0.3)

# ── Panel 1: Vc vs Weight ───────────────────────────────────────────────────
ax1 = fig.add_subplot(gs[0, 0])
# Population scatter (color by sex)
males = pop[pop['gender'] == 1]
females = pop[pop['gender'] == 2]
ax1.scatter(males['weight'], males['vc'], alpha=0.3, s=12, c='steelblue', label='Males', zorder=2)
ax1.scatter(females['weight'], females['vc'], alpha=0.3, s=12, c='salmon', label='Females', zorder=2)
# Eleveld exact
ax1.plot(wgts_fine, vc_eleveld_fine, 'k-', linewidth=2, label='Eleveld: 6.28·σ(W,33.6,1)/σ(70,33.6,1)', zorder=3)
# Hill fit extrapolation
vc_hill_fine = hill_rising(wgts_fine, *popt_vc)
ax1.plot(wgts_fine, vc_hill_fine, 'r--', linewidth=2, 
         label=f'Hill fit: base={popt_vc[0]:.2f} Emax={popt_vc[1]:.2f} W50={popt_vc[2]:.1f} γ={popt_vc[3]:.2f}', zorder=3)
# Data range shading
ax1.axvspan(pop['weight'].min(), pop['weight'].max(), alpha=0.08, color='green', zorder=0)
ax1.axvline(pop['weight'].median(), color='green', ls=':', alpha=0.5, label=f'Median W={pop["weight"].median():.0f} kg')
style_ax(ax1, 'Vc (Central Volume) vs Weight', 'Weight (kg)', 'Vc (L)')

# ── Panel 2: V2 age factor ──────────────────────────────────────────────────
ax2 = fig.add_subplot(gs[0, 1])
ax2.scatter(males['age'], males['v2_age_factor'], alpha=0.3, s=12, c='steelblue', label='Males', zorder=2)
ax2.scatter(females['age'], females['v2_age_factor'], alpha=0.3, s=12, c='salmon', label='Females', zorder=2)
ax2.plot(ages_fine, v2_age_eleveld, 'k-', linewidth=2, label='Eleveld: exp(-0.0156·(Age-35))', zorder=3)
v2_hill_fine = hill_falling(ages_fine, *popt_v2age)
ax2.plot(ages_fine, v2_hill_fine, 'r--', linewidth=2,
         label=f'Hill: top={popt_v2age[0]:.2f} drop={popt_v2age[1]:.2f} A50={popt_v2age[2]:.1f} γ={popt_v2age[3]:.2f}', zorder=3)
ax2.axvspan(pop['age'].min(), pop['age'].max(), alpha=0.08, color='green', zorder=0)
ax2.axvline(pop['age'].median(), color='green', ls=':', alpha=0.5, label=f'Median Age={pop["age"].median():.0f} yr')
style_ax(ax2, '🔴 V2 Age Factor (candidate for Hill replacement)', 'Age (years)', 'V2 age scaling factor')

# ── Panel 3: V3 age factor (opioid) ─────────────────────────────────────────
ax3 = fig.add_subplot(gs[1, 0])
ax3.scatter(opi['age'], opi['v3_age_factor'], alpha=0.3, s=12, c='darkorange', label='Opioid patients', zorder=2)
ax3.plot(ages_fine, v3_age_eleveld, 'k-', linewidth=2, label='Eleveld: exp(-0.0138·Age)', zorder=3)
v3_hill_fine = hill_falling(ages_fine, *popt_v3age)
ax3.plot(ages_fine, v3_hill_fine, 'r--', linewidth=2,
         label=f'Hill: top={popt_v3age[0]:.2f} drop={popt_v3age[1]:.2f} A50={popt_v3age[2]:.1f} γ={popt_v3age[3]:.2f}', zorder=3)
ax3.axvspan(opi['age'].min(), opi['age'].max(), alpha=0.08, color='green', zorder=0)
ax3.axvline(opi['age'].median(), color='green', ls=':', alpha=0.5, label=f'Median Age={opi["age"].median():.0f} yr')
style_ax(ax3, '🔴 V3 Age-Opioid Factor (candidate for Hill replacement)', 'Age (years)', 'V3 age scaling factor')

# ── Panel 4: CL maturation vs PMA ───────────────────────────────────────────
ax4 = fig.add_subplot(gs[1, 1])
ax4.scatter(pop['pma'], pop['cl_mat'], alpha=0.3, s=12, c='purple', label='All patients', zorder=2)
ax4.plot(pmas_fine, cl_mat_eleveld, 'k-', linewidth=2, label='Eleveld: σ(PMA, 42.3, 9.06)', zorder=3)
cl_mat_hill = hill_rising(pmas_fine, *popt_clmat)
ax4.plot(pmas_fine, cl_mat_hill, 'r--', linewidth=2,
         label=f'Hill fit (perfect): PMA50={popt_clmat[2]:.1f} γ={popt_clmat[3]:.2f}', zorder=3)
ax4.axvspan(pop['pma'].min(), pop['pma'].max(), alpha=0.08, color='green', zorder=0)
ax4.set_xlim(-50, 300)  # zoom into interesting region
style_ax(ax4, '✅ CL Maturation vs PMA (already a Hill — keep)', 'PMA (weeks)', 'CL maturation factor')

# ── Panel 5: Q3 maturation vs Age ───────────────────────────────────────────
ax5 = fig.add_subplot(gs[2, 0])
ax5.scatter(pop['age'], pop['q3_mat'], alpha=0.3, s=12, c='teal', label='All patients', zorder=2)
ax5.plot(ages_fine, q3_mat_eleveld, 'k-', linewidth=2, label='Eleveld: σ(AgeWeeks+40, 68.3, 1)', zorder=3)
q3_hill_fine = hill_rising(q3_input_fine + 40, *popt_q3)
ax5.plot(ages_fine, q3_hill_fine, 'r--', linewidth=2,
         label=f'Hill fit (perfect): x50={popt_q3[2]:.1f} γ={popt_q3[3]:.2f}', zorder=3)
ax5.axvspan(pop['age'].min(), pop['age'].max(), alpha=0.08, color='green', zorder=0)
style_ax(ax5, '✅ Q3 Maturation vs Age (already a Hill — keep)', 'Age (years)', 'Q3 maturation factor')

# ── Panel 6: CL1 age factor (opioid) ────────────────────────────────────────
ax6 = fig.add_subplot(gs[2, 1])
ax6.scatter(opi_cl['age'], opi_cl['cl1_age_factor'], alpha=0.3, s=12, c='firebrick', label='Opioid patients', zorder=2)
ax6.plot(ages_fine, cl1_age_eleveld, 'k-', linewidth=2, label='Eleveld: exp(-0.00286·Age)', zorder=3)
cl1_hill_fine = hill_falling(ages_fine, *popt_cl1age)
ax6.plot(ages_fine, cl1_hill_fine, 'r--', linewidth=2,
         label=f'Hill: top={popt_cl1age[0]:.2f} drop={popt_cl1age[1]:.2f} A50={popt_cl1age[2]:.1f} γ={popt_cl1age[3]:.2f}', zorder=3)
ax6.axvspan(opi_cl['age'].min(), opi_cl['age'].max(), alpha=0.08, color='green', zorder=0)
ax6.axvline(opi_cl['age'].median(), color='green', ls=':', alpha=0.5, label=f'Median Age={opi_cl["age"].median():.0f} yr')
style_ax(ax6, '🟡 CL1 Age Factor (mild candidate — small effect)', 'Age (years)', 'CL1 age scaling factor')

plt.savefig(OUT_PATH, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {OUT_PATH}")

# ── Print summary table ─────────────────────────────────────────────────────
print()
print("Hill fit parameters summary:")
print(f"  Vc(W):       base={popt_vc[0]:.4f}  Emax={popt_vc[1]:.4f}  W50={popt_vc[2]:.2f}  γ={popt_vc[3]:.3f}")
print(f"  V2(Age) ↓:   top={popt_v2age[0]:.4f}  drop={popt_v2age[1]:.4f}  A50={popt_v2age[2]:.2f}  γ={popt_v2age[3]:.3f}")
print(f"  V3(Age) ↓:   top={popt_v3age[0]:.4f}  drop={popt_v3age[1]:.4f}  A50={popt_v3age[2]:.2f}  γ={popt_v3age[3]:.3f}")
print(f"  CLmat(PMA):   base={popt_clmat[0]:.4f}  Emax={popt_clmat[1]:.4f}  PMA50={popt_clmat[2]:.2f}  γ={popt_clmat[3]:.3f}")
print(f"  Q3mat(Age):   base={popt_q3[0]:.4f}  Emax={popt_q3[1]:.4f}  x50={popt_q3[2]:.2f}  γ={popt_q3[3]:.3f}")
print(f"  CL1(Age) ↓:  top={popt_cl1age[0]:.4f}  drop={popt_cl1age[1]:.4f}  A50={popt_cl1age[2]:.2f}  γ={popt_cl1age[3]:.3f}")
