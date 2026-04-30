"""
plot_sigmoid_fits.py

Generates comparison figures of the original Eleveld exponential age-decay
functions vs the Hill sigmoidal replacements for V2 and V3, extrapolated
well beyond the training data range (age 0-100 years).

If optimized_theta.npy exists, also overlays the post-optimization curves.

Output:
    figures/sigmoid_comparison_V2.png
    figures/sigmoid_comparison_V3.png
    figures/sigmoid_comparison_combined.png
"""

import os
import sys
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(SCRIPT_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

sys.path.insert(0, SCRIPT_DIR)
from models.eleveld_updated import THETA_DEFAULT

# ── load optimized theta if available ────────────────────────────────────────
THETA_NPY = os.path.join(SCRIPT_DIR, 'optimized_theta.npy')
theta_opt = None
if os.path.isfile(THETA_NPY):
    theta_opt = np.load(THETA_NPY)
    print(f"Loaded optimized theta from {THETA_NPY}")
else:
    print("No optimized_theta.npy found — plotting default sigmoid only.")

# ── age axis: 0-105 years, highlight training range 0-88 ─────────────────────
AGE_WIDE   = np.linspace(0.1, 105, 500)    # extrapolation to 105 yr
AGE_TRAIN  = 88.0                           # max age in Eleveld dataset

# ── Eleveld exponential age functions ─────────────────────────────────────────
def eleveld_V2_age(age):
    return np.exp(-0.0156 * (age - 35))

def eleveld_V3_age(age):
    return np.exp(-0.0138 * age)

# ── Hill (sigmoidal) age functions ────────────────────────────────────────────
def hill_falling(age, top, drop, a50, gamma):
    return top - drop * age**gamma / (a50**gamma + age**gamma)

# ── V2 relative scaling: both functions normalised at age=35 ─────────────────
def get_v2_curves(theta):
    exp_curve  = eleveld_V2_age(AGE_WIDE)
    hill_curve = hill_falling(AGE_WIDE, theta[15], theta[16], theta[17], theta[18])
    # Normalise hill to match exp at age=35 for visual comparison
    hill_at35  = hill_falling(35, theta[15], theta[16], theta[17], theta[18])
    exp_at35   = eleveld_V2_age(35)
    scale      = exp_at35 / hill_at35 if hill_at35 > 0 else 1.0
    return exp_curve, hill_curve * scale

def get_v3_curves(theta):
    exp_curve  = eleveld_V3_age(AGE_WIDE)
    hill_curve = hill_falling(AGE_WIDE, theta[19], theta[20], theta[21], theta[22])
    hill_at35  = hill_falling(35, theta[19], theta[20], theta[21], theta[22])
    exp_at35   = eleveld_V3_age(35)
    scale      = exp_at35 / hill_at35 if hill_at35 > 0 else 1.0
    return exp_curve, hill_curve * scale

# ── absolute predictions for a reference patient: 70 kg, male ────────────────
WEIGHT_REF = 70.0
FFMREF = (
    (0.88 + (1 - 0.88) / (1 + (35 / 13.4) ** (-12.7)))
    * (9270 * 70) / (6680 + 216 * (70 / 1.7**2))
)

def abs_v2(theta, age_arr):
    return theta[3] * (WEIGHT_REF / 70) * hill_falling(age_arr, theta[15], theta[16], theta[17], theta[18])

def abs_v2_exp(age_arr):
    return 25.5 * (WEIGHT_REF / 70) * eleveld_V2_age(age_arr)

def abs_v3(theta, age_arr, ffm=FFMREF):
    return theta[4] * ffm / FFMREF * hill_falling(age_arr, theta[19], theta[20], theta[21], theta[22])

def abs_v3_exp(age_arr, ffm=FFMREF):
    return 273.0 * ffm / FFMREF * eleveld_V3_age(age_arr)

# ──────────────────────────────────────────────────────────────────────────────
# FIGURE : Combined 2×2 layout
# ──────────────────────────────────────────────────────────────────────────────

fig = plt.figure(figsize=(14, 10))
fig.patch.set_facecolor('#ffffff')
plt.rcParams.update({
    'text.color': 'black', 'axes.labelcolor': 'black',
    'xtick.color': 'black', 'ytick.color': 'black',
    'axes.titlecolor': 'black',
})

gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.40, wspace=0.35)
gs.update(left=0.08, right=0.97, top=0.90, bottom=0.09)

ALPHA_SHADE = 0.12
COLOR_EXP   = '#ff6b6b'   # red — original Eleveld
COLOR_SIG_D = '#4ecdc4'   # teal — default sigmoid (Eleveld-equiv init)
COLOR_SIG_O = '#ffd93d'   # gold — optimized sigmoid

def style_ax(ax):
    ax.set_facecolor('#ffffff')
    for spine in ax.spines.values():
        spine.set_edgecolor('#cccccc')
    ax.grid(True, alpha=0.3, color='gray', linewidth=0.5)
    ax.axvline(AGE_TRAIN, color='black', linewidth=1.5, linestyle=':', alpha=0.6,
               label='Max training age (88 yr)')


# ── PANEL A: Merged Relative Age Factors ──────────────────────────────────────
v2_exp, v2_sig_d = get_v2_curves(THETA_DEFAULT)
_, v2_sig_o = get_v2_curves(theta_opt) if theta_opt is not None else (None, np.zeros_like(AGE_WIDE))
v3_exp, v3_sig_d = get_v3_curves(THETA_DEFAULT)
_, v3_sig_o = get_v3_curves(theta_opt) if theta_opt is not None else (None, np.zeros_like(AGE_WIDE))

ax1 = fig.add_subplot(gs[0, :])
style_ax(ax1)
ax1.plot(AGE_WIDE, v2_exp, color='#2c3e50', lw=2.5, linestyle='-', label='V₂ Orig Exponential')
ax1.plot(AGE_WIDE, v2_sig_d, color='#2c3e50', lw=2.5, linestyle='--', alpha=0.6, label='V₂ Hill (Init)')
if theta_opt is not None:
    ax1.plot(AGE_WIDE, v2_sig_o, color='#16a085', lw=3.0, label='V₂ Hill (Optimized)')

ax1.plot(AGE_WIDE, v3_exp, color='#8e44ad', lw=2.5, linestyle='-', label='V₃ Orig Exponential (Opioid)')
ax1.plot(AGE_WIDE, v3_sig_d, color='#8e44ad', lw=2.5, linestyle='--', alpha=0.6, label='V₃ Hill (Init)')
if theta_opt is not None:
    ax1.plot(AGE_WIDE, v3_sig_o, color='#e67e22', lw=3.0, label='V₃ Hill (Optimized)')

ax1.fill_between(AGE_WIDE[AGE_WIDE <= AGE_TRAIN], 0,
                 np.ones(sum(AGE_WIDE <= AGE_TRAIN)) * 2.0, alpha=ALPHA_SHADE, color='black')
ax1.set_xlim(0, 105); ax1.set_ylim(0.0, 1.8)
ax1.set_xlabel('Age (years)', fontsize=11); ax1.set_ylabel('Age Decay Multiplier / Factor', fontsize=11)
ax1.set_title('A  |  Relative Age-Decay Factor for Structural Volumes (V₂ and V₃)', fontsize=13)
ax1.legend(loc='upper right', ncol=2, fontsize=10, labelcolor='black', framealpha=0.9)

# ── PANEL B: Absolute V2 ──────────────────────────────────────────────────────
ax3 = fig.add_subplot(gs[1, 0])
style_ax(ax3)
ax3.plot(AGE_WIDE, abs_v2_exp(AGE_WIDE), color='#2c3e50', lw=2.5, label='Orig Eleveld (Exponential)')
ax3.plot(AGE_WIDE, abs_v2(THETA_DEFAULT, AGE_WIDE), color='#2c3e50', lw=2.5, linestyle='--', alpha=0.6, label='Eleveld Updated (Init)')
if theta_opt is not None:
    ax3.plot(AGE_WIDE, abs_v2(theta_opt, AGE_WIDE), color='#16a085', lw=3.0, label='Eleveld Updated (Optimized)')
ax3.fill_between(AGE_WIDE[AGE_WIDE <= AGE_TRAIN], 0,
                 abs_v2_exp(AGE_WIDE[AGE_WIDE <= AGE_TRAIN]).max() * 1.1,
                 alpha=ALPHA_SHADE, color='black')
ax3.axvline(AGE_TRAIN, color='black', linewidth=1.5, linestyle=':', alpha=0.6)
ax3.set_xlim(0, 105); ax3.set_ylim(bottom=0)
ax3.set_xlabel('Age (years)'); ax3.set_ylabel('V₂ (L)')
ax3.set_title('C  |  Absolute V₂  (70 kg, male reference patient)')
ax3.legend(fontsize=7.5)

# ── PANEL C: Absolute V3 (opioid) ─────────────────────────────────────────────
ax4 = fig.add_subplot(gs[1, 1])
style_ax(ax4)
if v3_exp is not None:
    ax4.plot(AGE_WIDE, abs_v3_exp(AGE_WIDE), color='#8e44ad', lw=2.5, label='Orig Eleveld (Exponential)')
    ax4.plot(AGE_WIDE, abs_v3(THETA_DEFAULT, AGE_WIDE), color='#8e44ad', lw=2.5, linestyle='--', alpha=0.6, label='Eleveld Updated (Init)')
if theta_opt is not None:
    ax4.plot(AGE_WIDE, abs_v3(theta_opt, AGE_WIDE), color='#e67e22', lw=3.0, label='Eleveld Updated (Optimized)')
ax4.fill_between(AGE_WIDE[AGE_WIDE <= AGE_TRAIN], 0,
                 abs_v3_exp(AGE_WIDE[AGE_WIDE <= AGE_TRAIN]).max() * 1.1,
                 alpha=ALPHA_SHADE, color='black')
ax4.axvline(AGE_TRAIN, color='black', linewidth=1.5, linestyle=':', alpha=0.6)
ax4.set_xlim(0, 105); ax4.set_ylim(bottom=0)
ax4.set_xlabel('Age (years)'); ax4.set_ylabel('V₃ (L)')
ax4.set_title('D  |  Absolute V₃  (70 kg male, with opioid)')
ax4.legend(fontsize=7.5)

# ── Suptitle + annotation ─────────────────────────────────────────────────────
fig.suptitle(
    'Exponential vs. Hill Sigmoid Age-Decay Functions in the Eleveld PK Model\n'
    'Shaded region = training data range (age 0–88 yr);  dashed line = extrapolation boundary',
    fontsize=12, fontweight='bold', color='black', y=0.97
)

out_path = os.path.join(FIG_DIR, 'sigmoid_comparison_combined.png')
plt.savefig(out_path, dpi=180, bbox_inches='tight', facecolor=fig.get_facecolor())
plt.close()
print(f"Saved: {out_path}")
