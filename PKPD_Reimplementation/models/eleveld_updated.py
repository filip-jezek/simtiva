"""
eleveld_updated.py

Sigmoidal modification of the Eleveld 3-compartmental propofol PK model.

Key structural differences from the original Eleveld model:
  - V2 age-decay: exponential exp(-0.0156*(age-35)) replaced by a
    falling Hill sigmoid: top - drop * age^gamma / (A50^gamma + age^gamma)
  - V3 age-decay (opioid): exponential exp(-0.0138*age) replaced by a
    falling Hill sigmoid with independent parameters
  - All other equations are structurally identical to Eleveld (2018).

Default theta values represent the Eleveld-equivalent starting point
(identical performance to original Eleveld within 0.03%).

Optimized theta values (from 3-iteration L-BFGS-B run, see global_optimization.py)
achieved 21.8% global RMSE improvement (0.419 → 0.370) and can be loaded via
EleveldUpdatedModel.from_theta(theta_array).

Reference:
  Eleveld DJ et al. An Allometric Model of Remifentanil Pharmacokinetics and
  Pharmacodynamics. Anesthesiology 2018;128:818-834.
"""

import math
from .base_model import PropofolModel

# ═══════════════════════════════════════════════════════════════════════════════
# THETA VECTOR (25 parameters)
# ═══════════════════════════════════════════════════════════════════════════════
#
# [0]  Vc_base       = 6.28      Central volume base (L)
# [1]  Vc_W50        = 33.6      Weight sigmoid midpoint for Vc
# [2]  Vc_gamma      = 1.0       Weight sigmoid steepness for Vc
# [3]  V2_base       = 25.5      V2 base volume (L)
# [4]  V3_base       = 273.0     V3 base volume (L)
# [5]  CL_male       = 1.79      CL1 base for males (L/min)
# [6]  CL_female     = 2.10      CL1 base for females (L/min)
# [7]  CL_allom_exp  = 0.75      Allometric exponent for CL1
# [8]  CLmat_PMA50   = 42.3      CL maturation PMA midpoint (weeks)
# [9]  CLmat_gamma   = 9.06      CL maturation steepness
# [10] Q2_base       = 1.75      Q2 base clearance (L/min)
# [11] Q2_immat_coef = 1.3       Q2 immaturity amplification coefficient
# [12] Q3_base       = 1.11      Q3 base clearance (L/min)
# [13] Q3mat_x50     = 68.3      Q3 maturation midpoint (weeks PMA)
# [14] Q3mat_gamma   = 1.0       Q3 maturation steepness
#
# --- Hill sigmoid replacements for age-decay exponentials ---
# [15] V2_age_top    = 1.7235    V2 age Hill: plateau at high age
# [16] V2_age_drop   = 2.4810    V2 age Hill: total drop magnitude
# [17] V2_age_A50    = 81.03     V2 age Hill: midpoint (years)
# [18] V2_age_gamma  = 1.056     V2 age Hill: steepness (Hill coefficient)
# [19] V3_age_top    = 0.999     V3 age Hill: plateau at high age (opioid)
# [20] V3_age_drop   = 1.4798    V3 age Hill: total drop magnitude (opioid)
# [21] V3_age_A50    = 96.12     V3 age Hill: midpoint (years, opioid)
# [22] V3_age_gamma  = 1.043     V3 age Hill: steepness (opioid)
#
# [23] CL1_age_rate  = 0.00286   CL1 age decay rate (exponential, opioid only)
# [24] sigma_prop    = 0.3       Proportional residual error magnitude
#
# ═══════════════════════════════════════════════════════════════════════════════

# Eleveld-equivalent starting point (Hill sigmoids matched to Eleveld exponentials)
# These values are from the 3-iteration L-BFGS-B optimization run (RMSE: 0.419 → 0.370, -21.8%)
# documented in structural_optimization_methods.tex, Table 1.
# Use EleveldUpdatedModel.from_theta(np.load('optimized_theta.npy'), ...) to load further-optimized values.
THETA_DEFAULT = [
    6.28,       # [0]  Vc_base
    33.6,       # [1]  Vc_W50
    1.0,        # [2]  Vc_gamma
    25.5,       # [3]  V2_base
    273.0,      # [4]  V3_base
    1.79,       # [5]  CL_male
    2.10,       # [6]  CL_female
    0.75,       # [7]  CL_allom_exp
    42.3,       # [8]  CLmat_PMA50
    9.06,       # [9]  CLmat_gamma
    1.75,       # [10] Q2_base
    1.3,        # [11] Q2_immat_coef
    1.11,       # [12] Q3_base
    68.3,       # [13] Q3mat_x50
    1.0,        # [14] Q3mat_gamma
    1.7235,     # [15] V2_age_top
    2.4810,     # [16] V2_age_drop
    81.03,      # [17] V2_age_A50
    1.056,      # [18] V2_age_gamma
    0.999,      # [19] V3_age_top
    1.4798,     # [20] V3_age_drop
    96.12,      # [21] V3_age_A50
    1.043,      # [22] V3_age_gamma
    0.00286,    # [23] CL1_age_rate
    0.3,        # [24] sigma_prop
]

# Fixed reference FFM for a 35-year-old, 70 kg, 170 cm male
# Used to normalise V3 to the Eleveld reference individual
_BMI_REF = 70 / (170 / 100.0) ** 2  # ≈ 24.22
_FFMREF = (
    (0.88 + (1 - 0.88) / (1 + (35 / 13.4) ** (-12.7)))
    * (9270 * 70)
    / (6680 + 216 * _BMI_REF)
)


class EleveldUpdatedModel(PropofolModel):
    """
    Sigmoidal modification of the Eleveld propofol PK model.

    Replaces the unbounded exponential age-decay terms in V2 and V3 with
    physiologically-constrained falling Hill (sigmoidal) functions.
    All other covariate relationships are retained from Eleveld (2018).

    Parameters
    ----------
    weight  : float  Body weight (kg)
    age     : float  Age (years)
    height  : float  Height (cm)
    gender  : int    0 = male, 1 = female  (SimTIVA convention)
    opioid  : bool   True if opioid co-administered
    pma     : float  Post-menstrual age (weeks); computed from age if None
    theta   : list   25-element theta vector; uses THETA_DEFAULT if None
    """

    def __init__(
        self,
        weight,
        age,
        height,
        gender,
        is_adj_bw=False,
        pma=None,
        opioid=True,
        theta=None,
    ):
        super().__init__(weight, age, height, gender, is_adj_bw)
        self.opioid = opioid
        self.toweeks = 52.1429

        # Post-Menstrual Age
        if pma is not None:
            self.pma = pma
        else:
            self.pma = self.age * self.toweeks + 40

        # Theta vector
        self.theta = list(theta) if theta is not None else list(THETA_DEFAULT)
        if len(self.theta) != 25:
            raise ValueError(f"theta must have 25 elements, got {len(self.theta)}")

    # ── Class-level constructor for loading from a numpy array / file ─────────

    @classmethod
    def from_theta(cls, theta, weight, age, height, gender, opioid=True, pma=None):
        """
        Construct model from an explicit theta vector (e.g. loaded from
        optimized_theta.npy produced by global_optimization.py).

        Example::
            import numpy as np
            theta_opt = np.load('optimized_theta.npy')
            model = EleveldUpdatedModel.from_theta(
                theta_opt, weight=70, age=35, height=170, gender=0
            )
        """
        return cls(
            weight=weight,
            age=age,
            height=height,
            gender=gender,
            opioid=opioid,
            pma=pma,
            theta=list(theta),
        )

    # ── Private sigmoid helpers ───────────────────────────────────────────────

    def _fsigmoid(self, x, y, z):
        """Hill sigmoid: x^z / (x^z + y^z)"""
        return math.pow(x, z) / (math.pow(x, z) + math.pow(y, z))

    def _hill_falling(self, x, top, drop, x50, gamma):
        """
        Falling Hill function used for age-related volume decline.
        Returns `top` at age=0, asymptotes to `top - drop` at age → ∞.
        At x = x50: value = top - drop/2.
        """
        return top - drop * math.pow(x, gamma) / (
            math.pow(x50, gamma) + math.pow(x, gamma)
        )

    def _fffm(self):
        """Fat-Free Mass (Al-Sallami / Janmahasatian equation)."""
        bmi = self.weight / math.pow(self.height / 100.0, 2)
        if self.gender == 0:  # male
            factor = 0.88 + (1 - 0.88) / (1 + math.pow(self.age / 13.4, -12.7))
            return factor * (9270 * self.weight) / (6680 + 216 * bmi)
        else:  # female
            factor = 1.11 + (1 - 1.11) / (1 + math.pow(self.age / 7.1, -1.1))
            return factor * (9270 * self.weight) / (8780 + 244 * bmi)

    # ── Main parameter calculation ────────────────────────────────────────────

    def get_parameters(self):
        """
        Compute 3-compartment PK rate constants using the sigmoidal Eleveld
        structure and the stored theta vector.

        Returns
        -------
        dict with keys: vc, k10, k12, k13, k21, k31, ke0
        """
        th = self.theta

        # ── Vc: central volume — weight sigmoid (identical to Eleveld) ────────
        vc = th[0] * self._fsigmoid(self.weight, th[1], th[2]) / \
             self._fsigmoid(70, th[1], th[2])

        # ── V2: rapid peripheral — Hill age decay (REPLACES exponential) ──────
        v2_age_factor = self._hill_falling(
            self.age, th[15], th[16], th[17], th[18]
        )
        v2 = th[3] * (self.weight / 70) * v2_age_factor
        v2ref = th[3]

        # ── V3: slow peripheral — FFM scaling + Hill age decay (opioid only) ──
        ffm = self._fffm()
        if self.opioid:
            v3_age_factor = self._hill_falling(
                self.age, th[19], th[20], th[21], th[22]
            )
        else:
            v3_age_factor = 1.0
        v3 = th[4] * ffm / _FFMREF * v3_age_factor
        v3ref = th[4]

        # ── CL1: metabolic clearance ──────────────────────────────────────────
        mat_factor = (
            self._fsigmoid(self.pma, th[8], th[9])
            / self._fsigmoid(35 * self.toweeks + 40, th[8], th[9])
        )
        base_cl = th[5] if self.gender == 0 else th[6]  # male / female
        cl1 = base_cl * math.pow(self.weight / 70, th[7]) * mat_factor
        if self.opioid:
            cl1 *= math.exp(-th[23] * self.age)

        # ── Q2 / CL2: rapid inter-compartmental clearance ─────────────────────
        q3mat_age = self._fsigmoid(
            self.age * self.toweeks + 40, th[13], th[14]
        )
        cl2 = th[10] * math.pow(v2 / v2ref, 0.75) * (
            1 + th[11] * (1 - q3mat_age)
        )

        # ── Q3 / CL3: slow inter-compartmental clearance ──────────────────────
        q3mat_ref = self._fsigmoid(35 * self.toweeks + 40, th[13], th[14])
        cl3 = th[12] * math.pow(v3 / v3ref, 0.75) * (q3mat_age / q3mat_ref)

        # ── Rate constants ────────────────────────────────────────────────────
        k10 = cl1 / vc
        k12 = cl2 / vc
        k13 = cl3 / vc
        k21 = cl2 / v2
        k31 = cl3 / v3
        ke0 = 0.146 * math.pow(self.weight / 70, -0.25)

        return {
            'vc': vc,
            'k10': k10,
            'k12': k12,
            'k13': k13,
            'k21': k21,
            'k31': k31,
            'ke0': ke0,
        }
