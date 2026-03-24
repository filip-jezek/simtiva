import math
from .base_model import PropofolModel

class EleveldNonMem(PropofolModel):
    def __init__(self, weight, age, height, gender, pma=None, tech=2, a1v2=1, ev1=None, ev2=None, ev3=None, ecl=None, eq2=None, eq3=None):
        super().__init__(weight, age, height, gender, False)
        
        # NONMEM Covariates mapping
        self.wgt = weight
        self.age = age
        self.hgt = height
        self.m1f2 = gender # 1 = Male, 2 = Female
        self.tech = tech # 1 = Propofol alone, 2 = Opiates
        self.a1v2 = a1v2 # 1 = Arterial, 2 = Venous
        
        self.ev1 = ev1
        self.ev2 = ev2
        self.ev3 = ev3
        self.ecl = ecl
        self.eq2 = eq2
        self.eq3 = eq3
        
        # PMA is explicitly requested or derived
        if pma is not None:
            self.pma = pma
        else:
            self.pma = self.age * (52.1429/52.0) + 40/52.0 # Approximation

    def get_parameters(self):
        # Extract direct THETA values as strictly solved by NONMEM
        theta = {
            1: 1.837860e+00,
            2: 3.238730e+00,
            3: 5.608800e+00, # Note: THETA(3) in file is actually THETA(4) for v3ref? Let's check the index.
            # Wait, line 70 is (4, 5.608, 7) for V3. But earlier it uses THETA(3). I must align with the NONMEM block exactly.
            # I will manually align the THETAs based on the $THETA block
        }
        
        theta1 = 1.837860e+00
        theta2 = 3.238730e+00
        theta3 = 5.608800e+00 # v3ref
        theta4 = 5.819830e-01 # clref(male)
        theta5 = 5.596720e-01 # q2ref
        theta6 = 1.030460e-01 # q3ref
        theta7 = 1.913070e-01 # log-error
        theta8 = 3.744220e+00 # maturation CL E50
        theta9 = 2.203300e+00 # maturation CL slope
        theta10 = -1.563300e-02 # V2 declines with age
        theta11 = -2.857090e-03 # CL declines with age
        theta12 = 3.513130e+00 # V1 sigmoid E50
        theta13 = -1.381660e-02 # V3 declines with age
        theta14 = 4.223570e+00 # maturation Q3 E50
        theta15 = 7.420430e-01 # clref(female)
        theta16 = 2.656420e-01 # q2ref(child) factor
        theta17 = 3.498850e-01 # v1 extra venous
        theta18 = -3.849270e-01 # q2 less for venous

        # ETA values are effectively 0 for the generic population profile, but NONMEM datasets include Bayesian estimates explicitly.
        # However, for a generalized PK profile (like we are drawing via simulation), ETAs = 0.
        eta1 = eta2 = eta3 = eta4 = eta5 = eta6 = eta7 = 0.0

        # Al-sallami FFM
        ht2 = (self.hgt / 100.0) ** 2
        matm = 0.88 + ((1 - 0.88) / (1 + (self.age / 13.4)**(-12.7)))
        matf = 1.11 + ((1 - 1.11) / (1 + (self.age / 7.1)**(-1.1)))
        matr = 0.88 + ((1 - 0.88) / (1 + (35.0 / 13.4)**(-12.7)))
        
        ffmm = matm * 42.92 * ht2 * self.wgt / (30.93 * ht2 + self.wgt)
        ffmf = matf * 37.99 * ht2 * self.wgt / (35.98 * ht2 + self.wgt)
        ffmr = matr * 42.92 * (1.7 * 1.7) * 70.0 / (30.93 * (1.7 * 1.7) + 70.0)
        
        ffm = ffmm * (2 - self.m1f2) + ffmf * (self.m1f2 - 1)
        nffm = ffm / ffmr

        # maturation
        dv1 = 1
        dv2 = 1
        dv3 = 1
        
        # sigmoidal maturation of CL
        pmw = self.pma * 52.0
        pmr = (35.0 + 40.0 / 52.0) * 52.0
        me50 = math.exp(theta8)
        mgam = math.exp(theta9)
        mcl = (pmw**mgam) / (pmw**mgam + me50**mgam)
        rcl = (pmr**mgam) / (pmr**mgam + me50**mgam)
        dcl = mcl / rcl
        dq2 = 1

        # sigmoidal maturation of Q3 based on 40 weeks gestation
        pmew = self.age * 52.0 + 40.0
        pmer = 35.0 * 52.0 + 40.0
        qe50 = math.exp(theta14)
        mq3 = pmew / (pmew + qe50)
        rq3 = pmer / (pmer + qe50)
        dq3 = mq3 / rq3

        # aging
        kv1 = 1
        kv2 = math.exp(theta10 * (self.age - 35.0))
        kv3 = math.exp(theta13 * self.age * (self.tech - 1))
        kcl = math.exp(theta11 * self.age * (self.tech - 1))
        kq2 = 1
        kq3 = 1

        # covariate structure
        vv50 = math.exp(theta12)
        cv1 = self.wgt / (self.wgt + vv50)
        rv1 = 70.0 / (70.0 + vv50)
        
        m1 = (cv1 / rv1) * kv1 * dv1
        vcv1 = (self.a1v2 - 1) * (1 - cv1)
        v1 = math.exp(theta1 + eta1) * m1 * (1 + vcv1 * math.exp(theta17))
        
        m2 = (self.wgt / 70.0)**1 * kv2 * dv2
        v2 = math.exp(theta2 + eta2) * m2
        
        m3 = (nffm)**1 * kv3 * dv3
        v3 = math.exp(theta3 + eta3) * m3
        
        m4 = (self.wgt / 70.0)**0.75 * kcl * dcl
        cl_param = math.exp((2 - self.m1f2) * theta4 + (self.m1f2 - 1) * theta15 + eta4) * m4
        
        rv2 = math.exp(theta2)
        m5 = (v2 / rv2)**0.75 * kq2 * dq2
        km5 = 1 + math.exp(theta16) * (1 - mq3)
        q2 = math.exp(theta5 + eta5 + (self.a1v2 - 1) * theta18) * m5 * km5
        
        rv3 = math.exp(theta3)
        m6 = (v3 / rv3)**0.75 * kq3 * dq3
        q3 = math.exp(theta6 + eta6) * m6

        import pandas as pd
        if self.ev1 is not None and pd.notnull(self.ev1):
            v1 = self.ev1
            v2 = self.ev2
            v3 = self.ev3
            cl_param = self.ecl
            q2 = self.eq2
            q3 = self.eq3

        k10 = cl_param / v1
        k12 = q2 / v1
        k21 = q2 / v2
        k13 = q3 / v1
        k31 = q3 / v3
        ke0 = 0.146 * math.pow(self.wgt / 70.0, -0.25) # Standard Eleveld KE0

        return {
            'vc': v1,
            'k10': k10,
            'k12': k12,
            'k13': k13,
            'k21': k21,
            'k31': k31,
            'ke0': ke0
        }
