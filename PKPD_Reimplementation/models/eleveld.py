import math
from .base_model import PropofolModel

class EleveldModel(PropofolModel):

    def __init__(self, weight, age, height, gender, is_adj_bw=False, pma=None, opioid=True, ev1=None, ev2=None, ev3=None, ecl=None, eq2=None, eq3=None):
        super().__init__(weight, age, height, gender, is_adj_bw)
        self.opioid = opioid
        # Eleveld calculates age structurally using week scale for infants
        self.toweeks = 52.1429
        
        self.ev1 = ev1
        self.ev2 = ev2
        self.ev3 = ev3
        self.ecl = ecl
        self.eq2 = eq2
        self.eq3 = eq3
        
        # Post Menstrual Age - either explicitly passed or fallback to generic
        if pma is not None:
            self.pma = pma
        else:
            self.pma = self.age * self.toweeks + 40

    def _fsigmoid(self, x, y, z):
        return math.pow(x, z) / (math.pow(x, z) + math.pow(y, z))

    def _fcentral(self, x):
        return self._fsigmoid(x, 33.6, 1)

    def _fageing(self, x):
        return math.exp(x * (self.age - 35))

    def _fclmaturation(self, x):
        return self._fsigmoid(x, 42.3, 9.06)

    def _fq3maturation(self, x):
        return self._fsigmoid(x + 40, 68.3, 1)

    def _fffm(self):
        bmi = self.weight / math.pow(self.height / 100.0, 2)
        if self.gender == 0:
            factor = (0.88 + (1 - 0.88) / (1 + math.pow((self.age / 13.4), -12.7)))
            return factor * ((9270 * self.weight) / (6680 + 216 * bmi))
        else:
            factor = (1.11 + (1 - 1.11) / (1 + math.pow((self.age / 7.1), -1.1)))
            return factor * ((9270 * self.weight) / (8780 + 244 * bmi))

    def get_parameters(self):
        vc = 6.28 * self._fcentral(self.weight) / self._fcentral(70)
        v2 = 25.5 * (self.weight / 70) * self._fageing(-0.0156)
        v2ref = 25.5
        
        # Fixed reference point used mathematically in eleveld
        ffmref = (0.88 + (1 - 0.88) / (1 + math.pow((35 / 13.4), -12.7))) * ((9270 * 70) / (6680 + 216 * 24.22145))

        if self.opioid:
            v3 = 273 * self._fffm() / ffmref * math.exp(-0.0138 * self.age)
        else:
            v3 = 273 * self._fffm() / ffmref
        v3ref = 273

        # Maturation variables
        mat_factor = self._fclmaturation(self.pma) / self._fclmaturation(35 * self.toweeks + 40)
        age_decay = math.exp(-0.00286 * self.age)

        if self.gender == 0:
            if self.opioid:
                cl1 = 1.79 * math.pow(self.weight / 70, 0.75) * mat_factor * age_decay
            else:
                cl1 = 1.79 * math.pow(self.weight / 70, 0.75) * mat_factor
        else:
            if self.opioid:
                cl1 = 2.1 * math.pow(self.weight / 70, 0.75) * mat_factor * age_decay
            else:
                cl1 = 2.1 * math.pow(self.weight / 70, 0.75) * mat_factor

        cl2 = 1.75 * math.pow(v2 / v2ref, 0.75) * (1 + 1.3 * (1 - self._fq3maturation(self.age * self.toweeks)))
        cl3 = 1.11 * math.pow(v3 / v3ref, 0.75) * (self._fq3maturation(self.age * self.toweeks) / self._fq3maturation(35 * self.toweeks))

        import pandas as pd
        if hasattr(self, 'ev1') and self.ev1 is not None and pd.notnull(self.ev1):
            vc = self.ev1
            v2 = self.ev2
            v3 = self.ev3
            cl1 = self.ecl
            cl2 = self.eq2
            cl3 = self.eq3

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
            'ke0': ke0
        }
