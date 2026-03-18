import math
from .base_model import PropofolModel

class PaedfusorModel(PropofolModel):
    def get_parameters(self):
        k12 = 0.114
        k13 = 0.0419
        k21 = 0.055
        k31 = 0.0033

        age = self.age
        mass = self.weight

        if age < 13:
            vc = 0.4584 * mass
            k10 = 0.1527 * math.pow(mass, -0.3)
            var_const = 0.7751 * math.exp(-0.099 * age)
            ke0 = 0.0351 * math.log(mass) + var_const
        elif 13 <= age < 14:
            vc = 0.4 * mass
            k10 = 0.0678
            ke0 = 0.319
        elif 14 <= age < 15:
            vc = 0.342 * mass
            k10 = 0.0792
            ke0 = 0.286
        elif 15 <= age < 16:
            vc = 0.284 * mass
            k10 = 0.0954
            ke0 = 0.251
        else: # Adult equivalent
            vc = 0.229 * mass
            k10 = 0.119
            ke0 = 1.21  # Modified Marsh ke0

        return {
            'vc': vc,
            'k10': k10,
            'k12': k12,
            'k13': k13,
            'k21': k21,
            'k31': k31,
            'ke0': ke0
        }
