from .base_model import PropofolModel

class MarshModel(PropofolModel):
    def get_parameters(self):
        # Check for Adjusted Body Weight usage mapping SimTIVA behavior 
        if self.is_adj_bw:
            vc = 0.228 * self.adj_bw
        else:
            vc = 0.228 * self.weight

        # Pure Marsh constants
        k10 = 0.119
        k12 = 0.112
        k13 = 0.0419
        k21 = 0.055
        k31 = 0.0033
        ke0 = 1.21  # 'fast' ke0 

        return {
            'vc': vc,
            'k10': k10,
            'k12': k12,
            'k13': k13,
            'k21': k21,
            'k31': k31,
            'ke0': ke0
        }
