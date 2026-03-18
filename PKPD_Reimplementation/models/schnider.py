from .base_model import PropofolModel

class SchniderModel(PropofolModel):
    def get_parameters(self):
        vc = 4.27
        v2 = 18.9 - 0.391 * (self.age - 53)
        v3 = 238

        if not self.is_adj_bw:
            cl1 = 1.89 + 0.0456 * (self.weight - 77) - 0.0681 * (self.lbm - 59) + 0.0264 * (self.height - 177)
        else:
            if self.gender == 0:
                lbm2 = 1.1 * self.adj_bw - 128 * (self.adj_bw / self.height) ** 2
            else:
                lbm2 = 1.07 * self.adj_bw - 148 * (self.adj_bw / self.height) ** 2
            cl1 = 1.89 + 0.0456 * (self.adj_bw - 77) - 0.0681 * (lbm2 - 59) + 0.0264 * (self.height - 177)

        cl2 = 1.29 - 0.024 * (self.age - 53)
        cl3 = 0.836

        k10 = cl1 / vc
        k12 = cl2 / vc
        k13 = cl3 / vc
        k21 = cl2 / v2
        k31 = cl3 / v3
        ke0 = 0.456

        return {
            'vc': vc,
            'k10': k10,
            'k12': k12,
            'k13': k13,
            'k21': k21,
            'k31': k31,
            'ke0': ke0
        }
