from abc import ABC, abstractmethod

class PropofolModel(ABC):
    def __init__(self, weight, age, height, gender, is_adj_bw=False):
        self.weight = weight
        self.age = age
        self.height = height
        self.gender = gender # 0 for male, 1 for female (matching SimTIVA JS)
        self.is_adj_bw = is_adj_bw

        # derived standard variables matching SimTIVA
        self.lbm = self._calculate_lbm()
        self.adj_bw = self._calculate_adj_bw()

    def _calculate_lbm(self):
        # Default simple James equation (matching typical Schnider LBM fallback)
        if self.gender == 0:
            return 1.1 * self.weight - 128 * (self.weight / self.height) ** 2
        else:
            return 1.07 * self.weight - 148 * (self.weight / self.height) ** 2

    def _calculate_adj_bw(self):
        # Ideal body weight calculation (Robinson's approximation roughly)
        if self.gender == 0:
            ibw = 50 + 0.9 * (self.height - 152)
        else:
            ibw = 45.5 + 0.9 * (self.height - 152)
        return ibw + 0.4 * (self.weight - ibw)

    @abstractmethod
    def get_parameters(self):
        """
        Calculates and returns the specific volumes and clearance micro-rates
        for the given model based on the instance covariates.

        Returns a dictionary:
        {
            'vc': float,
            'k10': float,
            'k12': float,
            'k13': float,
            'k21': float,
            'k31': float,
            'ke0': float
        }
        """
        pass
