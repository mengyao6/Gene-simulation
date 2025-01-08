import numpy as np

class Genes:
    def __init__(self, values=None):
        # Initialize gene activity values to 0 as np.float64
        self.cox = np.float64(0)
        self.narG = np.float64(0)
        self.nirk = np.float64(0)
        self.nrf = np.float64(0)
        self.dsr = np.float64(0)
        self.amoA = np.float64(0)
        self.hzo = np.float64(0)
        self.nap = np.float64(0)
        self.nor = np.float64(0)
        self.sox = np.float64(0)

        if values is not None:
            self.update(values)

    def update(self, new_values):
        """Update gene activity values."""
        # 0: cox, 1: narG, 2: nirk, 3: nrf, 4: dsr, 5: amoA, 6: hzo, 7: nap, 8: nor, 9: sox
        (self.cox, self.narG, self.nirk, self.nrf, 
         self.dsr, self.amoA, self.hzo, self.nap, 
         self.nor, self.sox) = new_values

    def to_array(self):
        """Convert gene activity values to a numpy array."""
        return np.array([self.cox, self.narG, self.nirk, self.nrf, 
                         self.dsr, self.amoA, self.hzo, self.nap, 
                         self.nor, self.sox], dtype=np.float64)


class Mus:
    def __init__(self):
        # Initialize rate constants (unit: s^-1)
        self.cox = 0.28 / 86400
        self.narG = 0.151 / 86400
        self.nirk = 0.247 / 86400
        self.nrf = 0.162 / 86400
        self.dsr = 0.0636 / 86400
        self.amoA = 0.432 / 86400
        self.hzo = 0.864 / 86400
        self.nap = 0.864 / 86400
        self.nor = 0.432 / 86400
        self.sox = 0.864 / 86400

    def to_array(self):
        """Convert the rate constants to a NumPy array."""
        return np.array([
            self.cox, self.narG, self.nirk, self.nrf,
            self.dsr, self.amoA, self.hzo, self.nap,
            self.nor, self.sox
        ], dtype=np.float64)