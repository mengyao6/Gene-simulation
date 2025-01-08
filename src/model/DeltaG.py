import numpy as np


class DeltaGs:
    def __init__(self, value=None):
        # Default concentration values as np.float64(0)
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

    def to_array(self):
        return np.array([self.cox, self.narG, self.nirk, self.nrf, 
                        self.dsr, self.amoA, self.hzo, self.nap, 
                        self.nor, self.sox])


