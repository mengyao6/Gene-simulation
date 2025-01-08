import numpy as np


class Chemicals:
    def __init__(self, values=None):
        # Default concentration values as np.float64(0)
        self.C_co2 = np.float64(0)
        self.C_hco3 = np.float64(0)
        self.C_c6 = np.float64(0)
        self.C_H = np.float64(0)
        self.C_no2 = np.float64(0)
        self.C_no3 = np.float64(0)
        self.C_n2 = np.float64(0)
        self.C_nh4 = np.float64(0)
        self.C_h2s = np.float64(0)
        self.C_so4 = np.float64(0)
        self.C_DIC = np.float64(0)
        self.C_o2 = np.float64(0)
        self.C_POC = np.float64(0)
        
        if values is not None:
            self.update(values)

    def update(self, new_values):
        """Update concentration values."""
        # 0: C_co2, 1: C_hco3, 2: C_c6, 3: C_H, 4: C_no2, 5: C_no3, 6: C_n2, 7: C_nh4, 8: C_h2s, 9: C_so4, 10:C_DIC,  11:C_o2 ,    12:C_POC
        (self.C_co2, self.C_hco3, self.C_c6, self.C_H, self.C_no2, 
         self.C_no3, self.C_n2, self.C_nh4, self.C_h2s, self.C_so4, 
         self.C_DIC, self.C_o2, self.C_POC) = new_values

    def to_array(self):
        """Convert concentration values to a numpy array."""
        return np.array([self.C_co2, self.C_hco3, self.C_c6, self.C_H, 
                         self.C_no2, self.C_no3, self.C_n2, self.C_nh4, 
                         self.C_h2s, self.C_so4, self.C_DIC, self.C_o2, self.C_POC], dtype=np.float64)