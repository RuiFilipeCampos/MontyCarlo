#External Imports
from numpy import *
from scipy.interpolate import *
from scipy.integrate import *

#Internal Imports
#from .photon import Photon
from ....tools.CStools import TotalCrossSection
from ....tools.distributions import UnivariateDistribution
from ....tools.data import getAxis



class Pairproduction:
    
    def __init__(self, upper_dir, CS):
        self.CS = TotalCrossSection(upper_dir, CS)
        self.CS.final_init(upper_dir.density)
        #self.S  = FormFactor(F)

    def __add__(self, other):
        self.CS = self.CS + other.CS
        #self.F  = self.F + other.F
        return self
    
    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        self.CS = self.CS*other
        #self.F  = self.F*other
        return self
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def final_init(self, density):
        self.CS.final_init(density)
