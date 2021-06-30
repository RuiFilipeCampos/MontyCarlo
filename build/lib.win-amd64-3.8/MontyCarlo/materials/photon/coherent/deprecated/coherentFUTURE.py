# cython: profile=True

"""
NOTES


  2) The angular distribution of coherently scattered photons is,  100 1451   55
																   100 1451   56
	 S(E,MU) = 3*T/8*(1 + MU**2)*((FF(M)+F1(E))**2+(F2(E))**2)     100 1451   57
																   100 1451   58
	 S(E,MU) = Angular distributions (barns per unit cosine)       100 1451   59
	 T       = The Thomson cross section                           100 1451   60
	 MU      = The cosine of the scattering angle                  100 1451   61
	 FF(M)   = The atomic form factor                              100 1451   62
	 M       = Momentum transfer = SIN(THETA/2)/LAMBDA/4*PI        100 1451   63
	 F1(E)   = The real anomalous scattering factor                100 1451   64
	 F2(E)   = The imaginary anomalous scattering factor           100 1451   65
																   100 1451   66
	 F1(E) and F2(E) are isotropic.                                100 1451   67
																   100 1451   68
	 When F1(E) = F2(E) = 0 this reduces to the result obtained    100 1451   69
	 using only form factors (the older ENDF/B convention).        100 1451   70

"""

#External Imports
from numpy import *
from scipy.interpolate import *
from scipy.integrate import *
from scipy.optimize import curve_fit

from numba import jit, cfunc, njit

#Internal Imports
from ....tools.CStools import TotalCrossSection

from ....tools.data import getAxis
from ....tools.RITA import RationalInterpolation

#built in imports
import os


#dealing with paths >.<
from ....settings import __montecarlo__
__coherent__ = __montecarlo__/'materials'/'photon'/'coherent'

directory = str(__coherent__)

class Coherent:
	def __init__(self, formula, density):
		self.CS   = TotalCrossSection(formula, density)
		self.FF   = FormFactor(formula, density)
    		#self.FF = RationalInterpolation(self.FF._eval, 0, 300**2, True)




