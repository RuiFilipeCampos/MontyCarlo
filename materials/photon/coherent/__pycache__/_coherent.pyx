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
from ....tools.distributions import UnivariateDistribution
from ....tools.data import getAxis
from ....tools.RITA cimport RationalInterpolation, func

#built in imports
import os


#dealing with paths >.<
from ....settings import __montecarlo__
__coherent__ = __montecarlo__/'materials'/'photon'/'coherent'

directory = str(__coherent__)















cdef class Coherent:
	def __init__(self, 
				 Material upper_dir, 
				 tuple CS, 
			  	 tuple F, 
				 tuple R, 
				 tuple I, 
				 tuple E):


		self.upper_dir = upper_dir
		self.CS  = TotalCrossSection(upper_dir, CS)
		self._FF  = FormFactor(upper_dir, F)

	cdef public double thomsonDCS(self, double c):
		"""Thomson DCS."""
		return (1 + c**2)/2

	def __repr__(self):
		return "<" + self.upper_dir._name + ": /photon/coherent>"
 

	cpdef final_init(self, double density):
		#RITA method on squared form factor
		self.FF = RationalInterpolation.create(self._FF, 0, 300**2, True)
		self.CS.final_init(density)

	def __add__(self, Coherent other):
		CS  = self.CS + other.CS
		_FF  = self._FF + other._FF
		return CompoundCoherent(self.upper_dir, CS, _FF)
	
	def __mul__(self, float other):
		CS = self.CS*other
		_FF = self._FF*other
		return CompoundCoherent(self.upper_dir, CS, _FF)

		

cdef class CompoundCoherent(Coherent):
	def __init__(self, Material upper_dir, object CS, FormFactor _FF):
		self.CS = CS
		self._FF = _FF

























cdef class FormFactor(func):
	"""Holds the form factor."""

	
	def __init__(self, Material upper_dir, tuple CS):

		self.upper_dir = upper_dir
		self.Z = upper_dir.Z

		import os
		path = directory + "\pickles\\" + str(upper_dir.Z)
		if not os.path.isfile(path):
			print("> Writing form factor...")
			from .FormFactorWriter import FormFactorWriter
			FormFactorWriter(CS, upper_dir.Z)
			
		import dill as pickle
		with open(path, 'rb') as file:
			self.param = pickle.load(file)

		self.a1, self.a2, self.a3, self.a4, self.a5 = self.param



	
	cdef double _eval(self, double x2):
		"""
		Callable. Receives x**2, outputs |F(Z, x)|**2.
		"""
		cdef double A, B, f, alpha, a, Q, gamma, num, den, Fk, x
		

		x = sqrt(x2)
		A = 1 + self.a1*x**2 + self.a2 * x**3 + self.a3 * x**4
		B = 1 + self.a4*x**2 + self.a5 * x**4
		f = self.Z*A/B**2

		if self.Z < 11: return f**2
		if f < 2: return f**2
		else:
			alpha = 1/137.03599908421  #fine structure constant
			a = (Z - 5/6)*alpha

			Q = x/20.6074/2/a #changing variables to suit formula in paper

			gamma = (1 - a**2)**.5

			#calculating Fk
			num = sin(2*gamma*arctan(Q))
			den = gamma*Q*(1 + Q**2)**gamma
			Fk  = num/denom

			return max(f, Fk)**2

		
	def __call__(self, x):
		return self._eval(x)

	def __repr__(self):
		return "<" + self.upper_dir._name + ": /photon/coherent/FormFactor>"

		
	def __add__(self, other):
		return CompoundFormFactor(self, other)

	def __mul__(self, other):
		return ScaledFormFactor(self, other)




cdef class CompoundFormFactor(FormFactor):
	"""
	Container class. Holds two squared form factors, and outputs
	their sum when called.
	"""
	def __init__(self, FormFactor F1, FormFactor F2):
		self.F1 = F1
		self.F2 = F2

	def __call__(self, double x):
		return self.F1(x) + self.F2(x)

cdef class ScaledFormFactor(FormFactor):
	"""
	Container class. Holds one squared form factor and a scalar, outputs
	their product when called.
	"""
	def __init__(self, FormFactor F1, double scalar):
		self.scalar = scalar
		self.F1 = F1

	def __call__(self, double x):
		return self.scalar*self.F1(x)


