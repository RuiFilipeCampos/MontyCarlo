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

	def __init__(self, upper_dir, CS, F, R, I, E):
		self.upper_dir = upper_dir
		print("> *** photon/coherent: creating TotalCrossSection")
		self.CS   = TotalCrossSection(upper_dir, CS)
		self.CS.final_init(upper_dir.density)


		print("> *** photon/coherent: creating FormFactor")
		self.FF   = FormFactor(upper_dir, F)
		
		print("> *** photon/coherent/formfactor: creating RationalInterpolation")
		self.FF = RationalInterpolation(self.FF._eval, 0, 300**2, True)

	#finalize material
	def final_init(self, density):
		#RITA method on squared form factor
		self.FF = RationalInterpolation(self.FF._eval, 0, 300**2, True)
		self.CS.final_init(density)

	#New instances are given to parent photon.
	def __add__(self, other):
		new_self = Coherent.__new__(Coherent)
		new_self.CS  = self.CS + other.CS
		new_self.FF  = self.FF + other.FF

		return new_self
	
	def __mul__(self, other):
		new_self = Coherent.__new__(Coherent)
		new_self.CS = self.CS*other
		new_self.FF = self.FF*other
		return new_self

	#other stuff
	def __repr__(self):
		return "<" + self.upper_dir._name + ": /photon/coherent>"

	def showFrame(self, windows):
		#from matplotlib.figure import Figure
		
		w = tk.Frame(windows)

		#fig = Figure(seize=(6, 6))
		#subP = fig.add_subplot(111)
		#subP.plot(self.)


class FormFactor:
	"""Holds the form factor."""

	
	def __init__(self, upper_dir, CS):
		"""
		1) Checks if Form Factor fit has already been calculated.
		2) if there is no pickled fit, import FormFactorWriter
		3) proceed to read the result from fit
		4) using paramaters, create a callable representing |F(Z, x**2)|**2

		The callable is jitted (@njit), meaning that when the first call occurs
		the functions code will be compiled to machine code. First call
		is slow, but subsquent calls are fast. 
		"""

		self.upper_dir = upper_dir
		self.Z = upper_dir.Z

		import os
		path = directory + "\pickles\\" + str(upper_dir.Z)
		if not os.path.isfile(path):
			print("> *** photon/coherent: Writing form factor...")
			from .FormFactorWriter import FormFactorWriter
			FormFactorWriter(CS, upper_dir.Z)
			
		import dill as pickle
		with open(path, 'rb') as file:
			self.param = pickle.load(file)
			print(f"> *** photon/coherent: Unpickled paramaters of form factor of element {self.Z}")

		a1, a2, a3, a4, a5 = self.param
		Z = self.Z

		#depending on Z, the fit is calculated differently
		if Z < 11:

			@njit
			def _eval(x2):
				x = sqrt(x2)
				A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
				B = 1 + a4*x**2 + a5 * x**4
				f = Z*A/B**2
				return f**2

		else:
			@njit
			def _eval(x2):

				x = sqrt(x2)
				A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
				B = 1 + a4*x**2 + a5 * x**4
				f = Z*A/B**2

				if f > 2: return f**2
			
				alpha = 1/137.03599908421  #fine structure constant
				a = (Z - 5/6)*alpha

				Q = x/20.6074/2/a #changing variables to suit formula in paper

				gamma = (1 - a**2)**.5

				#calculating Fk
				num = sin(2*gamma*arctan(Q))
				den = gamma*Q*(1 + Q**2)**gamma
				Fk  = num/den

				return max(f, Fk)**2

		#whatever the result of fit, it's stored here
		self._eval = _eval

	def __repr__(self):
		return "<" + self.upper_dir._name + ": /photon/coherent/FormFactor>"

	def __add__(self, other):

		#constructing chain of operation
		eval1, eval2 = self._eval, other._eval

		@njit
		def new_fit(x2):
			return eval1(x2) + eval2(x2)

		#constructing new instance
		new_self = FormFactor.__new__(FormFactor)
		new_self._eval = new_fit
		return new_self



	def __mul__(self, other):
		#constructing chain of operation
		eval1 = self._eval

		@njit
		def new_fit(x2):
			return other*eval1(x2)

		#constructing new instance
		new_self = FormFactor.__new__(FormFactor)
		new_self._eval = new_fit
		return new_self

		




