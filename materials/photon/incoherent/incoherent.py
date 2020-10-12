# cython: profile=True

#External Imports
from numpy import *
from scipy.interpolate import *
from scipy.integrate import *
from numba import njit

#Internal Imports
#from .photon import Photon
from ....tools.CStools import TotalCrossSection
from ....tools.RITA import RationalInterpolation
from ....tools.data import getAxis

#dealing with paths >.<
from ....settings import __montecarlo__
__incoherent__ = __montecarlo__/'materials'/'photon'/'incoherent'

directory = str(__incoherent__)



class Incoherent:

	# INIT PHASE
	############################################################################
	def __init__(self, upper_dir, CS, S):
		print("> *** photon/incoherent: creating TotalCrossSection")
		self.CS  = TotalCrossSection(upper_dir, CS)
		self.CS.final_init(upper_dir.density)

		print("> *** photon/incoherent: creating IncoherentFormFactor")
		self.S   = IncoherentFormFactor(upper_dir, S)

	############################################################################






	# CONSTRUCTION PHASE
	############################################################################
	#New instances are given to parent photon.
	def __add__(self, other):
		new_self = Incoherent.__new__(Incoherent)
		new_self.CS  = self.CS + other.CS
		new_self.S  = self.S + other.S
		return new_self
	
	def __mul__(self, other):
		new_self = Incoherent.__new__(Incoherent)
		new_self.CS = self.CS*other
		new_self.S = self.S*other
		return new_self
	############################################################################




	# FINALIZE 
	############################################################################
	#finalize material
	def final_init(self, density):
		self.CS.final_init(density)
		#self.S = RationalInterpolation.create(self._S, 1e-5, 300, False) #gotta check this
	############################################################################















class IncoherentFormFactor:
	"""Holds the form factor."""

	# INIT PHASE
	############################################################################
	def __init__(self, upper_dir, S):
		self.upper_dir = upper_dir
		self.Z = upper_dir.Z

		import os

		path = directory + "\pickles\\" + str(upper_dir.Z)
		
		if not os.path.isfile(path):
			print("> Writing Incoherent Form Factor...")
			from .IncoherentFormFactorWriter import IncoherentFormFactorWriter
			IncoherentFormFactorWriter(S, upper_dir.Z)
			
		import dill as pickle
		with open(path, 'rb') as file:
			self.param = pickle.load(file)

		a1, a2, a3, a4, a5 = self.param
		Z = self.Z

		@njit
		def _eval(x):
			A = 1 + a1*x**2 + a2 * x**3 + a3 * x**4
			B = 1 + a4*x**2 + a5 * x**4
			return (1 - A/B**2)/Z

		self._eval = _eval
	############################################################################












	############################################################################
	#COMPOSITION PHASE
	#I am not sure about the rules for constructing the compound incoherent
	#form factor.
	############################################################################
	def __add__(self, other):
		eval1 = self._eval
		eval2 = other._eval

		@njit
		def new_eval(x):
			return eval1(x) + eval2(x)

		new_self = IncoherentFormFactor.__new__(IncoherentFormFactor)
		new_self._eval = new_eval

		return new_self
	
	def __mul__(self, other):
		eval_ = self._eval

		@njit
		def new_eval(x):
			return other*eval_(x)
		
		new_self = IncoherentFormFactor.__new__(IncoherentFormFactor)
		new_self._eval = new_eval
		return new_self
	############################################################################
	############################################################################









	#RUNTIME PHASE
	############################################################################
	def __call__(self, x):
		return self._eval(x)
	############################################################################
	############################################################################
	############################################################################