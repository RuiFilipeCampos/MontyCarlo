# cython: profile=True

#External Imports
from numpy import *
from scipy.interpolate import *
from scipy.integrate import *

#Internal Imports
#from .photon import Photon
from ....tools.CStools import TotalCrossSection
from ....tools.RITA cimport RationalInterpolation, func
from ....tools.data import getAxis

#dealing with paths >.<
from ....settings import __montecarlo__
__incoherent__ = __montecarlo__/'materials'/'photon'/'incoherent'

directory = str(__incoherent__)



cdef class Incoherent:
	def __init__(self, Material upper_dir, tuple CS, tuple S):
		self.CS  = TotalCrossSection(upper_dir, CS)
		self.S  = IncoherentFormFactor(upper_dir, S)

	def __add__(self, other):
		CS = self.CS + other.CS
		S  = self.S + other.S
		return CompoundIncoherent(CS, S)


	def __mul__(self, other):
		CS = self.CS*other
		S  = self.S*other
		return CompoundIncoherent(CS, S)

	def final_init(self, density):
		self.CS.final_init(density)
		#self.S = RationalInterpolation.create(self._S, 1e-5, 300, False) #gotta check this


cdef class CompoundIncoherent(Incoherent):
	def __init__(CS, S):
		self.CS = CS
		self.S = S






cdef class IncoherentFormFactor(func):
	"""Holds the form factor."""


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
		self.a1, self.a2, self.a3, self.a4, self.a5 = self.param


	cdef double _eval(self, double x):

		cdef double A, B


		A = 1 + self.a1*x**2 + self.a2 * x**3 + self.a3 * x**4
		B = 1 + self.a4*x**2 + self.a5 * x**4

		return 1 - A/B**2
			
	
	def __call__(self, x):
		return self._eval(x)

	def __add__(self, IncoherentFormFactor other):
		return CompoundIncoherentFF(self, other)

	def __mul__(self, double other):
		return ScaledIncoherentFF(self, other)

cdef class CompoundIncoherentFF(IncoherentFormFactor):
	def __init__(self, S1, S2):
		self.S1 = S1
		self.S2 = S2

	def __call__(self, double x):
		return self.S1(x) + self.S2(x)

cdef class ScaledIncoherentFF(IncoherentFormFactor):
	def __init__(self, S1, scalar):
		self.S1 = S1
		self.scalar = scalar

	def __call__(self, double x):
		return self.scalar*self.S1(x)

