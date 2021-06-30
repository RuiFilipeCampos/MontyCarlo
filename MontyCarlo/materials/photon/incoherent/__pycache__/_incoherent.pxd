from ...materials   cimport Material
from ....tools.RITA cimport RationalInterpolation, func


cdef class Incoherent:
	cdef public object CS
	cdef public IncoherentFormFactor  S


cdef class CompoundIncoherent(Incoherent):
	pass




cdef class IncoherentFormFactor(func):
	"""Holds the form factor."""
	cdef Material upper_dir
	cdef tuple param
	cdef double a1, a2, a3, a4, a5, Z

	#methods
	cdef double _eval(self, double x)

cdef class CompoundIncoherentFF(IncoherentFormFactor):
	cdef IncoherentFormFactor S1, S2


cdef class ScaledIncoherentFF(IncoherentFormFactor):
	cdef IncoherentFormFactor S1
	cdef double scalar

