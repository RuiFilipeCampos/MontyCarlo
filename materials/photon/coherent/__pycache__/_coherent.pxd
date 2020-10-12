
from ...materials   cimport Material
from ....tools.RITA cimport RationalInterpolation, func

cdef class Coherent:
	cdef public FormFactor _FF
	cdef public RationalInterpolation FF

	cdef public object CS
	cdef public Material upper_dir
	cdef public double thomsonDCS(self, double c)
	cpdef final_init(self, double density)

cdef class CompoundCoherent(Coherent):
	cdef bint x
















cdef class FormFactor(func):
	"""Holds the form factor."""

	cdef tuple param
	cdef Material upper_dir
	cdef double a1, a2, a3, a4, a5, Z
	cdef double _eval(self, double x)




cdef class CompoundFormFactor(FormFactor):
	"""
	Container class. Holds two squared form factors, and outputs
	their sum when called.
	"""
	cdef FormFactor F1, F2


cdef class ScaledFormFactor(FormFactor):
	"""
	Container class. Holds one squared form factor and a scalar, outputs
	their product when called.
	"""
	cdef double scalar
	cdef FormFactor F1




