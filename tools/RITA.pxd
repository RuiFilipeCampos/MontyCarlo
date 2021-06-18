



from .integration cimport Integrator


cdef class RationalInterpolation:
	"""
	Rational Interpolator designed to implement the RITA (Rational Inverse
	Transform with Aliasing) method. 

	Receives a low level callable (custom C type) representing
	an unnormalized probability distribution, as well as its domain.
	"""
	cdef double norm, error
	cdef public list _intervals,xAxis, PAxis
	cdef Integrator INTEGRATOR
	cdef object p
	cdef double prob(self, double x)
	cpdef invCum(self, double r)
	cpdef cumul(self, double x)


cdef class Interval(RationalInterpolation):
	"""
	Interval [u, v[
	"""
	cdef public double u,  v, Pu, Pv, pu, pv, dP, dx, a, b, h
	cdef RationalInterpolation upper
	@staticmethod
	cdef public Interval _create(RationalInterpolation upper,
								 double u, 
								 double v, 
								 double pu, 
								 double Pu)
	cpdef public double _eval(self, double x)
	cpdef double _invCum(self, double r)
	cdef double f(self, double k)
	cdef double _calculateError(self)
	cdef tuple _split(self)
	cdef double _cumul(self)
