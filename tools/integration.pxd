cdef class Integrator:
	#attributes
	cdef object f
	cdef public long double tol
	cdef public tuple unstable_intervals, stable_intervals
	cdef public long double result, error

	#methods
	#@staticmethod
	#cdef public Integrator create(object f, long double x0, long double xf)
	cdef long double sumResult(self)
	cdef splitAll(self)



cdef class Interval(Integrator):
	#attributes
	cdef public bint stable
	cdef public long double x0, xf

	#methods
	@staticmethod
	cdef public Interval create(object f, long double x0, long double xf)
	cdef public tuple split(self)
	cdef long double integrate(self, long double a, long double b)