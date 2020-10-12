cimport numpy as cnp

import numpy as np

#internal imports
from integration cimport Integrator

ctypedef double (*func)(double)

cdef class RationalInterpolation:
	"""
	Rational Interpolator designed to implement the RITA (Rational Inverse
	Transform with Aliasing) method. 

	Receives a low level callable (custom C type) representing
	an unnormalized probability distribution, as well as its domain.
	"""
	cdef double norm, error
	cdef dict explain
	cdef list _intervals,xAxis, PAxis	
	cdef Integrator INTEGRATOR
	cdef func p

	@staticmethod
	cdef public RationalInterpolation create(func p, 
											 long double x0, 
											 long double xf):

		#Manualy creating the instance of RationalInterpolation
		self = <RationalInterpolation>RationalInterpolation.__new__(RationalInterpolation)
		
		#the 'INTEGRATOR' stuff will change, going to make it a function
		#instead of a class
		self.p = p
		self.INTEGRATOR = Integrator.create(p, x0, xf)
		self.norm = self.INTEGRATOR.result


		
		cdef double dx = (xf - x0)/11 
		cdef Interval I = Interval._create(self, 
										   0.,  	      #x0 or u
										   dx, 			  #xf or v
										   self.prob(0.), #prob value calculated in previous iteration
										   0.)			  #cumul value calculated in previous iteration

		
		self._intervals = [I]

		cdef double v = dx

		cdef int i
		for i in range(10):
			v += dx
			I = Interval._create(self, 
								 I.v,  #first value of interval is last value of previous interval 
								 v,    #new last value
								 I.pv, #prob  value calculated in previous iteration
								 I.Pv) #cumul value calculated in previous iteration
			
			self._intervals += [I]

		
		cdef Interval I1, I2
		cdef list errors = [I.error for I in self]
		cdef double err
		cdef int k

		for _ in range(1000):  
			err = max(errors)
			k   = errors.index(err)

			I = self._intervals[k]

			I1, I2 = I._split()

			self._intervals[k] = I1
			self._intervals.insert(k+1, I2)

			errors = [I.error for I in self]

		self.PAxis = [I.Pu for I in self._intervals] + [self._intervals[-1].Pv]
		self.xAxis = [I.u for I in self._intervals]  + [self._intervals[-1].v]
		return self


	def __iter__(self):
		yield from self._intervals

	cdef double prob(self, double x):
		return self.p(x)/self.norm

	def __call__(self, double x):
		k = np.searchsorted(self.xAxis, x)
		return self._intervals[k-1]._eval(x)

	def scatter(self):
		

		import matplotlib.pyplot as plt
		plt.figure()
		ax = plt.gca()

		#x_interp = arange(min(self.xAxis), max(self.xAxis), .01)
		y_interp = [self(x) for x in self.xAxis]

		ax.plot(self.xAxis, y_interp)
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		plt.show()








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
								 double Pu):


		self = <Interval>Interval.__new__(Interval)
		self.upper = upper
		self.u,  self.v  = u,  v
		self.Pu, self.pu = Pu, pu


		cdef Integrator INTEGRATOR = Integrator.create(upper.p, 0., v)
		self.Pv = INTEGRATOR.result/upper.norm

		self.dx = self.v -  self.u
		self.dP = self.Pv - self.Pu

		#calculate new value of cumul
		self.pv = upper.prob(v)

		#best, copy, paste, ever
		self.b = 1 - (self.dP/self.dx)**2 / self.pu / self.pv
		self.a = (self.dP/self.dx) /self.pu - self.b - 1

		self.h = self.dx/50
		self.error = self._calculateError()
		return self



	cpdef public double _eval(self, double x):
		if self.a == -2. and self.b == 1.: return 0.
		if x == self.u:                    return self.pu
		if x == self.v:                    return self.pv #this should raise a value error, the interval is unbounded at v

		cdef double tau = (x - self.u)/self.dx

		cdef double B   = 4*self.b*tau**2
		cdef double A   = (1 + self.a + self.b - self.a*tau)/2/self.b/tau
		cdef double C   = (1 + self.a + self.b - self.a*tau)**2
		cdef double nu  = A * (1 - (1 - B/C)**.5)
		
		A = (1 + self.a*nu + self.b*nu**2)**2
		B = (1 + self.a + self.b)*(1 - self.b*nu**2)
		
		return A/B * self.dP/self.dx

	cpdef double _invCum(self, double r):
		cdef double nu = (r - self.Pu)/self.dP
		cdef double A  = (1 + self.a + self.b)*nu
		cdef double B  =  1 + self.a*nu + self.b*nu**2
		return self.u + A*(self.v - self.u)/B

	cdef double f(self, double k):
		cdef double x = self.u + k*self.h
		return (self.upper.prob(x) - self._eval(x))**2


	cdef double _calculateError(self):
		"""Estimate integral of prob(x)-interp(x) between u and v."""

		cdef double I1 = 0.
		for i in range(48):
			I1 += 4*self.f(i+1)

		cdef double I2 = 0.
		for i in range(46):
			I2 += 2*self.f(i+2)

		return self.h/50 * (self.f(0.) + 4.*I1 + 2.*I2 + self.f(50.))


	cdef tuple _split(self):

		cdef double mid_point = (self.u + self.v)/2

		cdef Interval I1 = Interval._create(self.upper, 
											self.u, 
											mid_point, 
											self.pu,
											self.Pu)

		cdef Interval I2 = Interval._create(self.upper, 
											mid_point, 
											self.v, 
											I1.pv, 
											I1.Pv)
		return I1, I2

	def __contains__(self, double x):
		return self.u <= x < self.v





cdef double f(double x):
	return x**2 + 1


print("start")
x = RationalInterpolation.create(f, 0, 1)
print("end")

x.scatter()


		