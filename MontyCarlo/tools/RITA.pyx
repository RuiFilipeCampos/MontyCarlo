__doc__ = """
"""

__author__ = "Rui Campos"


print("Importing .tools.RITA")



cimport numpy as cnp

import numpy as np

#internal imports
from .integration cimport Integrator


cdef class RationalInterpolation:
	"""
	Rational Interpolator designed to implement the RITA (Rational Inverse
	Transform with Aliasing) method. 

	Receives a low level callable (custom C type) representing
	an unnormalized probability distribution, as well as its domain.
	"""


	#@staticmethod
	def __init__(self, 
				object p, 
				long double x0,
				long double xf,
				bint normalize):

		#the 'INTEGRATOR' stuff will change, going to make it a function
		#instead of a class
		self.p = p
		if normalize:
			self.INTEGRATOR = Integrator.create(p, x0, xf)
			self.norm = self.INTEGRATOR.result
		else:
			self.norm = 1.


		
		cdef double dx = (xf - x0)/11 
		cdef Interval I = Interval._create(self, 
										   x0,  	      #x0 or u
										   x0 + dx, 	  #xf or v
										   self.prob(x0), #prob value calculated in previous iteration
										   0.)			  #cumul value calculated in previous iteration

		
		self._intervals = [I]

		cdef double v = x0 + dx

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
        
		#print( (x, y) for x, y in zip(self._intervals,errors) )
        
		cdef double err
		cdef int k

		for _ in range(300):  
			err = max(errors)
			k   = errors.index(err)

			I = self._intervals[k]
			I1, I2 = I._split()

			self._intervals[k] = I1
			self._intervals.insert(k+1, I2)

			errors = [I.error for I in self]

		self.PAxis = [I.Pu for I in self._intervals] + [self._intervals[-1].Pv]
		self.xAxis = [I.u  for I in self._intervals] + [self._intervals[-1].v ]
		

		#return self


	def __iter__(self):
		yield from self._intervals

	cdef double prob(self, double x):
		return self.p(x)/self.norm

	def __call__(self, double x):
		cdef Interval interval
		for interval in self:
			if x in interval:
				return interval._eval(x)

	cpdef cumul(self, double x):
		for interval in self:
			if x in interval:
				return interval.Pv


	cpdef invCum(self, double r):
		cdef Interval interval
		for interval in self:
			if interval.Pu <= r < interval.Pv:
				return interval._invCum(r)


	def scatter(self):
		

		import matplotlib.pyplot as plt
		plt.figure()
		ax = plt.gca()

		#x_interp = arange(min(self.xAxis), max(self.xAxis), .01)
		y_interp = [self(x) for x in self.xAxis]
		y_true   = [self.p(x)/self.norm for x in self.xAxis]




		#y_cumul = [self.invCum(x) for x in self.xAxis]

		ax.plot(self.xAxis, y_interp)
		ax.scatter(self.xAxis, y_true)
		ax.show()
		#ax.plot(self.xAxis, y_cumul)
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		plt.show()








cdef class Interval(RationalInterpolation):
	"""
	Interval [u, v[
	"""

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
		try:
			self.b = 1 - (self.dP/self.dx)**2 / self.pu / self.pv
			self.a = (self.dP/self.dx) /self.pu - self.b - 1

			self.h = self.dx/50
			self.error = self._calculateError()
		except ZeroDivisionError:
			print("dP - ", self.dP)
			print("dx - ", self.dx)
			print("pu - ", self.pu)
			print("pv - ", self.pv)
			print("u - ",  self.u)
			print("v - ",  self.v)
			raise ZeroDivisionError("IN RITA")
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

	def __str__(self):
		return f"[{self.u}, {self.v}["

	def __repr__(self):
		return f"[{self.u}, {self.v}["


	cdef double _cumul(self):

		return (self.Pu + self.Pv)*.5

		
