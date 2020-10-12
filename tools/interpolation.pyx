# cython: profile=True

cimport numpy as cnp
from .others import searchsorted

cdef class LinearInterpolation:
	"""
	Composition phase: holds yAxis for multiplication.
	"""
	cdef cnp.ndarray xAxis
	cdef cnp.ndarray yAxis
	cdef int Z
	cdef int N
	cdef double x_min, x_max
	cdef list intervals

	def __init__(self, 
				 int Z, 
				 cnp.ndarray[cnp.float_t, ndim=1] xAxis, 
				 cnp.ndarray[cnp.float_t, ndim=1] yAxis):

		self.xAxis = xAxis
		self.yAxis = yAxis
		self.Z = Z

	def __hash__(self):
		return self.Z

	def multiply(self, other):
		self.yAxis = self.yAxis*other


	def final_init(self):
		"""
		Creates intervals.
		"""

		self.N = len(self.xAxis) - 1
		self.x_max = self.xAxis[-1]
		self.x_min = self.xAxis[0]

		self.intervals = []

		for i in range(self.N):
			x0, xf = self.xAxis[i], self.xAxis[i+1]
			y0, yf = self.yAxis[i], self.yAxis[i+1]
			try:
				self.intervals += [Interval(x0, xf, y0, yf)]
			except ZeroDivisionError:
				print(f"Zero difivision error ignored: Interval({x0}, {xf}, {y0}, {yf})")
				pass


	def eval(self, x):
		return self._eval(x)

	cdef _eval(self, double x):
		if x <  self.x_min: return 0
		if x == self.x_min: return self.x_min
		if x == self.x_max: return self.x_max
		cdef int k 
		k = searchsorted(self.xAxis, x, 0, self.N)
		return self.intervals[k-1].eval(x)



cdef class Interval:
	cdef double x0, xf, yf, y0, m

	def __init__(self, double x0, double xf, double y0, double yf):
		self.m = (yf-y0)/(xf-x0)
		self.x0, self.xf = x0, xf
		self.y0, self.yf = y0, yf

	def eval(self, x):
		return self._eval(x)

	cdef public double _eval(self, double x):
		return self.y0 + self.m * (x - self.x0)

	def __repr__(self):
		return str((self.x0, self.yf))
		
	def __str__(self):
		return self.__repr__()