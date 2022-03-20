


cdef class hLinLinInterpolation:
    cdef double [:] _xAxis, _yAxis, _m
    cdef double _xMIN, _xMAX
    cdef int _N
    cdef int[:, :] lims
    cdef double[:] a, b
    cdef int i
    cdef int _sortedArrayDOUBLE(self, double l, int start, int end)
    cdef double _eval(self, double x)



cdef class LinLinInterpolation:
    cdef double [:] _xAxis, _yAxis, _m
    cdef double _xMIN, _xMAX
    cdef int _N
    
    cdef double _eval(self, double x)
    cdef int getINDEX(self, double x)
    
    cdef double evalbyINDEX(self, int i, double x)   
    

    
cdef class InvRationalInterpolation:
    cdef object p, c
    cdef double[:] X, C
    cdef list params
    cdef int Nparams, Npoints
    

    
    cdef double _eval(self, double r)
  
cdef class FastCubicSpline:
    cdef double xMIN, xMAX
    cdef int N
    cdef double[:] x
    cdef double[:, :] c 
    cdef double _eval(self, double x)



#from libc.math cimport log10

cdef class CSa:
    cdef:
        double xMIN, xMAX
        double[::1] X
        int[:, ::1] LIMS
        double[:, ::1] c
    cdef inline int find_index(self, double x)
    cdef double _eval(self, double x)
from MontyCarlo.tools cimport search


cdef class LogLinInterpolation:
    cdef double [:] _xAxis, _yAxis, _m
    cdef double _xMIN, _xMAX
    cdef int _N
  
    
    
    cdef double _eval(LogLinInterpolation self, double x)
    
    
    
  
    
# class LinearInterpolation:
# 	def __init__(self, xAxis, yAxis):

# 		self.xAxis = xAxis
# 		self.yAxis = yAxis

# 		N = len(xAxis)

# 		self.intervals = [] 

# 		for i in range(N-1):
# 			x0, xf = xAxis[i], xAxis[i+1]
# 			y0, yf = yAxis[i], yAxis[i+1]
# 			self.intervals += [Interval(x0, xf, y0, yf)]

# 	def __call__(self, x):
# 		if x == self.xAxis[0]: return self.xAxis[0]
# 		k = searchsorted(self.xAxis, x)
# 		return self.intervals[k-1](x)

# cdef class Interval:
# 	def __init__(self, x0, xf, y0, yf):
# 		self.m = (yf-y0)/(xf-x0)
# 		self.x0, self.xf = x0, xf
# 		self.y0, self.yf = y0, yf

# 	def __call__(self, x):
# 		return self.y0 + self.m * (x - self.x0)

# 	def __contains__(self, x):
# 		return self.x0 <= x < self.xf

# 	def __repr__(self):
# 		return str((self.x0, self.yf))
# 	def __str__(self):
# 		return self.__repr__()

