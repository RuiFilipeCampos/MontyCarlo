




from . import search
from numpy import *


class LinLinInterpolation:
    def __init__(self, xAxis, yAxis):
        self.xAxis, self.yAxis = map(array, (xAxis, yAxis))
        
        self.m = diff(self.yAxis)/diff(self.xAxis)
        
    def __call__(self, x):
        #print(self.xAxis)
        i = search.sortedArrayDOUBLE(self.xAxis, x) #
        #i = searchsorted(self.xAxis, x)
        return self.yAxis[i] + self.m[i]*(x - self.xAxis[i])
    
    
    
    
class LogLogInterpolation:
    def __init__(self, xAxis, yAxis):
        self.xAxis, self.yAxis = map(array, (xAxis, yAxis))
        self.xAxis, self.yAxis = map(log10, (xAxis, yAxis))
        
        self.m = diff(self.yAxis)/diff(self.xAxis)
        
    def __call__(self, x):
        #print(self.xAxis)
        i = searchsorted(self.xAxis, x)
        return 10**(self.yAxis[i] + self.m[i]*(x - self.xAxis[i]))


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

