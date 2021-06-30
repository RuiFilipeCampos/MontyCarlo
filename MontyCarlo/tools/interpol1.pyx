# cython: profile=False






class MAP(dict):
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError


from . cimport search

import numpy as np



cimport cython
#@cython.boundscheck(False)  # Deactivate bounds checking
#@cython.wraparound(False)   # Deactivate negative indexing.
def newLinLinInterpolation(xAxis, yAxis):
    return LinLinInterpolation(xAxis, yAxis)






def rebuildLinLinInterpolation(this):
    cdef LinLinInterpolation self
    self = <LinLinInterpolation> LinLinInterpolation.__new__(LinLinInterpolation)
    self._xAxis = this._xAxis
    self._yAxis = this._yAxis

    self._N = this._N

    self._m = this._m

    self._xMIN = this._xMIN
    self._xMAX = this._xMAX
    return self


## FastCubicSpline



@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef class LinLinInterpolation:
    #cdef double [:] xAxis, yAxis, m
    #cdef int N
    def __reduce__(self):
        this = MAP()
        from numpy import array

        this._xAxis = array(self._xAxis)
        this._yAxis = array(self._yAxis)

        this._N = self._N

        this._m = array(self._m)

        this._xMIN = array(self._xMIN)
        this._xMAX = array(self._xMAX)
        return rebuildLinLinInterpolation, (this, )
    
    def __init__(self, xAxis, yAxis):
        xAxis, yAxis = map(np.array, (xAxis, yAxis))
        
        DIFFS = np.diff(xAxis)
        
        INCLUDE = np.array([True] + list(DIFFS != 0))
        self._xAxis = xAxis[INCLUDE]
        self._yAxis = yAxis[INCLUDE]
        self._N = len(self._xAxis) - 1   
        self._m = np.diff(self._yAxis)/np.diff(self._xAxis)
        self._xMIN = np.min(self._xAxis)
        self._xMAX = np.max(self._xAxis)
        
        
        
    
    
    cdef double _eval(LinLinInterpolation self, double x):
        #print(self.xAxis)n
        
        if x < self._xMIN:
            return 0
        
        if x > self._xMAX:
            return 0
        
        
        cdef int i = search._sortedArrayDOUBLE(self._xAxis, x, 0, self._N)
        
        

        return self._yAxis[i] + self._m[i]*(x - self._xAxis[i])
    
    
    
    
    
    cdef int getINDEX(self, double x):
        if x < self._xMIN:
            return -2
        
        if x > self._xMAX:
            return -2
        
        return search._sortedArrayDOUBLE(self._xAxis, x, 0, self._N)
    
    cdef double evalbyINDEX(self, int i, double x):
        if i == -2:
            return 0
        return self._yAxis[i] + self._m[i]*(x - self._xAxis[i])
    
    
    def eval(self, double x):
        return self._eval(x)
    
    
    #PYTHON INTERFACE
    def __call__(self, double x):
        return self._eval(x)
    
    def test(self, double x, int N):
        cdef int i
        for i in range(N):
            self._eval(x)
    
    @property
    def xAxis(self):
        return self._xAxis
    
    @property
    def yAxis(self):
        return self._yAxis
    
    @property
    def N(self):
        return self._N
    
    @property
    def yMAX(self):
        return self._yMAX
    
    @property
    def xMAX(self):
        return self._xMAX
    
    
    
from libc.math cimport log10



import numpy as np


cimport cython

@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef class LogLinInterpolation:

    
    """
    xAxis is logspaced
    """
    def __init__(self, xAxis, yAxis):
        xAxis, yAxis = map(np.array, (xAxis, yAxis))
        
        DIFFS = np.diff(xAxis)
        
        INCLUDE = np.array([True] + list(DIFFS != 0)) 
        self._xAxis = xAxis[INCLUDE]
        self._yAxis = yAxis[INCLUDE]
        self._N = len(self._xAxis) - 1   
        self._m = np.diff(self._yAxis)/np.diff(self._xAxis)
        self._xMIN = np.min(self._xAxis)
        self._xMAX = np.max(self._xAxis)
        
        
        
    
    
    cdef double _eval(LogLinInterpolation self, double x):
        #print(self.xAxis)n
        x = log10(x)
        
        if x < self._xMIN:
            return 0
        
        if x > self._xMAX:
            return 0
        
        
        cdef int i = search._sortedArrayDOUBLE(self._xAxis, x, 0, self._N)
        
        
        #i = searchsorted(self.xAxis, x)
        
        if i > self._N - 1:
            return 0
        
        return self._yAxis[i] + self._m[i]*(x - self._xAxis[i])
    
    
    def eval(LogLinInterpolation self, double x):
        return self._eval(x)
    
    
    #PYTHON INTERFACE
    def __call__(LogLinInterpolation self, double x):
        return self._eval(x)
    
    
    @property
    def xAxis(self):
        return self._xAxis
    
    @property
    def yAxis(self):
        return self._yAxis
    
    @property
    def N(self):
        return self._N
    
    @property
    def yMAX(self):
        return self._yMAX
    
    @property
    def xMAX(self):
        return self._xMAX
    
    
    
def rebuildFastCubicSpline(this):
    cdef FastCubicSpline self
    self = <FastCubicSpline> FastCubicSpline.__new__(FastCubicSpline)
    self.xMIN = this.xMIN
    self.xMAX = this.xMAX
    self.x = this.x
    self.c = this.c

    self.N = this.N
    return self

    
    
    
from scipy.interpolate import CubicSpline
@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef class FastCubicSpline:

    def __reduce__(self):
        this = MAP()
        this.xMIN = self.xMIN
        this.xMAX = self.xMAX
        from numpy import array
        this.x = array(self.x)
        this.c = array(self.c)

        this.N = self.N
        return rebuildFastCubicSpline, (this,)
    
    def __init__(self, double[:] x, double[:] y):
        cdef object spline = CubicSpline(x, y)
        self.xMIN = min(x)
        self.xMAX = max(x)
        self.x = spline.x
        self.c = spline.c
        self.N = len(x) - 1
        
    cdef double _eval(self, double x):
        
        if x < self.xMIN:
            return 0
        
        if x > self.xMAX:
            return 0
        
        
        cdef int i = search._sortedArrayDOUBLE(self.x, x, 0, self.N) - 1
        
        if i == -1: 
            i = 0

        cdef double y = 0
        cdef double xi = self.x[i]
        cdef int k
        for k in range(0, 4):
            y += self.c[k, i]*(x - xi)**(3-k)
        return y
    
    def eval(self, double x):
        return self._eval(x)
    





cdef int get_exp(double x):
    cdef int exp;
    frexp(x, &exp);
    return exp;

    
cdef extern from "<math.h>" nogil:
    double frexp(double x, int* exponent)
    
    

def rebuildCSa(this):
    cdef CSa self
    self = <CSa> CSa.__new__(CSa)
    self.xMIN = this.xMIN
    self.xMAX = this.xMAX
    
    self.X = this.X
    self.c = this.c
    self.LIMS = this.LIMS
    return self

cimport numpy as cnp

cdef class CSa:

    def __reduce__(self):
        this = MAP()
        this.xMIN = self.xMIN
        this.xMAX = self.xMAX
        
        this.X = np.array(self.X)
        this.c = np.array(self.c)
        this.LIMS = np.array(self.LIMS)
        return rebuildCSa, (this,)

    def __init__(self, cnp.ndarray X, cnp.ndarray Y):
        cdef object spline = CubicSpline(X, Y)
        self.xMIN = min(X)
        self.xMAX = max(X)
        
        self.X = spline.x
        self.c = spline.c
       # self.N = len(x) - 1
        

        hashed =  np.array([get_exp(x) for x in self.X], dtype = int)
        indexes = np.arange(0, len(hashed))
        Imax = int(max(hashed))
        
        
        
        lims = [np.array([0, 0, 0], dtype = int)]
        
        
        cdef int i 
        for i in range(Imax + 1): #every possible value of the hash, index = hash
            selected = indexes[hashed == i]
            n = len(selected)
            if n == 0: # either out of bounds or no values in this range
                #if no values in this range, interpolate using last interval
                n_last = lims[-1][2]
                if n_last == 0: #out of bounds
                    lims.append(np.array([0, 0, 0], dtype = int))
                    continue
                i_last = lims[-1][1]
                lims.append(np.array([i_last, i_last, 1], dtype = int))
                continue
        
            lims.append(np.array([selected[0], selected[-1] , n], dtype = int))
        
        self.LIMS = np.array(lims[1:], dtype = int) ### memory view defined in pxd, cdef double[::1] EAX
        
        
        
        
        
        
        
        
        
        
        
        
        

    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)  # Deactivate bounds checking
    @cython.wraparound(False)   # Deactivate negative indexing.
    cdef double _eval(self, double x):
        
        if x < self.xMIN:
            return 0
        
        if x > self.xMAX:
            return 0
        
        cdef int i = self.find_index(x)

        cdef double y = 0
        cdef double xi = self.X[i]
        cdef int k
        for k in range(0, 4):
            y += self.c[k, i]*(x - xi)**(3-k)
        return y
    
    
    
    
    cdef inline int find_index(self, double x):
        cdef int i;
        frexp(x, &i);

        
        
        #cdef int i = get_exp(E)
        
        # if LIMS[i, 2] is 0:
        #     raise RuntimeError("OUT OF BOUNDS")
        cdef int k = self.LIMS[i, 2]
        if k is 1:
            return self.LIMS[i, 0]
        
        if k == 2:
            i = self.LIMS[i, 0]
            if x <= self.X[i + 1]: return i
            return i + 1
        
        if k is 3:
            i = self.LIMS[i, 0]
            if x <= self.X[i + 1]: return i
            if x <= self.X[i + 2]: return i + 1
            return i + 2
        
        if k is 4:
            i = self.LIMS[i, 0]
            if x <= self.X[i + 1]: return i
            if x <= self.X[i + 2]: return i + 1
            if x <= self.X[i + 2]: return i + 2
            return i + 3
        
        cdef int START, END, MID
        START = self.LIMS[i, 0]
        END   = self.LIMS[i, 1] 
        
        cdef double xMID
        
        #do binary search 
        while START <= END:
            #find middle
            MID = START + (END - START)//2  
            
            xMID = self.X[MID]
            
            if x is xMID: #found the value
                return MID
            
            if x < xMID: # discard right side
                END = MID - 1 # do not include mid
                continue
            
            START = MID + 1
        return END 

        






    def eval(self, double x):
        return self._eval(x)
    
    
    
    
    
    
    
    
    
    
    
    
    
class LogLogInterpolation:
    def __init__(self, xAxis, yAxis):
        self.xAxis, self.yAxis = map(np.array, (xAxis, yAxis))
        self.xAxis, self.yAxis = map(np.log10, (xAxis, yAxis))
        
        self.m = np.diff(self.yAxis)/np.diff(self.xAxis)
        
    def __call__(self, x):
        #print(self.xAxis)
        i = np.searchsorted(self.xAxis, x)
        return 10**(self.yAxis[i] + self.m[i]*(x - self.xAxis[i]))












from scipy.integrate import quad
from numpy.random import rand
from numpy import searchsorted
#from numpy import logspace, log10, array
import cython
from libc.math cimport log10



@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.boundscheck(True)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef class InvRationalInterpolation:

    
    def __init__(self, f:callable, x0:float, xf:float):
        """
        f -> cubic spline interpolation
        """
        
        A = f.integrate(x0, xf)
        self.p = lambda x: f(x)/A
        self.c = lambda x: f.integrate(x0, x)/A
        
        
        #self.p = normalize(f, x0, xf)
        #self.c = make_cumul(f, x0, xf)
        
        self.X = np.logspace(log10(x0), log10(xf), 500)
        self.C = np.array([self.c(x) for x in self.X])
        self.Npoints = len(self.C)
        
        #self.X = [x0 + rand()*(xf-x0) for _ in range(10)]
        
       
        
        self.params = []
        for i in range(len(self.C)-1):
            self.params.append( self.getab(i) )
            
        self.Nparams = len(self.params)
    
    
    def b(self, int i):
        C, X, p = self.C, self.X, self.p
        
        dC = C[i+1] - C[i]
        dx = X[i+1] - X[i]
        
        return 1 - (dC/dx)**2 / p(X[i+1]) / p(X[i])
    
    def getab(self, int i):
        C, X, p = self.C, self.X, self.p
        
        cdef double dC = C[i+1] - C[i]
        cdef double dx = X[i+1] - X[i]
        
        cdef double b = self.b(i)
        
        return dC/dx / p(X[i]) - b - 1, b
    
    cdef double _eval(self, double r):
        cdef int i = search._sortedArrayDOUBLE(self.C, r, 0, self.Npoints)
        
        if i > self.Nparams-1:
            return 1
        
        cdef double a, b
        a, b = self.params[i]
        
        cdef double nu = (r - self.C[i])/(self.C[i+1] - self.C[i])
        
        cdef double A = (1 + a + b)*nu
        cdef double B = 1 + a*nu + b*nu**2
        
        cdef double x0 = self.X[i]
        cdef double xf = self.X[i+1]
        
        return x0 + A/B * (xf-x0)        

    def eval(self, double r):
        return self._eval(r)
    
    def __call__(self, r):
        i = search._sortedArrayDOUBLE(self.C, r, 0, len(self.C))
        #i = searchsorted(self.C, r, side="right") + 1
        if i > len(self.params)-1:
            return 1
        
        
        a, b = self.params[i]
        
        nu = (r - self.C[i])/(self.C[i+1] - self.C[i])
        
        A = (1 + a + b)*nu
        B = 1 + a*nu + b*nu**2
        
        x0 = self.X[i]
        xf = self.X[i+1]
        
        return x0 + A/B * (xf-x0)


























# class LinearInterpolation:
#     def __init__(self, xAxis, yAxis):

#         self.xAxis = xAxis
#         self.yAxis = yAxis  

#         N = len(xAxis)

#         self.intervals = [] 

#         for i in range(N-1):
#             x0, xf = xAxis[i], xAxis[i+1]
#             y0, yf = yAxis[i], yAxis[i+1]
#             self.intervals += [Interval(x0, xf, y0, yf)]

#     def __call__(self, x):
#         if x == self.xAxis[0]: return self.xAxis[0]
#         k = searchsorted(self.xAxis, x)
#         return self.intervals[k-1](x)

# cdef class Interval:
#     def __init__(self, x0, xf, y0, yf):
#         self.m = (yf-y0)/(xf-x0)
#         self.x0, self.xf = x0, xf
#         self.y0, self.yf = y0, yf

#     def __call__(self, x):
#         return self.y0 + self.m * (x - self.x0)

#     def __contains__(self, x):
#         return self.x0 <= x < self.xf

#     def __repr__(self):
#         return str((self.x0, self.yf))
#     def __str__(self):
#         return self.__repr__()

from libc.math cimport log

@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef class hLinLinInterpolation:
    #cdef double [:] xAxis, yAxis, m
    #cdef int N
    
    
    def __init__(self, xAxis, yAxis):
        xAxis, yAxis = map(np.array, (xAxis, yAxis))
        
        DIFFS = np.diff(xAxis)
        
        INCLUDE = np.array([True] + list(DIFFS != 0))
        _xAxis = xAxis[INCLUDE]
        _yAxis = yAxis[INCLUDE]
        self._N = len(_xAxis) - 1   
        _m = np.diff(_yAxis)/np.diff(_xAxis)
        
        
        self._xMIN = np.min(_xAxis)
        self._xMAX = np.max(_xAxis)
        
        hashed = np.floor(np.log(_xAxis))
        indexes = np.arange(0, len(hashed))
        
        Imax = int(max(hashed))
        
        lims = []
        for i in range(Imax + 1):
            _indexes = indexes[hashed == i]
            
            n = len(_indexes)
            
            if n < 1:
                lims.append(np.array([0, 0, n], dtype = int))
                continue
            
            _lims = [_indexes[0], _indexes[n-1], n]
            lims.append(np.array(_lims, dtype = int))
        
        
        #print(lims)
        self.lims = np.array(lims)
        
        self.a = _yAxis[:-1] - _xAxis[:-1]*_m
        self.b = _m
        self._xAxis = _xAxis
        self._yAxis = _yAxis
        self._m = _m
        
        
        
        #self._yAxis[i] + self._m[i]*(x - self._xAxis[i])
        
    
    
    cdef double _eval( self, double x):
        #print(self.xAxis)n
        
        if x < self._xMIN:
            return 0
        
        if x > self._xMAX:
            return 0
        
    
        
        self.i = <int> log(x)

        
        if self.lims[self.i, 2] < 1:
            return self._yAxis[self.i] + self._m[self.i]*(x - self._xAxis[self.i])
        
            #return self.a[self.i] + self.b[self.i]*x
        
        #self.i = self._sortedArrayDOUBLE(x, 
        #                              self.lims[self.i, 0],
        #                              self.lims[self.i, 1])
        
        
        

        cdef int mid
        cdef int start = self.lims[self.i, 0]
        cdef int end = self.lims[self.i, 1]
        
        while start <= end:
            #print(f"Subarray: {list(arr[start:end+1])}")
            mid = start + (end - start)//2
            
            if x == self._xAxis[mid]:
                self.i = mid
                break
            
            if x < self._xAxis[mid]:
                end = mid - 1
            else:
                start = mid +1
        else:
            self.i = end
            
        if self.i >= self._N:
            self.i = self._N - 1

        #return self.a[self.i] + self.b[self.i]*x
        return self._yAxis[self.i] + self._m[self.i]*(x - self._xAxis[self.i])
    
    
    
    
    
    
    
    
    
    cdef int _sortedArrayDOUBLE(self, double l, int start, int end):
        if start > end:
            return -1
        
        cdef int mid_point = <int> (start + (end - start)/2)
        
        cdef int i
        
        if l < self._xAxis[mid_point]:
            i =  self._sortedArrayDOUBLE( l, start, mid_point-1)
            if i == -1: return mid_point - 1
            else:       return i
        
        if self._xAxis[mid_point] < l:
            i = self._sortedArrayDOUBLE( l, mid_point+1, end)
            if i == -1: return mid_point
            else:       return i
            
        else:
            return mid_point
        
    def test(self, double x, int N):
        cdef int i
        for i in range(N):
            self._eval(x)
    
    
    def eval(self, double x):
        return self._eval(x)
    
    
    #PYTHON INTERFACE
    def __call__(self, double x):
        return self._eval(x)
    
    
    @property
    def xAxis(self):
        return self._xAxis
    
    @property
    def yAxis(self):
        return self._yAxis
    
    @property
    def N(self):
        return self._N
    
    @property
    def yMAX(self):
        return self._yMAX
    
    @property
    def xMAX(self):
        return self._xMAX