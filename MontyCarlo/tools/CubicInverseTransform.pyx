

print("Importing .tools.CubicInverseTransform ...")


__doc__ = """
Cubic Inverse Transform :: An alternative to RITA

"""


# Should be imported from `.types`.
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



# EXTERNAL IMPORTS
import numpy as np
from scipy.integrate   import cumtrapz
from scipy.integrate   import trapz
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from numpy import array
from numpy import linspace
from numpy import  floor

cimport numpy as cnp
from numpy cimport ndarray
from libc.stdlib cimport rand
from libc.stdlib cimport RAND_MAX
cimport cython
#from libc.math cimport floor


# INTERNAL IMPORTS
from .main cimport remove_duplicates



@cython.cdivision(True)
cdef double urand():
    cdef double r = rand()
    return r / RAND_MAX



def getDist(X):   
    var = []
    prob = []
    for i in range(len(X) - 1):
        var.append(i)
        prob.append(X[i + 1] - X[i])
    var = array(var)
    prob = array(prob)
    return var, prob/sum(prob)


def makeAlias(X, Y):
    X, Y = X.copy(), Y.copy()
    N = len(Y)
    Y = N*Y
    points = []
    while len(Y) > 0:
        
        i_min, i_max = np.argmin(Y), np.argmax(Y)
        
        
        ymax = Y[i_max]
        ymin = Y[i_min]
        
        xmin = X[i_min]
        xmax = X[i_max]
        
        dy = 1 - ymin
        
        Y[i_max] -= dy
        
        point = [xmin, ymin, xmax]
        points.append(point)
        
        X = np.delete(X, i_min)
        Y = np.delete(Y, i_min)

    
    points = sorted(points)
    return array(points)


def fromCallable(f, a, b, num = 100, endpoints = True):
    if endpoints:
        X = linspace(a, b, num = num)
    else:
        X = linspace(a, b, num = num + 2)
        X = X[1:-1]

    
    Y = [f(x) for x in X]
    Y = array(Y)
    
    
    A = trapz(Y, X)
    Y = Y/A
    
    cumul = cumtrapz(Y, X, initial = 0)
    
    indexes, prob = getDist(cumul)
    A = makeAlias(indexes, prob)
    
    return aFastCubicSpline(cumul,X , aliases = A)


cdef aFastCubicSpline fromSample(x, y):
    y = y/trapz(y, x)
    
    cumul = cumtrapz(y, x, initial = 0)
    
    cumul, x = remove_duplicates(cumul, x)
    indexes, prob = getDist(cumul)
    A = makeAlias(indexes, prob)        
    
    return aFastCubicSpline(cumul, x, aliases = A)

    





def rebuildaFastCubicSpline(this):
    cdef aFastCubicSpline self
    self    = <aFastCubicSpline> aFastCubicSpline.__new__(aFastCubicSpline)
    self.c  = this.c
    self.x  = this.x
    self.DX = this.DX

    self.cut_offs    = this.cut_offs
    self.rej_indexes = this.rej_indexes

    self.N = this.N
    return self




@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef class aFastCubicSpline:
    def __reduce__(self):
        this = MAP()
        from numpy import array
        this.c = array(self.c)
        this.x = array(self.x)
        this.DX = array(self.DX)
        
        this.cut_offs = array(self.cut_offs)
        this.rej_indexes = array(self.rej_indexes)

        this.N = self.N
        return rebuildaFastCubicSpline, (this, )



    def __init__(self, x, y, aliases = None):
        if aliases is None:
            raise RuntimeError("No aliases were given.")
            
            
        spline = CubicSpline(x, y)

        self.c = spline.c
        self.x = x
        self.DX = np.diff(x)/RAND_MAX
        
        self.cut_offs = aliases[:, 1]
        self.rej_indexes = array([int(x) for x in aliases[:, 2]], dtype = np.int32)

        self.N =  (len(self.x) - 2) 
        self.N = self.N / RAND_MAX
        
        

    # cdef double _cumul(self, double x):
    #     pass
    
    
    # cdef double _invCumul(self, double r):
        
    
    cdef double _sample(self):
        self.R = self.N * rand()
        self.i = <int> self.R

        if   self.R - self.i < self.cut_offs[self.i]:
            
            self.dx =  rand() * self.DX[self.i]
            
            self.y = self.c[3, self.i]
            
            self.y += self.c[2, self.i]*self.dx

            self.dx *= self.dx
            self.y += self.c[1, self.i]*self.dx

            self.dx *= self.dx
            self.y += self.c[0, self.i]*self.dx
            
            return self.y      
        
        
        self.i = self.rej_indexes[self.i]

        
        self.y = self.c[3, self.i]
        self.dx = rand()*self.DX[self.i]
        
        self.y += self.c[2, self.i]*self.dx

        self.dx *= self.dx
        self.y += self.c[1, self.i]*self.dx

        self.dx *= self.dx
        self.y += self.c[0, self.i]*self.dx
        
        return self.y 
    
    
 
    def sample(self, int N):
        cdef int i
        for i in range(N):
            self._sample()
            
    def getsample(self, int N):
        cdef int i
        cdef double[:] sample = np.zeros(N)
        for i in range(N):
            sample[i] = self._sample()
        return sample
        
        
    
