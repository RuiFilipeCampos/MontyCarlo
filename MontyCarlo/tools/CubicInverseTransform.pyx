print("Importing .tools.CubicInverseTransform ...")


__doc__ = """

Cubic Inverse Transform :: An alternative to RITA

"""


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
from scipy.integrate   import cumtrapz, trapz
from scipy.interpolate import CubicSpline

import numpy as np
from numpy import array, linspace, floor
from numpy cimport ndarray
#from libc.math cimport floor
from libc.stdlib cimport rand, RAND_MAX
cimport cython


# INTERNAL IMPORTS



@cython.cdivision(True)
cdef double urand():
    cdef double r = rand()
    return r / RAND_MAX





cdef object remove_duplicates(ndarray x, ndarray Y):
    cdef ndarray u, c, dup
    u, c = np.unique(x, return_counts=True)
    dup = u[c > 1]
    
    cdef bint keep = True
    cdef list new_y = []
    cdef int i
    cdef double y
    
    for i, y in enumerate(Y):
        
        if x[i] in dup:
            if keep:
                new_y.append(y)
                keep = False
                continue
            else: continue
        
        if keep is False:
            keep = True
        new_y.append(y)
        
    for y1, y2 in zip(Y, new_y): print(y1, y2)
    import time
    time.sleep(10000)
    return u, np.array(new_y)



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


from scipy.integrate import quad
    
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

    
cimport numpy as cnp





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
        self.rej_indexes = array([int(x) for x in aliases[:, 2]], dtype = np.int)

        self.N =  (len(self.x) -2) 
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
        
        
    