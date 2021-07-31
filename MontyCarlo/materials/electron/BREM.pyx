# cython: profile=False
# cython: annotate=True
# distutils: language = c++ 
# distutils: extra_compile_args = -std=c++11

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



from ..._random.interface cimport mixmax_engine



cimport cython
from libc.math cimport log, fmax
import numpy as np
from MontyCarlo.tools cimport search
from MontyCarlo.tools.interpol1 cimport LinLinInterpolation


#cdef double Wcr = 0.001







# @cython.cdivision(True)
# cdef double urand():
#     cdef double r = rand()
#     return r / RAND_MAX

# def purand():
#     print(urand())
    
    
def rebuildsampler(this):
    cdef sampler self
    self = <sampler> sampler.__new__(sampler)
    self.Zeff = this.Zeff

    self.k    = this.k       # molecule.k
    self.X    = this.X       # np.array(XX)    
    self.Xmax = this.Xmax    # np.array(Xmax)
    self.E    = this.E       # np.array(molecule.E)*1e6
    self.kcr  = this.kcr     # molecule.Wcr/(molecule.E*1e6)
    self.logE = this.logE    # np.log(molecule.E*1e6)

    self.En = this.En
    return self


@cython.boundscheck(False)
@cython.wraparound(False) 
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef class sampler:


    def __reduce__(self):
        this = MAP()
        this.Zeff = self.Zeff

        from numpy import array
        this.k    = array(self.k   )    # molecule.k
        this.X    = array(self.X   )    # np.array(XX)    
        this.Xmax = array(self.Xmax)    # np.array(Xmax)
        this.E    = array(self.E   )    # np.array(molecule.E)*1e6
        this.kcr  = array(self.kcr )    # molecule.Wcr/(molecule.E*1e6)
        this.logE = array(self.logE)    # np.log(molecule.E*1e6)

        this.En = self.En
        return rebuildsampler, (this, )


    
    
    def __init__(self, object molecule):
        self.Zeff = molecule.Zeff
        cdef double[:] X
        self.k = molecule.k
        
        XX = []
        Xmax = []
        for X in molecule.ds:
            Xmax.append(max(X))
            XX.append(LinLinInterpolation(self.k, X))
        
        self.X = np.array(XX)    
        self.Xmax = np.array(Xmax)
        
        
        #molecule.E = np.array(molecule.E)*1e6
        
        self.E = np.array(molecule.E)*1e6

        
        assert len(self.X) == len(self.E)
        self.kcr = molecule.Wcr/(molecule.E*1e6)
        self.logE = np.log(molecule.E*1e6)
        
        self.En = len(self.E) - 1
        
        
    cdef (double, double) full_sample(self, double E,  mixmax_engine *genPTR):
        """Sample electrons fractional energy loss (k) and the emitted photons polar angular with respect to the direction of the electrons movement (theta).
        """
        
        cdef double k = self._sample(E,genPTR)
        #cdef double theta = sample_theta(E, self.Zeff, k)
        return (k , sample_theta(E, self.Zeff, k, genPTR))
        
    cdef double _sample(self, double E, mixmax_engine *genPTR):
        """Sample the electrons fractional energy loss (k).
        """

        self.i = search._sortedArrayDOUBLE(self.E, E, 0, self.En)
        
        if self.E[self.i] == E:
            return self.sample_ds(genPTR)
        
        cdef double logE1, logE2, logE
        assert self.E[self.i] <= E < self.E[self.i+1]
        
        logE1, logE2 = self.logE[self.i], self.logE[self.i+1]
        logE = log(E)
  
        
        #cdef double pi_1, pi_2
        #pi_1 = 
        #pi_2 = (logE - logE1)/log_diff
        
        

        
        if genPTR.get_next_float() > ( logE2 - logE ) / (logE2 - logE1): 
            self.i +=1
        
        return self.sample_ds(genPTR)
        
    cdef double sample_ds(self,  mixmax_engine *genPTR):
        """Sample the electrons fractional energy loss (k) from the chosen X-Section.
        """
        cdef double w
        cdef double kcr = self.kcr[self.i] 
        cdef LLI XX = self.X[self.i]
        cdef double Xmax = self.Xmax[self.i]
        while 1:
            w = kcr**genPTR.get_next_float()
            if genPTR.get_next_float()*Xmax < XX._eval(w):
                return w
    

        
    
   # def sample(self, double E):
     #   return self._sample(E)   
    
  #  def py_full_sample(self, double E):
   #     return self.full_sample(E)
    
    
from libc.math cimport pi, log

@cython.cdivision(True)
cdef double f(double x, double E0):
    return (1 + 1/(pi*E0**2))/(x + 1)**2

@cython.cdivision(True)
cdef double m(double x, double E0, double Z, double r):
    return ((1-r)/(2*E0*r))**2 + (Z**(1/3)/(111*(x+1)))**2

@cython.cdivision(True)
cdef double g(double x, double E0, double Z, double r):

    cdef double A = 4 + log(m(x, E0, Z, r))
    cdef double B = (1 + r**2) - 4*x*r /(x+1)**2
    return 2*r - 3*(1 + r**2) - A*B

@cython.cdivision(True)
cdef double getTheta(double E0, mixmax_engine *genPTR):
    cdef double r = genPTR.get_next_float()
    #cdef double A = r/(1 - r + 1/(pi*E0)**2)
    return (r/(1 - r + 1/(pi*E0)**2))**.5/E0

@cython.cdivision(True)
cdef double getNr( double E0, double Z, double r):
    cdef double g0 = g(0, E0, Z, r)
    cdef double g1 = g(1, E0, Z, r)
    cdef double gg = g((pi*E0)**2, E0, Z, r)
    return max(g0, g1, gg)**-1



@cython.cdivision(True)
cdef double sample_theta(double E0, double Z, double k, mixmax_engine * genPTR):
    cdef double theta, Nr, g_test, r
    Nr = max(g(0, E0, Z, k), 
             g(1, E0, Z, k), 
             g((pi*E0)**2, E0, Z, k))**-1
    while 1:
        r = genPTR.get_next_float()
        theta = (r/(1 - r + 1/(pi*E0)**2))**.5/E0
 
        if genPTR.get_next_float() < Nr*g((E0*theta)**2, E0, Z, k):
            return theta
        
        
        
