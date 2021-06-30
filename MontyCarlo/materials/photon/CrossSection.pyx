# cython: profile=True


# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 14:41:41 2020

@author: Rui Campos
"""
print(">>>> IMPORTING material.photon.CrossSection.pyx")








from ...tools.data import getAxis
from ...tools.interpol1 cimport LinLinInterpolation

from numba import njit


def getCS(tuple ID, int Z):
    from ..database import EPDL, EADL
    
    cdef object TABLE = EPDL[Z-1]

    E, CS = TABLE[ID].X, TABLE[ID].Y #MeV and Barn
    
    return E*1e6, CS*1e-24  #changing units to cm^2

#from ..electron.main import eax as _eax

from ..._init import eax


import numpy as np
# cdef class IMFP:
    
#     def __init__(self, LinLinInterpolation[:] interpolators, list coefs, double Nk):
#         self.interpolators = np.array(interpolators)
#         self.coefs = np.array(coefs)
#         self.N = len(coefs)
#         self.Nk = Nk
    
#     cdef double _eval(IMFP self, double E):
#         cdef double sum = 0
#         cdef int i
#         cdef LinLinInterpolation interp
#         for i in range(self.N):
#             interp = self.interpolators[i]
#             sum += interp._eval(E)*self.coefs[i]
        
#         return sum*self.Nk
    
#     def eval(IMFP self, double E):
#         return self._eval(E)
    
#     def __call__(self, double E):
#         return self._eval(E)
        



# cdef IMFP getMFP(tuple ID, dict formula, double density):
#     from ..database import EADL, Na
    
#     Am = 0
        
#     for Z in formula:
#         #print(EADL[Z-1]['Aw'])
#         Am += formula[Z]*EADL[Z-1]['Aw']
        
#     #print(Am)
        
#     Nk  = density * Na / Am #/ 1.660540/1e-24
    
#     #K = 1/N

#     cdef LinLinInterpolation[:] interpolators = np.array([getCS(ID, Z) for Z in formula])
#     coefs = list(formula.values())



#     return IMFP(interpolators, coefs, Nk)

def getLinLin(x, y):
    m = np.diff(y)/np.diff(x)
    
    #y  = m*x - m*x[i] + y[i]
    #m*(x - x[i]) + y[i]
    return - m*x[:-1] + y[:-1], m


from scipy.interpolate import CubicSpline


cimport cython
def REBUILD_CSLOGIC(CLS, imfpA, imfpB):
    cdef CSLOGIC new = CLS({}, 0, pickle = True)
    new.imfpA = imfpA
    new.imfpB = imfpB
    return new

@cython.auto_pickle(True)
cdef class CSLOGIC:
    def __init__(self, tuple ID, formula, double density):
        


        imfpA = []
        imfpB = []
        
        for Z, Ni in formula.items():
            E, CS = getCS(ID, Z)
            imfp = formula.N * CS * Ni                   #1/cm
            splines = CubicSpline(E, imfp, extrapolate = False)
            
            imfp = splines(eax)
            test = np.isinf(imfp)
            
            if np.any(test):
                print("err")
                raise ValueError("CSLOGIC: infinite values in one of the cross sections")
            
            np.nan_to_num(imfp, nan=0, copy = False )
            
            a, b = getLinLin(eax, imfp)
            
            imfpA.append(a)
            imfpB.append(b)
            
        self.imfpA = sum(imfpA)
        self.imfpB = sum(imfpB)
        
        
    def __reduce__(self):
        d = dict()
        d['imfpA'] = np.asarray(self.imfpA)
        d['imfpB'] = np.asarray(self.imfpB)
        
        return REBUILD_CSLOGIC, (self.__class__, d['imfpA'], d['imfpB']) 

        
        
    def __setstate__(self, d):
        self.imfpA = d['imfpA']
        self.imfpB = d['imfpB']
        
        
        
        
        
        
        #self.mfp = getMFP(ID, formula, density)


























































# cdef class CSLOGIC:
#     cdef public IMFP mfp
#     """
#     PEMELOPE SECTION: https://drive.google.com/file/d/15gqkmhli03I00n9vdcN5ygDSbBhaWPbg/
#     """
#     def __init__(self, ID, formula, density):
        
#         #from . import CrossSection
#         self.mfp = getMFP(ID, formula, density)

    # @njit
    # def mu(E):
    #     sum = 0
    #     for i in range(N):
    #         sum += interpolators[i](E)*coefs[i]
    #     return sum*Nk
            
        
        
        
    #     #return Nk*sum(interpolators[i](E)*coefs[i] for i in range(N)) 
    
    # return mu


# def getMFP(ID, formula, density):
    
#     from ..database import EADL, Na
    
#     Am = 0
        
#     for Z in formula:
#         Am += formula[Z]*EADL[Z-1]['Aw']
        
        
#     N  = density * Na / Am
#     CS = composeCS(ID, formula)
    
#     def mfp(E):
#         return 1/N/CS(E)
    
#     return mfp