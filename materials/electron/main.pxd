from ..materials cimport Atom as _Atom

from ..._random.interface cimport mixmax_engine

from numpy cimport ndarray


from scipy.interpolate import CubicSpline
from .. import database as db
from ...tools.interpol1 cimport InvRationalInterpolation, LinLinInterpolation, FastCubicSpline, hLinLinInterpolation
#from libcpp.vector cimport vector

from ...tools cimport search, CubicInverseTransform

#getDCS = lambda Z: db.EEDL(Z)[(9, 8, 0, 0, 9, 22)]

#from numpy import array
#cdef double[::1] _eax
ctypedef CubicInverseTransform.aFastCubicSpline CITA

#cdef int[:, ::1] LIMS
cdef class Electron:
    cdef ndarray imfp
    cdef double[::1] imfpA, imfpB
    
    cdef ndarray softSP
    cdef double[::1]  softSPA, softSPB
    
    cdef ndarray softSTRAGG
    cdef double[::1] softSTRAGGA, softSTRAGGB
    
    
    
    
    cdef Inelastic inelastic
    cdef Brem brem
    cdef Elastic elastic
    cdef CITA gauss
    
    cdef double[::1] Itable
    cdef FastCubicSpline integral 
    cdef LinLinInterpolation invI
    
    cdef double find_wmax(self, double smax, double E0)
    cdef int find_index(self, double E)

from . cimport GOS
from . cimport GOSfinal
cdef class Inelastic:
    cdef double[::1] sIMFP1, sIMFP2
    cdef _Atom[::1] arr_atoms
    cdef double[::1] imfpA, imfpB
    cdef double[:, ::1] arr
    cdef int lenposs
    cdef GOSfinal.gosMolecule gosMOLECULE
    cdef GOS.CMolecule GOSmodel
    cdef ndarray fullSP, softSP, imfp, softSTRAGG, fullSTRAGG
    #cdef hLinLinInterpolation softSP, softSTRAGG, imfp
        
    cpdef (double, double, double) get(Inelastic self, double E, bint record)
    cpdef double SP(Inelastic self, double E, bint record)

    
    
from . cimport BREM
cdef class Brem:
    cdef double[::1] imfpA, imfpB
    cdef ndarray imfp, softSP, softSTRAGG, fullSP, fullSTRAGG

    cdef BREM.sampler sampler
    
    cpdef (double, double, double) get(self, double E)




cdef class DIST:
    cdef double rc, mu_c, T1, T2
    cdef double sample(self, mixmax_engine *genPTR)

cdef class Elastic:
    cdef double[:] imfpA, imfpB
    cdef ndarray imfp
    cdef str path
    cdef DIST[::1] DISTRIBUTIONS
    #cdef hLinLinInterpolation imfp0, imfp1, imfp
    cpdef double getSample(self, double E)
    cdef double sMFP(self, double E)
    cdef double[::1] sIMFP1A, sIMFP1B, sIMFP2A, sIMFP2B

    cdef object remove_duplicates(self, ndarray x, ndarray Y)














# cdef class DCS:
#     cdef list __container__
#     cdef vector[double] E
#     cdef LinLinInterpolation T
#     cdef int N
    
#     def __init__(self, Z):
#         self.__container__ = []
        
#         cdef dict xE, pE
#         el = self.getDCS(Z)
#         xE, pE = el.Y1, el.Y2
        
#         cdef list E = list(xE.keys())
#         self.N = len(E)
#         cdef double e

            
        

#         self.__container__ = []
#         cdef double[:] x
#         for e in E:
#             self.E.push_back(e)
#             x = xE[e]
            
            
#             p = CubicSpline(x, pE[e])
#             sampler = InvRationalInterpolation(p, min(x), max(x))
#             self.__container__.append(sampler)
        
        
#         TCS = self.getTCS(Z)
#         self.T = LinLinInterpolation(TCS.X, TCS.Y)
        
#     cdef InvRationalInterpolation _eval(self, double E):
#         i = search._sortedListDOUBLE(self.E, E, 0, self.N)
#         return self.__container__[i]
    
#     def eval(self, double E):
#         return self._eval(E)
    

        
    
#     @classmethod
#     def getDCS(cls, Z):
#         return db.EEDL(Z)[(9, 8, 0, 0, 9, 22)]
    
#     @classmethod
#     def getTCS(cls, Z):
#         return db.EEDL(Z)[(9, 8, 0, 0, 0, 0)]
    





# import numpy as np
# from numpy.random import rand



# cdef class Elastic:
#     cdef vector[double] coefs
#     cdef DCS[:] elements
#     cdef int N
#     cdef vector[double] cache
#     cdef double number_dens
#     cdef LinLinInterpolation[:] transport
    
#     def __init__(self, formula):
        
#         a = list(formula.values())
        
#         for x in a:
#             self.coefs.push_back(x)
#             self.cache.push_back(0.)
            

        
        
#         self.elements = array([DCS(Z) for Z in formula])
        
#         transport = []
        
#         for Z in formula:
#             data = db.EEDL(Z)[()]
#             transport.append(LinLinInterpolation(data.X, data.Y))
            
#         self.transport = array(transport)
        
#         self.N = len(self.elements)
#         self.number_dens = 1
        
        
#     cdef double _imfp(self, double E, double S):
#         cdef double result = 0
#         cdef int i
#         cdef double sigma_el, sigma_el1, chosen
        
        
#         for i in range(self.N):
            
#             sigma_el  = self.elements[i].T._eval(E) 
#             sigma_el1 = self.transport[i]._eval(E)
            
#             chosen = min(sigma_el, max(sigma_el1/C1, S/E / C2)) * self.coefs[i]
            
#             self.cache[i] = chosen
#             result += chosen
            
            
#         return result*self.number_density

#     def choose_element(self):
#         return self._choose_element()

    
#     def imfp(self, E):
#         return self._imfp(E)
        
    
#     cpdef _choose_element(self):
#         cdef double r = rand()*self.cache[self.N-1]
#         print(r)
#         cdef int i
#         for i in range(self.N):
#             if self.cache[i] < r:
#                 return self.elements[i]





# water = Elastic({1:2, 8:1})
# print(water.imfp(1))