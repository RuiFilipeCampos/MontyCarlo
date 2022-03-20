# cython: profile=True

#External Imports

#Internal Imports
from ...tools.RITA cimport RationalInterpolation

cimport cython

from .CrossSection cimport CSLOGIC


from ...tools.interpol1 cimport LinLinInterpolation


from numpy cimport ndarray



cdef class Photon:
    cdef double density
    cdef Coherent coherent
    cdef Incoherent incoherent
    cdef Pairproduction pairproduction
    cdef Tripletproduction tripletproduction
   




cdef class Coherent(CSLOGIC):
    cdef RationalInterpolation FF
    cdef:
        double[:, ::1] xSPLINE
        double[:, ::1] ySPLINE
        double[::1] X, Y
        int[:, ::1] xLIMS, yLIMS
        double xMAX, xMIN
        int xADDER, yADDER
        
        
    cdef int find_index_X(self, double x)
    cdef int find_index_Y(self, double x)

    cdef double evalY(self, double x)
    cdef double evalX(self, double r)



cdef class Incoherent(CSLOGIC):
    cdef object S




cdef class Pairproduction(CSLOGIC):
    cdef:
        double factor
        double CONST
        double alpha
        double a
        double Zeq
        double fC
    



    cpdef double F0(self, double k)


  #      return  (-1.774 - 12.10*a  + 11.18*a2) * (k2)**.5    \
   #             + (8.523 +  73.26*a   − 44.41*a2) * (k2)        \
    #            - (13.52  + 121.1*a  − 96.41*a2) * (k2)**(3/2) \
     #           + (8.946  + 62.05*a  − 64.41*a2) * (k2)**2


    cpdef double g1(self, double b)

    cpdef double g2(self, double b)

    cpdef double g0(self, double k)

    cpdef double b(self, double eps, double k)

    cpdef  (double, double) getPhis(self, double eps, double k)




cdef class Tripletproduction(CSLOGIC):
    cdef double[:, ::1] ALIAS
    



















