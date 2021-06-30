
from ...tools.interpol1 cimport LinLinInterpolation




# cdef LinLinInterpolation getCS(tuple ID, int Z)


# cdef class IMFP:
#     cdef LinLinInterpolation[:] interpolators
#     cdef int[:] coefs
#     cdef int N
#     cdef double Nk
#     cdef double _eval(self, double E)

# cdef IMFP getMFP(tuple ID, dict formula, double density)

cimport cython
@cython.auto_pickle(True)
cdef class CSLOGIC:
    cdef double[:] imfpA, imfpB

