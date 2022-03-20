cdef extern from "<math.h>" nogil:
    double frexp(double x, int* exponent)

cdef double EAX[2695] 
cdef int[:, ::1] LIMS
