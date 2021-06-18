


cimport cython

cdef aFastCubicSpline fromSample(x, y)

cdef class aFastCubicSpline:
    cdef public object spline
    cdef public double[:, :] c
    cdef public double[:] x, cut_offs
    cdef public int [:] rej_indexes
    cdef public double y
    cdef public double dx
    cdef public double N
    cdef public double[:] DX

    cdef double R
    cdef int  i
    cdef public double[:] YY
    
 
   
    cdef double _sample(self)
    
    
        
        
        
        
        
        
    