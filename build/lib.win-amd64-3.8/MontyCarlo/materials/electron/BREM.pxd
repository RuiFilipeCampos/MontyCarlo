


from MontyCarlo.tools cimport search
from MontyCarlo.tools.interpol1 cimport LinLinInterpolation


from ..._random.interface cimport mixmax_engine


ctypedef LinLinInterpolation LLI


cdef class sampler:
    cdef LLI[:] X
    cdef double[:] Xmax
    cdef double[:] E, logE, k
    cdef int En
    cdef double[:] kcr
    cdef double Zeff
    cdef int i
    
    

        
    cdef double _sample(self, double E,  mixmax_engine *genPTR)
    cdef double sample_ds(self,  mixmax_engine *genPTR)
    cdef (double, double) full_sample(self, double E,  mixmax_engine *genPTR)
        
    
 
        
        
        