from ..materials cimport Material

from .coherent.coherent     cimport Coherent
#from .incoherent.incoherent cimport Incoherent
#from .photoelectric         cimport Photoelectric
#from .pairproduction        cimport Pairproduction
#from .tripletproduction     cimport Tripletproduction

cdef class Photon:
	cdef public Material upper_dir
	cdef public object coherent
	cdef public object incoherent
	cdef public object photoelectric
	cdef public object pairproduction
	cdef public object tripletproduction
	cpdef final_init(self, double density)

cdef class CompoundPhoton(Photon):
	cdef bint x
