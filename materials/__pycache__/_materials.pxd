from .photon.photon cimport Photon

cdef str msg

cdef class Material:
	cdef public double A, Z
	cdef public Photon photon
	cdef public bint electron
	cdef public str _name
	cdef public bint created
	cpdef create(self, double density)
