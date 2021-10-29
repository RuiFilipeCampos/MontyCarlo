




ctypedef (*map_t)(double3& pos, double[:] transformation)


cdef class Primitive(CSGvol):
    cdef bint has_transform
    cdef map_t apply_transform
	cdef double[16] direct_transform
    cdef double[16] inverse_transform