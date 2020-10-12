from ..tools.vectors cimport Vector
from ..geometry.primitives cimport Volume



cdef object choose(list cumul, list items)


cdef class Particle:
	cdef public double E, theta, phi
	cdef public bint simulate_secondary
	cdef public Vector pos, ex, ey, ez
	cdef public Volume space
	cdef public Volume current_region
	cdef int propagate(self) except -1
	cdef change_direction(self, double cos, double phi)
	cdef public update_coefs(self)
