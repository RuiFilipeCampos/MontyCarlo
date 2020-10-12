cdef class Vector:
	cdef public double x, y, z, real
	cdef public double norm(self)
	cdef public Vector normalize(self)
	cdef public Vector _mult(self, Vector q)
	cdef public Vector _conj(self)
	cdef public Vector rotateCos(self, Vector axis, double cos)
