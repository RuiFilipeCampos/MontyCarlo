from vectors cimport Vector

cdef double do(Vector v2):
	cdef Vector v1
	v1 = Vector(1., 1., 1.)
	cdef double a = v1.norm()
	return v2.norm()

print(do(Vector(2, 2, 2)))