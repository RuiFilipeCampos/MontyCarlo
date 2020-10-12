from ..tools.vectors cimport Vector
#from ..materials.materials cimport Material


ctypedef bint (*bool_func)(bint, bint)


cdef class Volume:
	cdef public int N
	cdef bool_func condition
	cdef public Volume outer
	cdef public object material
	cdef public bint _eval(self, double x, double y, double z)
	cdef public list inner
	cpdef Volume fill(self, object material)
	cdef public double getOuterCrossing(self, Vector pos, Vector ez, double L)
	cdef public tuple getCrossings(self, Vector pos, Vector ez, double L)



cdef class VolumeNode(Volume):
	cdef public Volume vol1, vol2
	cdef public int imp
	cpdef list intercept(self, Vector pos, Vector ez, double L)

cdef class EnclosingVolume(Volume):
	cdef public Volume shape
	cdef public list interior

cdef class UnboundedVolume(Volume):
	cdef public Volume ROI
	cdef public int imp
	cpdef list intercept(self, Vector pos, Vector ez, double L)

cdef class Sphere(Volume):
	cdef public double radius, x0, y0, z0
	cdef public Vector center
	cdef public int imp
	#cpdef bint condition(self, double x, double y, double z)
	cpdef list intercept(self, Vector pos, Vector ez, double L)