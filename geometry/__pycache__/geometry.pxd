cdef class Volume:
	cdef public int nothing

cdef class VolumeNode:
	cdef public Volume vol1, vol2
	cdef public object condition
	cdef public int imp

cdef class EnclosingVolume:
	cdef public Volume shape
	cdef public list interior
	cdef public object condition

cdef class UnboundedVolume:
	cdef public int imp
	cdef public list volumes
	cdef public object condition

cdef class Sphere:
	cdef public float radius, x0, y0, z0
	cdef public Vector center
	cdef public int imp