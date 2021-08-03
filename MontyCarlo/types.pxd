from .external.mixmax_interface cimport mixmax_engine

cdef double nan

cdef struct double3:
    double x, y, z

cdef struct STATE:
	mixmax_engine* genPTR
	void *current_region
	double3 pos
	double3 dire
	double3 axis
	double E
	double L 
	double last_displacement

cdef class py_state:
	cdef STATE