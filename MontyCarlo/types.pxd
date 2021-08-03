from .external.mixmax_interface cimport mixmax_engine

from numpy cimport ndarray

cdef double nan

cdef struct double3:
    double x, y, z

cdef struct STATE:
    mixmax_engine *genPTR
    void *current_region
    double3 pos
    double3 dire
    double3 axis
    double E
    double L 
    double last_displacement

cdef class py_state:
    cdef STATE state
    cdef mixmax_engine gen # need to keep the generator somewhere...

    # This guys need to be acessible from python.
    cdef public ndarray[ndim=1] pos
    cdef public ndarray[ndim=1] dire
    cdef public ndarray[ndim=1] axis
    cdef public double E
    cdef public double L 
    cdef public double last_displacement
    cdef public long int seed

    # Note: I can't just make `state` public since it contains c types
    # like the `mixmax_engine` and `void *` pointer. 
