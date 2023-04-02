from ..tools.vectors cimport Vector
from ..geometry.main cimport Volume
from ..materials.materials cimport Material, Atom, Shell
cimport numpy as cnp
#cdef double rand()
#cdef object choose(list cumul, list items)
from ..external.mixmax_interface cimport mixmax_engine
from ..types cimport double3

DEF FULL_RECORD = True

from ..types cimport STATE

from libcpp.vector cimport vector
cdef class Particle:
    cdef STATE state
    cpdef get_record_pos(self)
    IF FULL_RECORD:
        cdef vector[STATE] state_record
    ELSE:
        cdef vector[double3] pos_record
        cdef vector[double] E_record

    cdef void invert_axis(self)


    cdef void invert_dire(self)

    cdef void rotateTHETAvers2(self, double cos)
    cdef double ENERGY(self)


    cdef void deposit(self)

    cdef int nSECONDARY
    cdef public object FILE
    cdef long double imfp_T
    cdef object secondary
    

    cdef void rotateAZIMUTH(self, double cos)
    
    
    cpdef add_to_cell(self, object old_cell, object old_points,  int Npoints)

    

    
    
    cdef void throwAZIMUTH(self)
    
    cdef double rsqrt(self, double x)
    
    cdef void normalize(self)
    
    cdef void rotateTHETA(self, long double cos)
    
    # cdef int propagate(self) except -1
    #cdef int change_direction(self, double cos, double phi)
    cdef void update_references(self)
    
    cdef void move(self, double L)

    #cdef double rsqrt(self, double x)

    cdef void record(self)
    cdef void _run(self, mixmax_engine* genPTR)