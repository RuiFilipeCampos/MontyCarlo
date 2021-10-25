# distutils: language = c++

from ..materials.materials cimport Material
from ..particles.particle cimport Particle
from ..types cimport double3
from ..types cimport STATE

cdef struct TRANSFORM:
    double3 a
    double3 b


cimport numpy as cnp
cdef class Volume:
    cdef public object mesh
    cdef double x, y, z
    cdef bint opaque
    cdef Volume outer
    cdef Material material
    cdef cnp.ndarray voxels
    cdef void depositUNIFORM(self, STATE& state, double SP)
    cdef void depositLOCAL(self, double3& pos, double E)
    cdef void depositRANDOM(self, STATE& state, double E, double tau)
    cdef bint move(self, STATE& state, double SP)
    cdef double SDF(self, double3 pos)
    cdef bint is_inside(self, double3& pos)
    cdef void exit(self)
    cdef void depositLocaly(self, double3 pos, double E)


        

            
        

