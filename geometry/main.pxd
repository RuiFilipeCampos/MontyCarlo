# distutils: language = c++


from ..materials.materials cimport Material
from ..particles.particle cimport Particle

from ..types cimport double3



cdef struct TRANSFORM:
    double3 a
    double3 b

# cdef class CSGvol


# cdef class AND(CSGvol):
#     cdef CSGvol vol1, vol2
    
#     cdef Interval getINTERSECTION(self, origin, direction)
#     cdef double getSDF(double x, double y, double z)

# cdef class OR(CSGvol):
#     cdef CSGvol vol1, vol2
    
#     cdef Interval getINTERSECTION(self, origin, direction)
#     cdef double getSDF(double x, double y, double z)


# cdef class NOT(CSGvol):
#     cdef CSGvol vol1
    
#     cdef Interval getINTERSECTION(self, origin, direction)
#     cdef double getSDF(double x, double y, double z)


from ..particles.particle cimport STATE

cimport numpy as cnp
cdef class Volume:
    #cdef TRANSFORM tr
    cdef public object mesh
    cdef double x, y, z
    cdef bint opaque
    #cdef int closest
    #cdef Particle p 
    #cdef int Ninner
    #cdef void** inner
    #cdef list tmp_inner
    cdef Volume outer
    cdef Material material
    cdef cnp.ndarray voxels
    cdef void depositUNIFORM(self, STATE& state, double SP)
    cdef void depositLOCAL(self, double3& pos, double E)
    cdef void depositRANDOM(self, STATE& state, double E, double tau)
    #cdef public int imp
    #cdef bint NEW_VOL
    #cdef double SDF(self, double x, double y, double z)
    #cdef bint move(self, double L)
    #cdef void* find_inner(self, double x, double y, double z)
    #cdef cnp.ndarray _generate_points(self, int Nrays)
    #cdef bint trace(self, double L)
   # cdef void move_to_surface(self)
    # cdef Point _getIntersection(self, Vector origin, Vector direction)
    # cpdef Interval _getINTERSECTION(self, Vector origin, Vector direction)
    # cdef Volume cross(self, Vector pos)
    cdef bint move(self, STATE& state, double SP)
    cdef double SDF(self, double3& pos)
    cdef bint is_inside(self, double3& pos)
    cdef void exit(self)
    cdef void depositLocaly(self, double3 pos, double E)


# cdef class Interval:
#     cdef long double t1, t2
#     cdef bint EMPTY
#     cdef Volume REGION
    
#     @staticmethod
#     cdef Interval _new(long double t1, long double t2)
#     cdef Interval _NOT(Interval self)
#     cdef Interval _AND(Interval self, Interval other)
#     cdef Interval _OR(Interval self, Interval other)
        

            
        

