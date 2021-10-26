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

    # OBLIGATORY
    cdef bint opaque # TRUE = particles get instantly absorbed
    cdef Volume outer # oof, this is meant for the BVH right?
    cdef Material material




    cdef public object mesh
    cdef cnp.ndarray voxels
    cdef void depositUNIFORM(self, STATE& state, double SP)
    cdef void depositLOCAL(self, double3& pos, double E)
    cdef void depositRANDOM(self, STATE& state, double E, double tau)
    cdef bint move(self, STATE& state, double SP)
    cdef double SDF(self, double3 pos)
    cdef bint is_inside(self, double3& pos)
    cdef void exit(self)
    cdef void depositLocaly(self, double3 pos, double E)


        

cdef class BVH(Volume):
	# Workspace
	cdef int Nws
	cdef list tmp_ws
	cdef void** ws
	cdef void** original_ws
	cdef bint has_name


	# Boundary Crossing
	cdef int position_in_outer       #position in outers work space
	cdef bint keep                   #which intersected volume will keep its intersections cached for the next iteration

	# user related
	cdef bint lock 					#prevent user from modifying volume after exit code
	cdef str name
	cdef bint render

	# ray marching
	cdef double sdf                  # nearest distance to this volumes surface

	# Ray Tracing
	cdef bint cache                  # is this volume storing cached intersections?
	cdef intIterator cross           # custom c++ iterator for aiding in simulation with cached intersections
	cdef double particle_position;   # must keep track of particles position along the ray 
        


	cdef void exit(self)
	cdef bint move(self, STATE& state, double SP)
	cdef void depositUNIFORM(self, STATE& state, double SP)
	cdef void depositDISCRETE(self, STATE& state)	
	cdef void depositLOCAL(self, double3& pos, double E)
	cdef void depositRANDOM(self, STATE& state, double E, double tau)

	cdef double main_intersect(self, STATE& state)
	cdef double SDF(self, double3 pos)
	cdef bint is_inside(self, double3 pos)