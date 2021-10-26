# distutils: language = c++
# cython: profile=False

__doc__ = """Contains an abstract class that serves as template for any geometry system.
"""

__author__ = "Rui Campos"

print("Importing `.geometry.main`")



# External Imports
import numpy as np 

from libc.math     cimport fmin, fmax, sqrt , cos, pi
from libcpp.vector cimport vector
from numpy.math    cimport INFINITY
from numpy         cimport ndarray
from libcpp.list   cimport list as cpplist
from libc.stdlib   cimport malloc, free

from cython.operator import dereference, preincrement
from numpy.random import randint

# Internal Imports
#from ..tools.vectors cimport Vector
#from ..random.mixmax.interface cimport mixmax_engine
from ..types cimport STATE
from ..types cimport double3

cdef double MINdr = 0 # deprecated 



# LEGACY:
cdef struct Ray:
    double3 O
    double3 D

cdef struct intersection:
    double3 pos
    double d
####



cdef class Volume:
    """Abstract class that dictates the minimum functionality that any geometry implementation needs
    to be compatible with the physics engine.
    """

    def __init__(self, material = None):
        if material is None:
            raise RuntimeError("This volume has not been filled with a material.")

        self.material = material


        
    
    
    def fill(self, mat):
        """Sets attribute `material`.
        
        Note:
            This actually needs some work. It is important to check if the user has already
            filled the volume and raise an exception if so. Also, if it hasn't filled the
            volume. 
        """
        
        self.material = mat

    cdef bint move(self, STATE& state, double SP):
        """Performs particle transport.
        """
        
        raise RuntimeError("`move` called from virtual.")

    cdef bint is_inside(self, double3& pos):
        """Checks if a given position `pos` is contained by the volume.
        """
        
        raise RuntimeError("`is_inside` called from virtual in Volume")

    cdef double SDF(self, double3 pos):
        """Evaluate the Signed Distance Function (SDF) at position `pos`. 
        """

        raise RuntimeError("SDF FROM VOL WAS CALLED")
        
        
    cdef void depositLOCAL(self, double3& pos, double E):
        """Deposit the energy `E` locally in position `pos`. 
        
        Note: 
            The meaning of locality depends on the particulars of the geometry implentation. For
            solid volumes in the CSG implementation, the energy will be deposited in the entire
            volume. This will be different if the volume in question is a tally.
            
            This method is supposed to be overwritten. If called, it will raise a RuntimeError.
        """
        
        raise RuntimeError("'depositLOCAL' called from its virtual in 'Volume' ")

    cdef void depositUNIFORM(self, STATE& state, double SP):
        """Deposit the energy uniformly along its step lenght.
        """
        
        
        raise RuntimeError("depositUNIFORM called from Volume (virtual)")
        import time
        time.sleep(10_000)
    
    cdef void exit(self):
        """A method for reseting state variables in a convenient manner. 
        """
        
        raise RuntimeError("EXIT CALLED FROM VIRTUAL 'Volume' ")


    cdef void depositRANDOM(self, STATE& state, double E, double tau):
        raise RuntimeError("depositRANDOM CALLED FROM VIRTUAL 'Volume' ")

    def __gt__(self, other):
        """Deprecated.
        """
        self.tr.b.x = -other[0]
        self.tr.b.y = -other[1]
        self.tr.b.z = -other[2]
        
        self.mesh.translate(other)
        return True
        
    def __lshift__(self, material):
        """Deprecated.
        """
        self.fill(material)
    
    def __enter__(self): 
        return self

    cdef void depositLocaly(self, double3 pos, double E):
        """Deprecated.
        """
        pass


    def set_array(self, arr):
        """Deprecated.
        """
        self.voxels = arr
    
    def get_arr(self):
        """Deprecated.
        """
        return self.voxels
    

    

        





cdef class BVH(Volume):

	def __init__(self, *args, **kwargs):
		super(BVH, self).__init__(material = kwargs['material'])

		self.Nws = len(args) + 1
		self.ws = <void**> malloc(self.Nws * sizeof(void*))

		self.ws[0] = <void*> self
		for i, volume in enumerate(args):
			self.ws[i+1] = <void*> volume
			(<BVH> volume).set_outer(self, i + 1)

		if kwargs['render'] == True:
			@plt_geo.sdf3
			def this():
				def SDF(double[:,:] P):
					cdef double3 p
					cdef int N = len(P)
					cdef cnp.ndarray sd = np.zeros(N)
					cdef int i
					for i in range(N):
						p.x = P[i, 0]
						p.y = P[i, 1]
						p.z = P[i, 2]
						sd[i] = self.SDF(p)
					return sd
				return SDF

			generator = this()
			generator.save(f"geo/{kwargs['name']}.stl")


		self.cache = False
		self.lock = True



	cpdef set_outer(self, BVH other, int index):
		"""
		other -> outer volume
		index -> self's position in outers workspace
		"""
		self.outer = other
		self.position_in_outer = index


	def __contains__(self, other):
		cdef double3 pos
		pos.x, pos.y, pos.z = other
		return self.is_inside(pos)


	def get_mesh(self):
		import pyvista as pv
		return pv.read(f"geo/{self.name}.stl")

	def plot(self):
		import pyvista as pv
		mesh = pv.read(f"geo/{self.name}.stl")
		mesh.plot()


	cdef void* searchO(self, STATE& state):
		cdef int i

		for i in range(1, self.Nws):
			if self.ws[i] == state.current_region: continue

			if (<BVH> self.ws[i]).is_inside(state.pos):
				IF DEBUG_MODE: print(i)
				return self.ws[i]
		IF DEBUG_MODE: print(0)
		return <void*> self

	cdef void exit(self):
		cdef int i
		for i in range(self.Nws):
			#if (<BVH> self.ws[i]).keep: continue
			(<BVH> self.ws[i]).cache = False


	cdef bint move(self, STATE& state, double SP):
		raise RuntimeError("'move' called from its virtual in 'Volume.BVH' ")

	cdef void depositUNIFORM(self, STATE& state, double SP):
		raise RuntimeError("depositUNIFORM called from BVH (virtual)")
		print("depositUNIFORM called from BVH (virtual)")
		import time
		time.sleep(10_000)


	cdef void depositDISCRETE(self, STATE& state):
		raise RuntimeError("'depositDISCRETE' called from its virtual in 'Volume.BVH' ")
	
	cdef void depositLOCAL(self, double3& pos, double E):
		raise RuntimeError("'depositLOCAL' called from its virtual in 'Volume.BVH' ")

	cdef void depositRANDOM(self, STATE& state, double E, double tau):
		raise RuntimeError("'depositRANDOM' called from its virtual in 'Volume.BVH' ")

	cdef double main_intersect(self, STATE& state):
		raise RuntimeError("'main_intersect' called from its virtual in 'Volume.BVH' ")

	cdef void localSDF(self, STATE& state):
		raise RuntimeError("'localSDF' called from its virtual in 'Volume.BVH' ")


	cdef double SDF(self, double3 pos):
		raise RuntimeError("'SDF' called from its virtual in 'Volume.BVH' ")

	cdef bint is_inside(self, double3 pos):
		raise RuntimeError("'is_inside' called from its virtual in 'Volume.BVH' ")





