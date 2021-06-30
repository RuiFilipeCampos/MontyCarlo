# distutils: language = c++
# cython: profile=False

print(">>>>>   IMPORTING GEOMETRY")

from libc.math cimport fmin, fmax, sqrt , cos, pi
from libcpp.vector cimport vector

from numpy.math cimport INFINITY
from ..particles.particle cimport STATE

#from ..tools.vectors cimport Vector

import numpy as np 

#from ..random.mixmax.interface cimport mixmax_engine

from numpy.random import randint


from numpy cimport ndarray


from libc.stdlib cimport malloc, free
cdef double MINdr = 0

from libcpp.list cimport list as cpplist
    
from cython.operator import dereference, preincrement



from ..types cimport double3
    
cdef struct Ray:
    double3 O
    double3 D

cdef struct intersection:
    double3 pos
    double d



cdef class Volume:

    cdef void depositLOCAL(self, double3& pos, double E):
        raise RuntimeError("'depositLOCAL' called from its virtual in 'Volume' ")

    cdef void depositUNIFORM(self, STATE& state, double SP):
        raise RuntimeError("depositUNIFORM called from Volume (virtual)")
        import time
        time.sleep(10_000)
    
    cdef void exit(self):
        raise RuntimeError("EXIT CALLED FROM VIRTUAL 'Volume' ")


    cdef void depositRANDOM(self, STATE& state, double E, double tau):
        raise RuntimeError("depositRANDOM CALLED FROM VIRTUAL 'Volume' ")

    def __gt__(self, other):
        self.tr.b.x = -other[0]
        self.tr.b.y = -other[1]
        self.tr.b.z = -other[2]
        
        self.mesh.translate(other)
        return True
        
    def __lshift__(self, material):
        self.fill(material)
    
    def __enter__(self): return self

    cdef void depositLocaly(self, double3 pos, double E):
        pass



    def set_array(self, arr):
        self.voxels = arr
    
    def get_arr(self):
        return self.voxels
    

    
    def fill(self, mat):
        self.material = mat 
        
    cdef bint move(self, STATE& state, double SP):
        raise RuntimeError("move called from virtual ")

    cdef bint is_inside(self, double3& pos):
        raise RuntimeError("is_inside called from virtual in Volume")

    cdef double SDF(self, double3& pos):
        raise RuntimeError("SDF FROM VOL WAS CALLED")