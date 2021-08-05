# cython: annotate=False


# remove ... 
from ..materials import database as db
from ..settings import __photonCUTOFF__


from ..external.mixmax_interface cimport mixmax_engine

from .particle  cimport Particle
from .electrons cimport Electron
from .photons   cimport Photon

from ..types cimport STATE
from ..geometry.main   cimport Volume
from ..tools.vectors   cimport Vector
from ..tools.interpol1 cimport hLinLinInterpolation
from ..materials.materials cimport Material
from ..materials.materials cimport Atom as _Atom
from ..materials.positron.main     cimport Brem
from ..materials.positron.main     cimport Inelastic
from ..materials.positron.main     cimport Elastic
from ..materials.positron.main     cimport Anihilation
from ..materials.positron.main     cimport Positron as MATpositron
from ..materials.positron.GOS      cimport CMolecule
from ..materials.positron.GOSfinal cimport gosMolecule

from libc.math cimport sin
from libc.math cimport cos
from libc.math cimport log
from libc.math cimport sqrt
from libc.math cimport pi
from libc.math cimport exp
from libc.math cimport fmin
from libc.math cimport fmax

cdef struct IFMPcumul:
    long double C0
    long double C1
    long double C2
    long double C3
    long double C4
    long double C5

ctypedef Material MAT

cdef class Positron(Particle):
    #cdef CMolecule GOS
    cdef gosMolecule GOS
    cdef IFMPcumul IMFP_CUMUL

    cdef double cos
    cdef double Esec
    cdef double cos_sec
    cdef double avgW
    cdef double varW

    cdef hLinLinInterpolation _imfp

    cdef MATpositron positron 

    cdef Material current_material

    cdef Elastic     elastic
    cdef Inelastic   inelastic
    cdef Brem        brem
    cdef Anihilation anih

    cdef double imfp_max
    cdef double SP
    cdef double STRAGG
    cdef double T1
    cdef double T2
    cdef double imfp0

    cdef double s
    cdef double s_max
    cdef double w
    cdef double mu
    cdef double rc
    

    cdef void _run(self, mixmax_engine* genPTR)

    # Constructors
    @staticmethod
    cdef Positron _new(STATE& state)
    
    @staticmethod
    cdef Positron _newISOTROPIC(STATE& state)

    # Updates
    cdef void update_references(self) 
    cdef void update_imfp_cumul(self)
    cdef void update_imfp(self) 


    # Interactions
    cdef void _anihilation(self)
    cdef inline void _elastic(self) 
    cdef inline void _brem(self) 
    cdef inline void _inelastic(self) 
    cdef inline void _delta(self)

    # Condensed History
    cdef inline void sample_w(self, double tau)
    cdef inline void do_hinge(self)

    # Utils
    cdef inline int find_index(self)
    
    