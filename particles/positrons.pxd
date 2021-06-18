# cython: annotate=True

from ..materials.materials cimport Atom as _Atom




#          _____          
#         /\    \         
#        /::\    \        
#       /::::\    \       
#      /::::::\    \      
#     /:::/\:::\    \     
#    /:::/__\:::\    \    
#   /::::\   \:::\    \   
#  /::::::\   \:::\    \  
# /:::/\:::\   \:::\    \ 
#/:::/__\:::\   \:::\____\
#\:::\   \:::\   \::/    /
# \:::\   \:::\   \/____/ 
#  \:::\   \:::\    \     
#   \:::\   \:::\____\    
#    \:::\   \::/    /    
#     \:::\   \/____/     
#      \:::\    \         
#       \:::\____\        
#        \::/    /        
#         \/____/ 











#Error messages (to be moved to its own module)
errorMSG1 = "Exhausted allowed number of iterations for rejection sampling."

from .._random.mixmax.interface cimport mixmax_engine

## PYTHON IMPORTS
#Local Imports

from .particle cimport STATE
from ..materials import database as db
from ..settings import __photonCUTOFF__

from .electrons cimport Electron
from .photons cimport Photon
## CYTHON IMPORTS
#Local Imports
from .particle cimport Particle
from ..geometry.main cimport Volume
from ..tools.vectors cimport Vector
from ..tools.interpol1 cimport hLinLinInterpolation
from ..materials.materials cimport Material
from ..materials.positron.main cimport Brem, Inelastic, Elastic, Anihilation
from ..materials.positron.main cimport Positron as MATpositron

from ..materials.positron.GOS cimport CMolecule
from ..materials.positron.GOSfinal cimport gosMolecule


# External Imports
from libc.math cimport sin, cos, log, sqrt, pi, exp, fmin, fmax
cdef struct IFMPcumul:
    long double C0, C1, C2, C3, C4, C5


cdef class Positron(Particle):
    #cdef CMolecule GOS
    cdef gosMolecule GOS
    cdef double cos, Esec, cos_sec
    cdef IFMPcumul IMFP_CUMUL
    cdef double avgW, varW
    cdef hLinLinInterpolation _imfp
    cdef MATpositron positron 
    cdef Material current_material
    cdef Elastic elastic
    cdef Inelastic inelastic
    cdef Brem brem
    cdef Anihilation anih
    cdef double imfp_max, SP, STRAGG
    cdef double T1, T2, imfp0

    cdef double s, s_max, w, mu
    cdef double rc
    
    # @staticmethod #custom constructor for speed
    # cdef Electron _new(Volume current_region,
    #                    double E, 
    #                    Vector pos, 
    #                    Vector ey,
    #                    Vector ez,
    #                    double s_max)
    
    @staticmethod
    cdef Positron _new(STATE& state)
    
    @staticmethod
    cdef Positron _newISOTROPIC(STATE& state)
    
    
    
    # @staticmethod #custom constructor for speed
    # cdef Electron _newISOTROPIC(Volume current_region,
    #                    double E, 
    #                    Vector pos)
        

    cdef void _run(self, mixmax_engine* genPTR)
    
    cdef void _anihilation(self)
    cdef void update_references(self) 



    cdef void update_imfp_cumul(self)
       
        
    cdef void update_imfp(self) 

        
        
        


    
    
    
    cdef inline void sample_w(self, double tau)
    
    
    cdef inline void do_hinge(self)
        
    
        
        
    
    
    
    cdef inline void _elastic(self) 
        
        
    cdef inline void _brem(self) 
        
    cdef inline void _inelastic(self) 
    cdef inline void _delta(self)
    
    cdef inline int find_index(self)
    
    