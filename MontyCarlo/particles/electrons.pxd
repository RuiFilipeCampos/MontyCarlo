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

from ..external.mixmax_interface cimport mixmax_engine

## PYTHON IMPORTS
#Local Imports
from .particle import StopSimulation
from ..materials import database as db
from ..settings import __photonCUTOFF__


## CYTHON IMPORTS
#Local Imports
from .particle cimport Particle
from ..geometry.main cimport Volume
from ..tools.vectors cimport Vector
from ..tools.interpol1 cimport hLinLinInterpolation
from ..materials.materials cimport Material
from ..materials.electron.main cimport Brem, Inelastic, Elastic
from ..materials.electron.main cimport Electron as MATelectron
from ..materials.electron.GOS cimport CMolecule
from ..materials.electron.GOSfinal cimport gosMolecule


# External Imports
from libc.math cimport sin, cos, log, sqrt, pi, exp, fmin, fmax
cdef struct IFMPcumul:
    long double C0, C1, C2, C3, C4, C5

from .particle cimport STATE


cdef class Electron(Particle):
    cdef double wmax
    cdef object MU
    #cdef CMolecule GOS
    cdef gosMolecule GOS
    cdef double cos, Esec, cos_sec
    cdef IFMPcumul IMFP_CUMUL
    cdef double avgW, varW
    cdef hLinLinInterpolation _imfp
    cdef MATelectron electron 
    cdef Material current_material
    cdef Elastic elastic
    cdef Inelastic inelastic
    cdef Brem brem
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
    cdef Electron _new(STATE& state)
    
    @staticmethod
    cdef Electron _newISOTROPIC(STATE& state)
    
    
    
    # @staticmethod #custom constructor for speed
    # cdef Electron _newISOTROPIC(Volume current_region,
    #                    double E, 
    #                    Vector pos)
        

    cdef void _run(Electron self, mixmax_engine* genPTR)
    

    cdef void update_references(self) 



    cdef void update_imfp_cumul(Electron self)
       
        
    cdef void update_imfp(Electron self) 

        
        
        


    
    
    
    cdef inline void sample_w(self, double tau)
    
    
    cdef inline void do_hinge(self)
        
    
        
        
    
    
    
    cdef inline void _elastic(self) 
        
        
    cdef inline void _brem(self) 
        
    cdef inline void _inelastic(self) 
    cdef inline void _delta(self)
    
    cdef inline int find_index(self)
    
    