


#          /\    \                 /\    \         
#         /::\    \               /::\____\        
#        /::::\    \             /:::/    /        
#       /::::::\    \           /:::/    /         
#      /:::/\:::\    \         /:::/    /          
#     /:::/__\:::\    \       /:::/____/           
#    /::::\   \:::\    \     /::::\    \           
#   /::::::\   \:::\    \   /::::::\    \   _____  
#  /:::/\:::\   \:::\____\ /:::/\:::\    \ /\    \ 
# /:::/  \:::\   \:::|    /:::/  \:::\    /::\____\
# \::/    \:::\  /:::|____\::/    \:::\  /:::/    /
#  \/_____/\:::\/:::/    / \/____/ \:::\/:::/    / 
#           \::::::/    /           \::::::/    /  
#            \::::/    /             \::::/    /   
#             \::/____/              /:::/    /    
#              ~~                   /:::/    /     
#                                  /:::/    /      
#                                 /:::/    /       
#                                 \::/    /        
#                                  \/____/         
                                             


#Local Imports
from ..materials import database as db

#from .particle cimport STATE
from .particle cimport Particle
from ..types cimport STATE
from ..types cimport PySTATE
from ..tools.vectors cimport Vector
from ..external.mixmax_interface cimport mixmax_engine
from ..geometry.main cimport Volume
from ..materials.materials cimport Molecule
from ..materials.materials cimport Atom
from ..materials.materials cimport Shell
from ..materials.materials cimport Material
from ..materials.photon.photon cimport Coherent
from ..materials.photon.photon cimport Incoherent
from ..materials.photon.photon cimport Pairproduction
from ..materials.photon.photon cimport Tripletproduction
from ..materials.photon.photon cimport Photon as MPhoton



#External Imports
from libc.math cimport sin
from libc.math cimport cos
from libc.math cimport log
from libc.math cimport sqrt
from libc.math cimport pi 

cdef struct IFMPcumul:
    long double C0
    long double C1
    long double C2
    long double C3
    long double C4
    long double C5

cdef IFMPcumul IMFP_CUMUL

ctypedef Coherent Coh
ctypedef Incoherent inCoh
ctypedef Pairproduction PP
ctypedef Tripletproduction PPP
ctypedef Molecule Mol
ctypedef MPhoton Ph
ctypedef Volume V
ctypedef Material Mat
ctypedef Material M

#cimport cython

cdef class Photon(Particle):
    cdef double k
    
    cdef int N_coh
    cdef int N_incoh
    cdef int N_photo
    cdef int N_pair
    cdef int N_trip

    cdef void* current_material
    cdef void* coherent
    cdef void* incoherent
    cdef void* pairproduction
    cdef void* tripletproduction

    cdef object S

    cdef IFMPcumul IMFP_CUMUL
    cdef void* current_molecule



    ### CONSTRUCTORS
    @staticmethod
    cdef Photon _new(STATE& state)

    @staticmethod
    cdef Photon _newISOTROPIC(STATE& state)
    
    cdef void _run(Photon self, mixmax_engine* genPTR)

    cdef void record(self)
    cdef inline int find_index(self)
    
    #### UPDATE METHODS
    cdef void update_references(self)
    cdef void update_imfp(Photon self)
    
    
    #### INTERACTION METHODS
    cdef void _coherent(Photon self)
    cdef void _incoherent(Photon self)
    cdef void _pairproduction(Photon self) 
    cdef void _photoelectric(Photon self)
    cdef void _tripletproduction(Photon self)