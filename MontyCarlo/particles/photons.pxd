




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
                                             










#Error messages (to be moved to its own module)
errorMSG1 = "Exhausted allowed number of iterations for rejection sampling."

#External Imports
#from numpy import *
#from numpy.random import rand, randint
#import pickle -> probly not needed any more?



from ..external.mixmax_interface cimport mixmax_engine


#Local Imports
from .particle import StopSimulation
from .particle cimport Particle
from ..tools.vectors cimport Vector




from ..materials import database as db
from ..materials.materials cimport Molecule, Atom, Shell

# --  -- from . import electrons as e
from libc.math cimport sin, cos, log, sqrt, pi 
from ..geometry.main cimport Volume


cdef struct IFMPcumul:
    long double C0, C1, C2, C3, C4, C5
    
cdef IFMPcumul IMFP_CUMUL
from ..materials.materials cimport Material
from ..materials.photon.photon cimport Coherent, Incoherent, Pairproduction, Tripletproduction
from ..materials.photon.photon cimport Photon as MPhoton

from .particle cimport STATE



#from ..materials.photon.photons cimport Electron as MATelectron

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
    
    cdef int N_coh, N_incoh, N_photo, N_pair, N_trip


    cdef void* current_material
    cdef void* coherent
    cdef void* incoherent
    cdef void* pairproduction
    cdef void* tripletproduction

    #cdef Material current_material
    #cdef Coherent coherent
    #cdef Incoherent incoherent
    #cdef Pairproduction pairproduction
    #cdef Tripletproduction tripletproduction
    cdef object S

    #cdef vector[double] ZZ
    cdef IFMPcumul IMFP_CUMUL
    cdef void* current_molecule



    ### CONSTRUCTORS
    @staticmethod
    cdef Photon _new(STATE& state)

    @staticmethod
    cdef Photon _newISOTROPIC(STATE& state)
    
    
    ### RUN SIMULATION
    cdef void _run(Photon self, mixmax_engine* genPTR)
    # cdef void _runFORCEcompton(self)
    # cdef void _runSPLIT(self)
    



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






