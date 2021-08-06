# cython: annotate = False
# cython: profile = False
# distutils: language = c++


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
                                             


print("Importing `.particles.photons`")



# Conditional Compilation for debugging.
DEF TEST = True
DEF _DEBUG_BASIC = False
DEF _DEBUG = False
DEF _DEBUGincoh = False
DEF RECORD = True

DEF _COH = True
DEF _INCOH = True
DEF _PP = True
DEF _TP = True
DEF _PH = True


# Internal Imports
from .._init import eax
from ..materials import database as db
from ..settings import __photonCUTOFF__
from ..settings import __electronCUTOFF__
from ..settings import DEBUG


#from ..materials.electron.main import eax as _eax
from ..external.mixmax_interface cimport mixmax_engine
from .._init cimport EAX
from .._init cimport LIMS
from .particle  cimport Particle
from .particle  cimport STATE
from .electrons cimport Electron
from .positrons cimport Positron
from ..geometry.main cimport Volume
from ..materials.pyRelax cimport Atom as RAtom
from ..materials.cppRelaxAPI cimport PARTICLES



# External Imports
from collections import deque
import numpy as np

from libcpp.vector cimport vector
from libc.math cimport sin
from libc.math cimport cos
from libc.math cimport log
from libc.math cimport sqrt
from libc.math cimport pi
from libc.math cimport acos
from libc.math cimport exp
cimport cython


cdef extern from "<math.h>" nogil:
    double frexp(double x, int* exponent)


#Error messages (to be moved to its own module)
errorMSG1 = "Exhausted allowed number of iterations for rejection sampling."



# CONSTANTS AND GLOBALS - this needs to be better for sure
cdef double Eel0_MeV = 0.510998950000
cdef double Eel0_eV = Eel0_MeV*1e6
cdef double k_cutoff = __photonCUTOFF__/Eel0_eV
#cdef double CUTOFF = __photonCUTOFF__
cdef double CUTOFFel = __electronCUTOFF__
cdef double photonCUTOFF = __photonCUTOFF__
cdef double electronCUTOFF = __electronCUTOFF__
cdef double minCUTOFF = min(photonCUTOFF, electronCUTOFF)

IMFP_CUMUL.C0 = 0.

cdef struct INCOHERENT:
    double t1
    double t2
    double tau_min
    double tau, cos, N, D, sin2, x , T
    double k




@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef class Photon(Particle):

    cdef double ENERGY(self):
        return self.k*Eel0_eV
    
    @staticmethod
    cdef Photon _new(STATE& state):
        cdef Photon self
        self = <Photon>Photon.__new__(Photon)
        self.state = state
        return self
    
    @staticmethod
    cdef Photon _newISOTROPIC(STATE& state):
        cdef Photon self
        self = <Photon>Photon.__new__(Photon)

        self.state = state

        cdef double a
        while True:
            x = 2*self.state.genPTR.get_next_float() - 1
            y = 2*self.state.genPTR.get_next_float() - 1

            a = x**2 + y**2

            if a < 1:
                break

        self.state.dire.x = 1 - 2*a 
        a = 2 * sqrt(1 - a)
        self.state.dire.y = x*a
        self.state.dire.z = y*a
        
        #azimuth is thrown in next interaction <- reconfirm
        return self


    ####################################################################################
    ########                           RUN                                      ########
    ########                           RUN                                      ########
    ####################################################################################

    
    cdef void _run(Photon self, mixmax_engine* genPTR):
        IF _DEBUG_BASIC: print("> PHOTON")

        #cdef double r
        self.secondary = deque()
        self.nSECONDARY = 0
        self.k = self.state.E/Eel0_eV
        if self.k < k_cutoff:
            (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
            return
        
        self.state.genPTR = genPTR

        self.update_references()

        IF RECORD: self.record()
        
        cdef double r
        cdef double L 

        IF _DEBUG: print("STARTED")
        
        #print(<M> self.current_material)
        while True:


            self.state.L = -log(1e-9 + (1-1e-9)*self.state.genPTR.get_next_float())/self.imfp_T
            #if self.state.L > 1e3: print("PHOTONS")


            if (<V> self.state.current_region).move(self.state, 0.):
                if (<V> self.state.current_region).opaque:
                    (<V> self.state.current_region).exit()
                    (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
                    return


                self.update_references()
                #print(<M> self.current_material)

                IF RECORD: self.record()

                continue
            
            if self.state.pos.x**2 + self.state.pos.y**2 + self.state.pos.z**2 > 10_000**2:
                return
            #    print("\n\n\n\n\n\n")
            #    print("*****")
#
            #    import time
            #    time.sleep(10000)
            
            IF RECORD: self.record()
            
            r = self.state.genPTR.get_next_float()*self.IMFP_CUMUL.C5



            if   r < self.IMFP_CUMUL.C1:
                self._coherent()

            elif r < self.IMFP_CUMUL.C2: 
                self._incoherent()

                if self.k < k_cutoff:
                    (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
                    (<V> self.state.current_region).exit()
                    return

                self.update_imfp()

            elif r < self.IMFP_CUMUL.C3:
                self._pairproduction()
                (<V> self.state.current_region).exit()
                return

            elif r < self.IMFP_CUMUL.C4:
                self._tripletproduction()
                (<V> self.state.current_region).exit()
                return

            else:     
                self._photoelectric()
                (<V> self.state.current_region).exit()
                return




    ####################################################################################
    ########                           UPDATE                                   ########
    ########                           METHODS                                  ########
    ####################################################################################

    cdef void update_references(self):
        """
        Updates all references. Called when there is a region crossing.
        """
        
        #getting material from current region
        cdef void* handler = self.state.current_region

        #input("current_mat")
        self.current_material =  <void*> (<V> handler).material
        #input("current_molecule")
        self.current_molecule =  <void*> (<M> (<V> handler).material).molecule
        #input("photon")
        cdef void* photon = <void*> (<M> self.current_material).photon
        #input("others")
        #these references are used by their corresponding _fooInteraction method
        self.coherent          = <void*>  (<Ph> photon).coherent
        #input("others")
        self.incoherent        = <void*>  (<Ph> photon).incoherent
        #self.photoelectric     = self.current_material.photoelectric
        #input("others")
        self.pairproduction    = <void*>  (<Ph> photon).pairproduction
        #input("others")
        self.tripletproduction = <void*>  (<Ph> photon).tripletproduction
        #input("all of them worked")

        #self.S = self.incoherent.S

        #since region crossing has ocurred, update the inverse mean free paths
        self.update_imfp()


    cdef void update_imfp(Photon self):
        """
        Updates inverse mean free paths. Constructs cumul.
        Called when there is a region crossing or energy of photon has changed.
        """

        self.state.E = self.k * Eel0_eV
        cdef int i = self.find_index()
        #self.state.Energy.push_back(self.state.E)
        #IMFP_CUMUL.C0 = 0.

        #print("energy:", self.state.E)
        self.IMFP_CUMUL.C1 = (<Coh> self.coherent).imfpA[i] + self.state.E*(<Coh> self.coherent).imfpB[i]
        #print("coh:", self.IMFP_CUMUL.C1)
        self.IMFP_CUMUL.C2 = self.IMFP_CUMUL.C1 + (<inCoh> self.incoherent).imfpA[i]      + self.state.E*(<inCoh> self.incoherent).imfpB[i]
        #print("incoh", (<inCoh> self.incoherent).imfpA[i]      + self.state.E*(<inCoh> self.incoherent).imfpB[i])
        self.IMFP_CUMUL.C3 = self.IMFP_CUMUL.C2 + (<PP> self.pairproduction).imfpA[i]     + self.state.E*(<PP> self.pairproduction).imfpB[i]
        #print("pp: ", (<PP> self.pairproduction).imfpA[i]     + self.state.E*(<PP> self.pairproduction).imfpB[i])
        self.IMFP_CUMUL.C4 = self.IMFP_CUMUL.C3 + (<PPP> self.tripletproduction).imfpA[i] + self.state.E*(<PPP> self.tripletproduction).imfpB[i]
        #print("ppp:", (<PPP> self.tripletproduction).imfpA[i] + self.state.E*(<PPP> self.tripletproduction).imfpB[i])
        self.IMFP_CUMUL.C5 = self.IMFP_CUMUL.C4 + (<Mol> self.current_molecule).PHELa[i]  + self.state.E*(<Mol> self.current_molecule).PHELb[i]
        #print("photo", (<Mol> self.current_molecule).PHELa[i]  + self.state.E*(<Mol> self.current_molecule).PHELb[i])
        #a = input("continue?")
        self.imfp_T = self.IMFP_CUMUL.C5
        #print("imfp_cumul WORKING")


    ####################################################################################
    ########                          INTERACTION                               ########
    ########                           SAMPLERS                                 ########
    ####################################################################################


    cdef void _coherent(Photon self):
        """
        Rayleigh sampling!
        """
        IF not _COH: return
        IF _DEBUG: print("(( ._coherent")

        #self.N_coh += 1
        
        
        
        cdef double k2 = self.k*self.k
        cdef double qmax2 = 2*k2
        cdef double cumulMAX = (<Coh> self.coherent).evalY(qmax2)
        
        cdef double  x2, cos
        
        cdef double r
        
        #print(cumulMAX)
        
        while 1:
            #r = self.state.dire.x.get_next_float()*cumulMAX
            r = self.state.genPTR.get_next_float()*cumulMAX
            #print(r)
            x2 = (<Coh> self.coherent).evalX(r)
            cos = 1 - x2/k2 #x2/k2 = 1 - cos
            if self.state.genPTR.get_next_float()*2 < 1 + cos*cos:
                break

        self.throwAZIMUTH()
        self.rotateTHETA(cos)
        
       #  cdef double r, x2, cos
       #  #cdef double x_max = 20.6074*2*self.k
       #  cdef double x_max = 41.2148*self.k
       #  cdef double r_max = self.coherent.FF.cumul(x_max) #internals of this needs work
       # # print(x_max, x_max**2)
         
       #  while 1:
            
       #      #Sample x**2 from squared form factor (limited in (0, x_max**2))
       #      r  = self.genPTR.get_next_float()*r_max

       #      x2 = self.coherent.FF.invCum(r)
       #      #print(self.coherent.FF.invCum(r))

       #      #Get cos(theta) from x**2 and k = E/0.511MeV
       #      cos  = 1 - 0.5 * x2 / (20.6074*self.k)**2
       #      #cos  = 1 - 0.0011773988362909597 * x2 / self.k**2

       #      #Sample from thomson scattering.
       #      if 2*self.genPTR.get_next_float() < 1+cos**2:
       #          #print(">>>ray", cos)
       #          self.throwAZIMUTH()
       #          self.rotateTHETA(cos)
       #          #if cos < -1: print("cos", cos)
       #          #if cos > 1: print("cos", cos)
                
       #          #self.change_direction(cos, 2*pi*self.genPTR.get_next_float())
       #          IF _DEBUG: print(" ._coherent  ))")
       #          break

    

    # cdef void _incoherentFREE(self):
    #     pass
    
    # cdef void _incoherentPENELOPE(self):
    #     pass
    
    # cdef void _incoherentLIVERMORE(self):
    #     pass

    cdef void _incoherent(Photon self):
        """
        Compton sampling.
        
        This sampling of this interaction consists in ? steps:
            

            
            (1) Select an atom at random. The probability of choosing the ith atom 
            is proportional to ni*Zi. Where,
                
                ni -> number of atoms Zi in molecule
                Zi -> number of electrons in the atom
            
            (2) Sample from cos(theta) from Klein-Nishina differential cross 
            section. Reject results based on the scattering function of the
            chosen atom.
            
            (3) Select a shell based on occupancy numbers. Reject shells whose
            binding energy is greater than the photons current energy.
            
            (4) Sample doppler broadning according. 
            
                (4.1) Sample the targeted electrons momentum p_z from the
                compton profile.
                
                (4.2) 
                
        
        """
        IF not _INCOH: return
        IF TEST:
            cdef double _E = self.state.E
        
        IF _DEBUG: print("(( ._incoherent")
        IF _DEBUG: print(f"STARTING: N = {self.N_incoh} | k =  {self.k} | E = {self.k*Eel0_eV}")
        
        
        # (1) Choose an atom with probabilities based on Z
     #   print("-------START COMPTON")

        cdef Atom active_atom = (<Mol> self.current_molecule).choose_atom(self.state.genPTR)
        cdef double Uk

        # (2) Sample cos(theta) from KN DCS
        
        cdef double eps_min = 1/(1 + 2*self.k) 
        cdef double alpha1 = log(1 + 2*self.k)
        cdef double alpha2 = (1 - eps_min*eps_min)/2
        cdef double eps, oneMINUScos
        cdef double g, cos, sin2
        cdef double gmax = eps_min + 1/eps_min
        
        while 1:
            if (alpha1 + alpha2)*self.state.genPTR.get_next_float() <= alpha1:
                eps = eps_min*exp(alpha1 * self.state.genPTR.get_next_float())
            else:
                eps = sqrt(eps_min*eps_min + 2*alpha2*self.state.genPTR.get_next_float() )
    
            oneMINUScos = (1/eps - 1)/self.k
            sin2 = 1 - (1 - oneMINUScos)**2
            
            #g = (1/eps + eps - sin2)/gmax
            g = 1 - eps*sin2/(1 + eps*eps)
            if self.state.genPTR.get_next_float() < g*active_atom.S._eval(oneMINUScos*(self.k)**2)  :
                break
            
            
        cos = (1 - oneMINUScos)
        # print(eps_min)
        # print(eps_min)
        # print(eps_min)
        # print(eps)
        # print("")
            
            
            
            
        #print(active_atom.S._eval(oneMINUScos*(eps*self.k)**2))
     #   print("g", g)
        cdef Shell active_shell
        
        cdef double R
        cdef int N
        cdef double pz
        cdef double proposed_k
       # print(self.k)
       # print("COS", cos)
        
        
        cdef Electron el
        cdef double new_k, Eel
        
        cdef int count1
        for count1 in range(1_000):

            
            ### CHOOSE SHELL
            while 1:
                R = self.state.genPTR.get_next_float()*active_atom.Nsh
                N = <int> R
                if R - N < active_atom.ALIAS[N, 1]:
                    active_shell = active_atom.arrSHELLS[<int> active_atom.ALIAS[N, 0]]
                    if active_shell.binding_energy < self.k*Eel0_eV: break
                    else: continue
                
                active_shell = active_atom.arrSHELLS[<int> active_atom.ALIAS[N, 2]]
                
                if active_shell.binding_energy < self.k*Eel0_eV: break
                else: continue
                
            ### SAMPLE ELECTRON MOMENTUM
            Uk = active_shell.binding_energy/Eel0_eV
            alpha1 = self.k*(self.k - Uk)*oneMINUScos #reusing declared double 
            pz = active_shell.sample_compton_profile(self.state.genPTR, 
                                                     (alpha1 - Uk)/sqrt(2*alpha1 - Uk*Uk) #= pz_max
                                                     )
            
            
            proposed_k = (1 - pz*pz*eps*cos + pz*sqrt( 1 - 2*eps*cos + eps*eps*(1 - pz*pz*sin2))  )/(1 - pz*pz*eps*eps)

            if 0 < proposed_k*eps < self.state.genPTR.get_next_float():
                new_k = proposed_k*eps*self.k
                break
    
        else: new_k = eps*self.k

        
        Eel = (self.k - new_k - Uk)*Eel0_eV
        
        self.throwAZIMUTH()
        
        
        if Eel > electronCUTOFF:

            el = Electron._new(self.state)

            el.state.E = Eel
            el.state.axis.x *= -1
            el.state.axis.y *= -1
            el.state.axis.z *= -1

            el.rotateTHETA((self.k - new_k*cos)/sqrt(new_k*new_k + self.k*self.k - 2*new_k*self.k*cos))
            self.nSECONDARY += 1
            self.secondary.append(el)
        else:
            (<V> self.state.current_region).depositLOCAL(self.state.pos, Eel)
            
        self.k = new_k
        self.state.E = self.k*Eel0_eV

        self.rotateTHETA(cos)
        #print("comp_cos ", cos)
        
        
        if active_shell.binding_energy < minCUTOFF:
            (<V> self.state.current_region).depositLOCAL(self.state.pos, active_shell.binding_energy)
            return
        
        ## relaxation
        cdef PARTICLES particles
        active_atom.ionize(active_shell.index, self.state.genPTR, &particles)
        
        
        cdef double E

        cdef double Etot = active_shell.binding_energy
        
        cdef Photon ph
        for i in range(particles.PHOTONS.size()):
            
            E = particles.PHOTONS.back()
            particles.PHOTONS.pop_back()
            if E < photonCUTOFF:
                continue
            
            
            ph = Photon._newISOTROPIC(self.state)
            ph.state.E = E
            Etot -= E

            self.nSECONDARY += 1
            self.secondary.append(ph)
        

        for i in range(particles.ELECTRONS.size()):
            E = particles.ELECTRONS.back()
            particles.ELECTRONS.pop_back()
            if E < CUTOFFel: continue
            el = Electron._newISOTROPIC(self.state)
            el.state.E = E
            Etot -= E
            self.nSECONDARY += 1
            self.secondary.append(el)
        
        (<V> self.state.current_region).depositLOCAL(self.state.pos, Etot)
        
        return




        
    
        
        # cdef Electron p
        # cdef INCOHERENT INCOH
        # INCOH.k = self.k
        # INCOH.t1 = log(1. + 2.*INCOH.k)
        # INCOH.t2 = INCOH.t1 + 2.*INCOH.k*(1. + INCOH.k)/(1. + 2.*INCOH.k)**2
        # #INCOH.t2 = 2.*INCOH.k*(1. + INCOH.k)/(1. + 2.*INCOH.k)**2
        # #INCOH.t1 = INCOH.t1/(INCOH.t1 + INCOH.t2)
        # INCOH.tau_min = 1./(1. + 2.*INCOH.k)

        # cdef double r
        # cdef double cos, E_new, tausq
        # tausq = INCOH.tau_min*INCOH.tau_min
        # #cdef int _
        # IF _DEBUGincoh: 
        #     print(f"Entering loop in 3 seconds")
        #     import time
        #     time.sleep(3)
        
        # cdef double S
        # cdef bint COND
        # #cdef double i= 0
        # while True:

        #     #i += 1
        #     #Sample fractional energy loss(tau)
        #     #if self.genPTR.get_next_float() < INCOH.t1/(INCOH.t1 + INCOH.t2):
                
        #    # COND = self.genPTR.get_next_float()*INCOH.t2 < INCOH.t1
            
            
        #    # r = self.genPTR.get_next_float()
        #     #INCOH.tau  = COND*(INCOH.tau_min**r) + (not COND)*sqrt(tausq + r*(1. - tausq))
                
        #     #print("choosing INCOH.tau")
        #     #time.sleep(2)
            
            
            
        #     if self.genPTR.get_next_float()*INCOH.t2 < INCOH.t1: 
        #         INCOH.tau = INCOH.tau_min**self.genPTR.get_next_float()
        #     else: 
        #         INCOH.tau = sqrt(tausq + self.genPTR.get_next_float()*(1. - tausq))

        #     #print("a bunch of calculations")
        #     #time.sleep(2)

        #     #Calculate cosine from current energy and tau
        #     #INCOH.cos = (INCOH.k + 1. - 1./INCOH.tau)/INCOH.k
        #     INCOH.cos = (1 - (1-INCOH.tau)/(INCOH.tau*INCOH.k))
        #     #constructing eqn in braces
        #     #INCOH.N = (1. - INCOH.tau)*( (2. * INCOH.k + 1.) * INCOH.tau - 1.)
        #     #INCOH.D = INCOH.k**2. * INCOH.tau * (1. + INCOH.tau**2.)

        #     #constructing argument for Incoherent Form Factor
        #     #INCOH.sin2 = sqrt(.5 * (1. - INCOH.cos))
        #     #INCOH.x = INCOH.k*INCOH.sin2*41.2148#2*20.6074
        #    # INCOH.x = INCOH.k*sqrt(.5 * (1. - INCOH.cos))*41.2148#2*20.6074
            
        #     #calculating proposal
        #     #INCOH.T = (1. - INCOH.N/INCOH.D)
        #     INCOH.T = (1. - (1. - INCOH.tau)*( (2. * INCOH.k + 1.) * INCOH.tau - 1.)/(INCOH.k**2. * INCOH.tau * (1. + INCOH.tau*INCOH.tau)))
        #     #S = SS(INCOH.k*sqrt(.5 * (1. - INCOH.cos))*41.2148)
            
            
            
            
            
        #     #INCOH.T *= SS(INCOH.k*sqrt(.5 * (1. - INCOH.cos))*41.2148)
            
            
        #     INCOH.T *= active_atom.S._eval(INCOH.k*sqrt(.5 * (1. - INCOH.cos))*41.2148)
            
            
        #     #print("accept ??")
        #     #time.sleep(2)
            
        #     #print("")
            

        #     if self.genPTR.get_next_float() < INCOH.T:

                
                
        #         IF _DEBUGincoh: print(f"ACCEPTED SAMPLE tau = {INCOH.tau}")
        #         #0.510998950
        #         #self.state.E = E0_el*self.k
                


        #         self.state.E = self.k*Eel0_eV
                
        #         self.k = INCOH.tau*INCOH.k
                
                
                
                
        #         #E_new = E0_el*self.k
        #         E_new = self.k*Eel0_eV
                
        #         cos = (self.state.E - E_new*INCOH.cos)/sqrt(self.state.E**2 + E_new**2 - 2*self.state.E*E_new*INCOH.cos)
                
                
        #         IF _DEBUGincoh: print(f"CHANGING DIRECTION: cos = {INCOH.cos} ")
        #         self.throwAZIMUTH()
                
            
        #         IF _DEBUGincoh: print(f"passed azimuth 1")
                
                
                
                
                
                
                
        #         E_new = (self.state.E - E_new)
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
        #         if E_new < CUTOFFel:
        #             self.rotateTHETA(INCOH.cos)
        #             return
                
                
        #         IF _DEBUGincoh: print("--- creating secondary ")
                
        #         p = Electron._new(E_new, self.x, self.y, self.z,
        #                                       self.state.Eyx, self.state.Eyy,self.state.Eyz,
        #                                       self.state.Ezx, self.state.Ezy,self.state.Ezz, 
        #                                       self.current_region)
                
        #         self.rotateTHETA(INCOH.cos)
                
        #         IF _DEBUGincoh: print("particle created")
        #         p.genPTR = self.genPTR
        #         p.throwAZIMUTH()
                
                
        #         IF _DEBUGincoh: 
        #             import time
        #             print(f"passed azimuth 2")
        #             time.sleep(1)
        #         p.rotateTHETA(cos)
        #         self.nSECONDARY += 1
        #         self.secondary.append(p)
                
                
                
        #         return
                
                
                
                
                
        #         # axis = ez0
        #         # ey = ey0.rotateAngle(axis, 2*pi*self.genPTR.get_next_float())
        
        #         # axis = ey
        #         # ez = ez0.rotateCos(axis, cos)
        #         # self.secondary.append(Electron._new(self.current_region,
        #         #               E_new, 
        #         #               self.pos,
        #         #               ey,
        #         #               ez,
        #         #               100))
                
                
                
        #         # self.nSECONDARY += 1
                

                
                
                    
        #         # self.change_direction(INCOH.cos, 2.*pi*self.genPTR.get_next_float())
                
                






    cdef void _pairproduction(Photon self):
        """
        TESTING:
            Conservation of energy: CONFIRMED
        """




        IF not _PP: return
        IF _DEBUG: print("(( ._pairproduction")
        
        #self.N_pair += 1
        
        
        #SAMPLING ENERGY OF POSITRON AND ELECTRON
        
        cdef double u1, u2, phiHalf_1, phiHalf_2, phi1, phi2
        cdef int i
        phiHalf_1, phiHalf_2 = (<PP> self.pairproduction).getPhis(.5, self.k)

        u1 = phiHalf_1 * (2/3) * (.5 - 1/self.k)**2
        u2 = phiHalf_2
        
        
        while 1:

            
            if True if (u1 + u2)*self.state.genPTR.get_next_float() < u1 else False:
                eps  = .5 + (.5 - 1/self.k)*(2*self.state.genPTR.get_next_float() - 1)**(1/3)
                phi1, phi2 = (<PP> self.pairproduction).getPhis(eps, self.k)

                if self.state.genPTR.get_next_float() <= phi1/phiHalf_1:
                    break

            else:
                eps = 1/self.k + (.5 - 1/self.k)*2*self.state.genPTR.get_next_float()
                phi1, phi2 = (<PP> self.pairproduction).getPhis(eps, self.k)

                if self.state.genPTR.get_next_float() <= phi2/phiHalf_2:
                    break
        #else:
            # import time
            # print(">>>> rej samp")
            # print("self.k = ", self.k)
            # print(phiHalf_1, phiHalf_2, u1, u2, u1 + u2)
            
            # time.sleep(1000)
#('self.k = ', 6.6522136006430825)
#(0.02498406417606347, -0.4748447902853812, 0.0, -0.4748447902853812, -0.4748447902853812)
        #SAMPLE THEIR DIRECTION
        # azimuth of both is unif distributed and independent
        
        cdef double E = self.k*Eel0_eV
        
        cdef double Eminus = eps*E - Eel0_eV
        cdef double Eplus  =  E - Eminus - Eel0_eV - Eel0_eV
        
        
        cdef double beta_p = sqrt(Eplus * (Eplus + 2*Eel0_eV))   / (Eplus + Eel0_eV)
        cdef double beta_m = sqrt(Eminus * (Eminus + 2*Eel0_eV)) / (Eminus + Eel0_eV)
        
        cdef double r = 2*self.state.genPTR.get_next_float() - 1
        cdef double cos_p = (r + beta_p)/(r*beta_p + 1)

        r = 2*self.state.genPTR.get_next_float() - 1
        cdef double cos_m = (r + beta_m)/(r*beta_m + 1)
        
        #cdef double phi_p = 2*pi*self.genPTR.get_next_float()
       # cdef double phi_m = 2*pi*self.genPTR.get_next_float()
        
        

        
        self.throwAZIMUTH()
        cdef Electron p
        #Eminus *= 1e6
        if Eminus > CUTOFFel:
            IF _DEBUG: print(" ._pairproduction))")

        
        
            p = Electron._new(self.state)
            p.state.E = Eminus

            
            #p.throwAZIMUTH()
            p.rotateTHETA(cos_m)
            self.nSECONDARY += 1
            self.secondary.append(p)
        else:
            (<V> self.state.current_region).depositLOCAL(self.state.pos, Eminus)
        
        cdef Positron pp
        if Eplus > CUTOFFel:
            pp = Positron._new(self.state)
            
            pp.state.E = Eplus
            pp.state.axis.x *= -1
            pp.state.axis.y *= -1
            pp.state.axis.z *= -1


            #p.throwAZIMUTH()
            pp.rotateTHETA(cos_p)
            self.nSECONDARY += 1
            self.secondary.append(pp)
        else:
            (<V> self.state.current_region).depositLOCAL(self.state.pos, Eplus)

            
        
        IF _DEBUG: print(" ._pairproduction))")
        return

        
        


    cdef void _photoelectric(Photon self):
        """Simulate the photoelectric interaction.
        
        The final result of this interaction is the absorption of the photon. That is,
        the photon simulation ends here.
        
        The algorithm is conceptually very simple:
          (1) Choose the shell containing the electron that will absorb the photon;
          (2) Store the electrons partial state (its energy: E_photon - binding_energy(shell) ); 
          (3) Introduce vacancy in that shell. (instruct pyRelax)
          (4) Run the relaxation model.
          (5) Collect all resulting particles partial states.
          (6) Emmit every particle (resulting from relaxation effects) in a random direction.
          
         All this, while accounting for the energy that is lost during the process. 
         This energy is deposited locally.
         
         NOTE: The electron that left a vacancy is not emited isotropically.
        """
      
      
      
        IF not _PH: return
        IF _DEBUG: print("(( .photoelectric ")
        #self.N_photo += 1
        
        self.state.E = self.k * Eel0_eV
        #cdef int i = self.find_index()
        cdef PARTICLES particles
        #self.curent_material.molecule.PHELionize(self.k*Eel0_eV, &particles)
        (<Mol> self.current_molecule).PHELionize(self.find_index(), self.state.E,  self.state.genPTR, &particles)
        
        
        # first energy value of electrons is the first ejected electron
        cdef double E
        cdef double Etot = self.state.E
        
        
        cdef Photon ph
        for i in range(particles.PHOTONS.size()):
            
            E = particles.PHOTONS.back()
            particles.PHOTONS.pop_back()
            if E < photonCUTOFF: continue
            Etot -= E
            
            ph = Photon._newISOTROPIC(self.state)
            ph.state.E = E
            self.nSECONDARY += 1
            self.secondary.append(ph)
        
        

            
        E = particles.ELECTRONS.back()
        
        
        
        
        if E < CUTOFFel: 
            # E_el = binding_energy, if this guy is bellow cutoff, all electrons produced by the relaxation will be bellow cut off
            (<V> self.state.current_region).depositLOCAL(self.state.pos, Etot)

            return
        particles.ELECTRONS.pop_back()
        
        # simulate ejected electron
        
        
        
        
        # create electrons from relaxation
        cdef Electron el
        
        
        cdef double v
        cdef double A, A2
        cdef double gamma = 1 + E/Eel0_eV
        A = 1/gamma - 1 
        A2 = A +2
        
        cdef double beta = sqrt(E*(E + Eel0_eV))/(E + Eel0_eV)
        
        
        
        cdef double C = .5*beta*gamma*(gamma-1)*(gamma-2)
        
        
        cdef double g0 = 2*(1/A   + C)
        
        cdef double r
        
        
        #(-0.17284654490603235, 1.8271534550939676, 1.2089655115886524, -11.65404076144662, -0.0415417034165238)
        #print(1)
        
        #print(A, A2, gamma, g0, C)
        
        while True:
            
            r = self.state.genPTR.get_next_float()
            v = 2*A/(A2**2 - r) * (2*r + A2*sqrt(r))

            if self.state.genPTR.get_next_float()*g0 < (2 - v)*(1/(A + v)  + C):
                
                
                el = Electron._new(self.state)

                el.state.E = E
                Etot -= E
                
                
                el.throwAZIMUTH()
                
                el.rotateTHETA(1-v)
                break
        
 
        self.secondary.append(el)
        self.nSECONDARY += 1


        for i in range(particles.ELECTRONS.size()):
            E = particles.ELECTRONS.back()
            particles.ELECTRONS.pop_back()
            if E < CUTOFFel: continue
            el = Electron._newISOTROPIC(self.state)
            Etot -= E
            el.state.E = E
            self.nSECONDARY += 1
            self.secondary.append(el)
        
        (<V> self.state.current_region).depositLOCAL(self.state.pos, Etot)
        



    cdef void _tripletproduction(Photon self):
        IF not _TP: return
        #self.N_trip += 1
        
        cdef PARTICLES particles
        (<Mol> self.current_molecule).ionize(self.state.genPTR, &particles)
        cdef Photon ph
        for i in range(particles.PHOTONS.size()):
            
            E = particles.PHOTONS.back()
            particles.PHOTONS.pop_back()
            if E < photonCUTOFF:
                (<V> self.state.current_region).depositLOCAL(self.state.pos, E)
                continue
            
            
            ph = Photon._newISOTROPIC(self.state)
            ph.state.E = E
            self.nSECONDARY += 1
            self.secondary.append(ph)
        
        
        
        cdef Electron el
        for i in range(particles.ELECTRONS.size()):
            E = particles.ELECTRONS.back()
            particles.ELECTRONS.pop_back()
            if E < electronCUTOFF:
                (<V> self.state.current_region).depositLOCAL(self.state.pos, E)
                continue

            el = Electron._newISOTROPIC(self.state)
            el.state.E = E
            self.nSECONDARY += 1
            self.secondary.append(el)

        
        
        self._pairproduction()
        
        
        
        
        

            
        
        
        







    cdef inline int find_index(self):
        cdef int i;
        frexp(self.state.E, &i);

        
        
        #cdef int i = get_exp(E)
        
        # if LIMS[i, 2] is 0:
        #     raise RuntimeError("OUT OF BOUNDS")
        cdef int k = LIMS[i, 2]
        if k is 1:
            return LIMS[i, 0]
        
        if k == 2:
            i = LIMS[i, 0]
            if self.state.E <= EAX[i + 1]: return i
            return i + 1
        
        if k is 3:
            i = LIMS[i, 0]
            if self.state.E <= EAX[i + 1]: return i
            if self.state.E <= EAX[i + 2]: return i + 1
            return i + 2
        
        if k is 4:
            i = LIMS[i, 0]
            if self.state.E <= EAX[i + 1]: return i
            if self.state.E <= EAX[i + 2]: return i + 1
            if self.state.E <= EAX[i + 2]: return i + 2
            return i + 3
        
        cdef int START, END, MID
        START = LIMS[i, 0]
        END   = LIMS[i, 1] 
        
        cdef double xMID
        
        #do binary search 
        while START <= END:
            #find middle
            MID = START + (END - START)//2 #prevents overflow somehow 
            
            xMID = EAX[MID]
            
            if self.state.E is xMID: #found the value
                return MID
            
            if self.state.E < xMID: # discard right side
                END = MID - 1 # do not include mid
                continue
            
            START = MID + 1
        return END 




cdef mixmax_engine GEN # space to store a generator for the python_hooks.Photon

class python_hooks:
    class Photon(Photon):
        """
        # Counters:
        - [ ] cdef int N_coh # Coherent
        - [ ] cdef int N_incoh # Incoherent
        - [ ] cdef int N_photo # Photoelectric
        - [ ] cdef int N_pair # Pair Production 
        - [ ] cdef int N_trip # Triplet Production

        """

        def __init__(self, pos   = np.array([0, 0, 0],  dtype = float),
                           dire  = np.array([0, 0, 1],  dtype = float),
                           axis  = np.array([0, 1, 0],  dtype = float), 
                           double E = 1e6
                    ):

            self.state.pos.x = pos[0]
            self.state.pos.y = pos[1]
            self.state.pos.z = pos[2]

            self.state.dire.x = dire[0]
            self.state.dire.y = dire[1]
            self.state.dire.z = dire[2]

            self.state.axis.x = axis[0]
            self.state.axis.y = axis[1]
            self.state.axis.z = axis[2]

            self.state.E = E





        @staticmethod
        def set_seed(long int seed):
            """Creates a mixmax_engine instance using `seed` and stores it in the module level
            variable `GEN`.
            """
            global GEN
            GEN = mixmax_engine(0, 0, 0, seed)

        def _run(self):
            """Very thin wrapper for the `_run` method.

            Note: The interface between `_run` of this python hook and the actual
            extension type are different:

            ```Cython
            cdef void _run(mixmax_engine *genPTR)
            ```

            versus

            ```
            @staticmethod
            def _run()
            ```

            This is so that the pRNG can be set seperatly in another static method (see `set_seed()`).
            Thus making this wrapper as thin as possible so that it can be properly benchmarked when
            needed.
            """
            global GEN
            (<Photon> self)._run(&GEN)



        def __getattr__(self, attribute):
            if attribute == "E":                return          (<Photon> self).state.E
            if attribute == "IMFP_CUMUL":       return          (<Photon> self).IMFP_CUMUL
            if attribute == "current_material": return  <MAT> ( (<Photon> self).current_material )
            if attribute == "coherent":         return  <COH> ( (<Photon> self).coherent )
            if attribute == "incoherent":       return  <INC> ( (<Photon> self).incoherent )
            if attribute == "pairproduction":   return  <PP>  ( (<Photon> self).pairproduction )
            if attribute == "S":                return  <PPP> ( (<Photon> self).S )
            if attribute == "current_molecule": return  <MOL> ( (<Photon> self).current_molecule )
            if attribute == "k":                return (<Photon> self).k
            if attribute == "pos":
                pos = np.array([0, 0, 0], dtype = float)
                pos[0] = self.state.pos.x
                pos[1] = self.state.pos.y
                pos[2] = self.state.pos.z
            if attribute == "dire":
                dire = np.array([0, 0, 0], dtype = float)
                dire[0] = self.state.dire.x
                dire[1] = self.state.dire.y
                dire[2] = self.state.dire.z
            if attribute == "axis":
                axis = np.array([0, 0, 0], dtype = float)
                axis[0] = self.state.axis.x
                axis[1] = self.state.axis.y
                axis[2] = self.state.axis.z


            if attribute in self.__dict__:
                return self.__dict__[attribute]

            raise AttributeError(f"No attribute named {attribute}")

        def __setattr__(self, attribute, value):
            if   attribute == 'E':                 (<Photon> self).state.E = value
            elif attribute == "current_material":  (<Photon> self).current_material = <void*> value
            elif attribute == "current_region":    (<Photon> self).state.current_region   = <void*> value
            elif attribute == "coherent":          (<Photon> self).coherent = <void*> value
            elif attribute == "incoherent":        (<Photon> self).incoherent = <void*> value
            elif attribute == "pairproduction":    (<Photon> self).pairproduction = <void*> value
            elif attribute == "S":                 (<Photon> self).S = value
            elif attribute == "current_molecule":  (<Photon> self).current_molecule = <void*> value
            elif attribute == "k":                 (<Photon> self).k = value
            elif attribute == "pos":
                self.state.pos.x = value[0]
                self.state.pos.y = value[1]
                self.state.pos.z = value[2]
            elif attribute == "dire":
                self.state.dire.x = value[0]
                self.state.dire.y = value[1]
                self.state.dire.z = value[2]
            elif attribute == "axis":
                self.state.axis.x = value[0]
                self.state.axis.y = value[1]
                self.state.axis.z = value[2]

            else: self.__dict__[attribute] = value

        def _coherent(self):          (<Photon> self)._coherent()
        def _incoherent(self):        (<Photon> self)._incoherent()
        def _pairproduction(self):    (<Photon> self)._pairproduction()
        def _tripletproduction(self): (<Photon> self)._tripletproduction()
        def _incoherent(self):        (<Photon> self)._incoherent()
        def update_references(self):  (<Photon> self).update_references()
        def update_imfp(self):        (<Photon> self).update_imfp()
        def record(self):             (<Photon> self).record()

        def find_index(self): return (<Photon> self).find_index()

        def __repr__(self):
            return "<python_hook.Photon>"

        def __str__(self):
            return "RETURN DEBUG INFO"







# REMOVE WHEN python_hooks IS DONE:

#cdef mixmax_engine gen = mixmax_engine(0,0,0,123);
#
#cdef class PhotonProbe(Photon):
#    
#    
#    def __init__(self, *args, **kwargs):
#        raise RuntimeError("use .init")
#    
#    @staticmethod
#    def init(double E, object current_region):
#        # self = PhotonProbe._new(E,
#        #                         0, 0, 0,
#        #                         0, 1, 0,
#        #                         0, 0, 1,
#        #                         current_region)
#        
#        self = <PhotonProbe>PhotonProbe.__new__(PhotonProbe)
#        
#        self.x = 0
#        self.y = 0
#        self.z = 0
#        self.k = E/Eel0_eV
#        self.state.Eyx = 0
#        self.state.Eyy = 1
#        self.state.Eyz = 0
#        self.state.Ezx = 0
#        self.state.Ezy = 0
#        self.state.Ezz = 1
#        self.state.E = E
#        self.current_region = current_region
#        
#        self.setGEN()
#        self.update_references()
#        self.nSECONDARY = 0
#        self.secondary = deque()
#        return self
#
#    cdef void setGEN(self):
#        self.state.genPTR = &gen
#    
#    
#    def reset(self):
#        self.state.Eyx = 0
#        self.state.Eyy = 1
#        self.state.Eyz = 0
#        
#        self.state.Ezx = 0
#        self.state.Ezy = 0
#        self.state.Ezz = 1
#
#    def test_compton(self, N):
#        k = self.k
#        THETAS = []
#        ENERGY = []
#        for _ in range(N):
#            self._incoherent()
#            
#            ENERGY.append(self.k*Eel0_eV)
#            self.k = k
#            
#            theta = self.state.Ezz
#            
#            THETAS.append(theta)
#            
#            self.reset()
#        import numpy as np
#        THETAS = np.array(THETAS)
#        ENERGY = np.array(ENERGY)
#        THETAS = np.arccos(THETAS)
#        THETAS = THETAS*180/pi
#        
#        return THETAS, ENERGY
#
#
#    def test_coherent(self, N):
#        THETAS = []
#        for _ in range(N):
#            self._coherent()
#            
#            
#            theta = self.state.Ezz
#            
#            THETAS.append(theta)
#            
#            self.reset()
#        import numpy as np
#        THETAS = np.array(THETAS)
#        THETAS = np.arccos(THETAS)
#        THETAS = THETAS*180/pi
#        
#        return THETAS#
