# cython: annotate=True
# distutils: language = c++
print(">>>>>   IMPORTING PHOTONS")


DEF _DEBUG_BASIC = False
DEF _DEBUG = False
DEF _DEBUGincoh = False
DEF RECORD = True

DEF _COH = True
DEF _INCOH = True
DEF _PP = True
DEF _TP = True
DEF _PH = True

from ..random.mixmax.interface cimport mixmax_engine



cdef extern from "<math.h>" nogil:
    double frexp(double x, int* exponent)


#from ..materials.electron.main import eax as _eax
from .._init import eax
from .._init cimport EAX
from .._init cimport LIMS




# cdef double eax[2695] 
# eax[:] = _eax

# from ..materials cimport LIMS 

from libcpp.vector cimport vector

from ..materials.cppRelaxAPI cimport PARTICLES


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
                                             







from collections import deque


#Error messages (to be moved to its own module)
errorMSG1 = "Exhausted allowed number of iterations for rejection sampling."

#External Imports
#from numpy import *
#from numpy.random import rand, randint
#import pickle -> probly not needed any more?





#Local Imports
from .particle cimport Particle
from .particle cimport STATE


from .electrons cimport Electron
from .positrons cimport Positron

from ..tools.vectors cimport Vector
from ..materials import database as db
# --  -- from . import electrons as e
from libc.math cimport sin, cos, log, sqrt, pi , acos, exp
from ..geometry.main cimport Volume
from ..materials.pyRelax cimport Atom as RAtom

#settings
from ..settings import __photonCUTOFF__, DEBUG, __electronCUTOFF__

#from ..materials.photon.CrossSection cimport IMFP






# CONSTANTS AND GLOBALS


cdef double Eel0_MeV = 0.51099895000
cdef double Eel0_eV = Eel0_MeV*1e6





cdef double k_cutoff = __photonCUTOFF__/Eel0_eV


#cdef double CUTOFF = __photonCUTOFF__
cdef double CUTOFFel = __electronCUTOFF__
cimport cython

cdef double photonCUTOFF = __photonCUTOFF__
cdef double electronCUTOFF = __electronCUTOFF__

cdef double minCUTOFF = min(photonCUTOFF, electronCUTOFF)


#from numpy.random import rand 


############################################################################
#       TOOLS

# cdef double rand():
#     cdef double r = randint()
#     return r/RAND_MAX







############################################################################
#       STRUCTURES



#cdef INCOHERENT INCOH




IMFP_CUMUL.C0 = 0.
############################################################################




#cimport cython

#@cython.boundscheck(False)  # Deactivate bounds checking
#@cython.wraparound(False)   # Deactivate negative indexing.


from numpy.random import rand

cdef object choose(list cumul, list items):

    cdef double r = rand()*cumul[-1]
    cdef int i = 0


    for i in range(len(cumul)):
        if cumul[i] > r:
            return items[i-1]


cdef int chooseI(list cumul, int N):
    #cdef int N = len(cumul)
    
    cdef double r = rand()*cumul[N-1]
    cdef int i = 0

    for i in range(N):
        if cumul[i] > r:
            return i-1


cdef struct INCOHERENT:
    double t1
    double t2
    double tau_min
    double tau, cos, N, D, sin2, x , T
    double k

AA = 0

REJECTED = 0
REJECTED2 = 0



@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef class Photon(Particle):

    cdef double ENERGY(self):
        return self.k*Eel0_eV
    
    
  #  cdef Photon _new(STATE state):
  #      cdef Photon self
 #       self = <Photon>Photon.__new__(Photon)
#
#
 #       self.state = state


    @staticmethod
    cdef Photon _new(double E,
                          double x, double y, double z,
                          double eyx, double eyy, double eyz,
                          double ezx, double ezy, double ezz, 
                          Volume current_region):
        
        cdef Photon self
        self = <Photon>Photon.__new__(Photon)

        self.current_region = current_region
        self.x = x# +1
        self.y = y
        self.z = z
        
        #print(self.x, self.current_region.p.x)
        

        
        self.k = E/Eel0_eV
        self.eyx = eyx
        self.eyy = eyy
        self.eyz = eyz
        self.ezx = ezx
        self.ezy = ezy
        self.ezz = ezz
        self.E = E
        
        self.index = <int> (10*self.z)
        
        
        return self
    
    @staticmethod
    cdef Photon _newISOTROPIC(double E, double x, double y, double z, Volume current_region, mixmax_engine* genPTR):
        
        cdef Photon self
        self = <Photon>Photon.__new__(Photon)
        self.current_region = current_region
        

        self.x = x
        self.y = y
        self.z = z
        self.E = E
        self.k = E/Eel0_eV
        
        cdef double a

        while True:
            x = 2*genPTR.get_next_float() - 1
            y = 2*genPTR.get_next_float() - 1
            a = x**2 + y**2
            if a < 1:
                break

        self.ezx = 1 - 2*a 
        a = 2 * sqrt(1 - a)
        self.ezy = x*a
        self.ezz = y*a
        #self.throwAZIMUTH()
        
        
     #   self.index = <int> (10*self.z)

        
        return self    





    # @staticmethod #custom constructor for speed
    # cdef Photon _new(Volume current_region,
    #                  double E, 
    #                  Vector pos,
    #                  Vector ey,
    #                  Vector ez):

    #     cdef Photon self
    #     self = <Photon>Photon.__new__(Photon)
  
        

        
    #     Particle.INIT(self, E, pos,
    #                      ey, ez,
    #                      current_region)
        
    #     #if DEBUG:
    #     #    self.FILE.write("""
    #     #                NEW PHOTON {self}
    #     #                    """)
        
    #     return self

    # @staticmethod #custom constructor for speed
    # cdef Photon _newISOTROPIC(Volume current_region,
    #                           double E, 
    #                           Vector pos):    
    #     cdef Photon self
    #     self = <Photon>Photon.__new__(Photon)
  
    #     self.k = E/E0_el

        
    #     axis = ez0
    #     ey = ey0.rotateAngle(axis, 2*pi*rand())
            
    #     axis = ey
    #     ez = ez0.rotateAngle(axis, pi*rand())
        
        
        
    #     Particle.INIT(self, E, pos,
    #                      ey, ez,
    #                      current_region)
        
    #     #if DEBUG:
    #     #    self.FILE.write("""
    #     #                NEW PHOTON {self}
    #     #                    """)
        
    #     return self        
    
    
    

    ####################################################################################
    ########                           RUN                                      ########
    ########                           RUN                                      ########
    ####################################################################################

    
    cdef void _run(Photon self, mixmax_engine* genPTR):
        IF _DEBUG_BASIC: print("> PHOTON")
        
        #cdef double r
        self.secondary = deque()
        self.nSECONDARY = 0
        if self.k < k_cutoff:
            return
        
        self.genPTR = genPTR
        
        self.update_references()
        
        IF RECORD: self.record()
        
        cdef double r
        cdef double L 
        #cdef int _
        IF _DEBUG: print("STARTED")
        

        while True:
            self.state.L = -log(1e-5 + (1-1e-5)*self.state.genPTR.get_next_float())/self.imfp_T

            if self.state.L > 1e3: print("PHOTONS")

            if (<V> self.state.current_region).move(self.state):
                if (<V> self.state.current_region).vaccum:
                    (<V> self.state.current_region).exit()
                    return

                self.update_references()

                IF RECORD: self.record()

                continue
            
            
            IF RECORD: self.record()
            
            r = self.state.genPTR.get_next_float()*self.IMFP_CUMUL.C5
            
            if   r < self.IMFP_CUMUL.C1:
                self._coherent()
                
                
            elif r < self.IMFP_CUMUL.C2: 
                self._incoherent()

                if self.k < k_cutoff:
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



   
 



    # cdef int propagate(self) except -1:
    #     cdef long double r = self.genPTR.get_next_float()
        
    #     if self.imfp_T == 0:
    #         raise StopSimulation
        
    #     cdef long double L = -log(r)/self.imfp_T
        
    #     if L > 1e10:
    #         raise StopSimulation
        
    #     if self.cross(L):
    #         self.update_references()
            
    #         if self.current_region.imp is 0:
    #             raise StopSimulation

    #         return self.propagate()
    
    #     self.move(L)
    #     return 0


    cdef void record(self):
        """
        Add current position to the track record.
        """
        
        self.X.push_back(self.x)
        self.Y.push_back(self.y)
        self.Z.push_back(self.z)
        
        #to_log = f"""
#{self.pos.x} |  {self.pos.y} | {self.pos.z}  |  {self.ez} | {self.current_region}  |  {self.current_material.molecule.formula}"""
        
        #if DEBUG:
        #    self.FILE.write(to_log)





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


        self.current_material =  <void*> (<V> handler).material
        self.current_molecule =  <void*> (<M> (<V> handler).material).molecule
        
        
        cdef void* photon = <void*> (<V> self.current_material).photon

        #these references are used by their corresponding _fooInteraction method
        self.coherent          = (<Ph> photon).coherent
        self.incoherent        = (<Ph> photon).incoherent
        #self.photoelectric     = self.current_material.photoelectric
        self.pairproduction    = (<Ph> photon).pairproduction
        self.tripletproduction = (<Ph> photon).tripletproduction

        #self.S = self.incoherent.S

        #since region crossing has ocurred, update the inverse mean free paths
        self.update_imfp()


    cdef void update_imfp(Photon self):
        """
        Updates inverse mean free paths. Constructs cumul.
        Called when there is a region crossing or energy of photon has changed.
        """

        self.STATE.E = self.k * Eel0_eV
        cdef int i = self.find_index()
        #self.energy.push_back(self.E)
        
        #IMFP_CUMUL.C0 = 0.
        self.IMFP_CUMUL.C1 = (<Coh> self.coherent).imfpA[i] + self.STATE.E*(<Coh> self.coherent).imfpB[i]
        
        self.IMFP_CUMUL.C2 = self.IMFP_CUMUL.C1 + (<inCoh> self.incoherent).imfpA[i]      + self.state.E*(<inCoh> self.incoherent).imfpB[i]
        self.IMFP_CUMUL.C3 = self.IMFP_CUMUL.C2 + (<PP> self.pairproduction).imfpA[i]     + self.state.E*(<PP> self.pairproduction).imfpB[i]
        self.IMFP_CUMUL.C4 = self.IMFP_CUMUL.C3 + (<PPP> self.tripletproduction).imfpA[i] + self.state.E*(<PPP> self.tripletproduction).imfpB[i]
        self.IMFP_CUMUL.C5 = self.IMFP_CUMUL.C4 + (<Mol> (<Coh> self.coherent)).PHELa[i]  + self.state.E*(<Mol> self.current_molecule).PHELb[i]

        self.imfp_T = self.IMFP_CUMUL.C5


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
            #r = self.genPTR.get_next_float()*cumulMAX
            r = self.state.genPTR.get_next_float()*cumulMAX
            #print(r)
            x2 = (<Coh> self.coherent).evalX(r)
            cos = 1 - x2/k2 #x2/k2 = 1 - cos
            if self.genPTR.get_next_float()*2 < 1 + cos*cos:
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

        
        IF _DEBUG: print("(( ._incoherent")
        IF _DEBUG: print(f"STARTING: N = {self.N_incoh} | k =  {self.k} | E = {self.k*Eel0_eV}")
        
        
        # (1) Choose an atom with probabilities based on Z
     #   print("-------START COMPTON")

        cdef Atom active_atom = (<Mol> self.current_molecule).choose_atom(self.genPTR)
        cdef double Uk

        # (2) Sample cos(theta) from KN DCS
        
        cdef double eps_min = 1/(1 + 2*self.k) 
        cdef double alpha1 = log(1 + 2*self.k)
        cdef double alpha2 = (1 - eps_min*eps_min)/2
        cdef double eps, oneMINUScos
        cdef double g, cos, sin2
        cdef double gmax = eps_min + 1/eps_min
        
        while 1:
            if (alpha1 + alpha2)*self.genPTR.get_next_float() <= alpha1:
                eps = eps_min*exp(alpha1 * self.genPTR.get_next_float())
            else:
                eps = sqrt(eps_min*eps_min + 2*alpha2*self.genPTR.get_next_float() )
    
            oneMINUScos = (1/eps - 1)/self.k
            sin2 = 1 - (1 - oneMINUScos)**2
            
            #g = (1/eps + eps - sin2)/gmax
            g = 1 - eps*sin2/(1 + eps*eps)
            if self.genPTR.get_next_float() < g*active_atom.S._eval(oneMINUScos*(self.k)**2)  :
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
                R = self.genPTR.get_next_float()*active_atom.Nsh
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
            pz = active_shell.sample_compton_profile(self.genPTR, 
                                                     (alpha1 - Uk)/sqrt(2*alpha1 - Uk*Uk) #= pz_max
                                                     )
            
            
            proposed_k = (1 - pz*pz*eps*cos + pz*sqrt( 1 - 2*eps*cos + eps*eps*(1 - pz*pz*sin2))  )/(1 - pz*pz*eps*eps)

            if 0 < proposed_k*eps < self.genPTR.get_next_float():
                new_k = proposed_k*eps*self.k
                break
    
        else: new_k = eps*self.k

        
        Eel = (self.k - new_k - Uk)*Eel0_eV
        
        self.throwAZIMUTH()
        
        
        if Eel > electronCUTOFF:
            el = Electron._new(Eel,
                         self.x, self.y, self.z,
                         -self.eyx, -self.eyy, -self.eyz,
                         self.ezx, self.ezy, self.ezz,
                         self.current_region)
            
            el.rotateTHETA((self.k - new_k*cos)/sqrt(new_k*new_k + self.k*self.k - 2*new_k*self.k*cos))
            self.nSECONDARY += 1
            self.secondary.append(el)
            
        self.k = new_k

        self.rotateTHETA(cos)
        
        
        if active_shell.binding_energy < minCUTOFF:
            return
        
        ## relaxation
        cdef PARTICLES particles
        active_atom.ionize(active_shell.index, self.genPTR, &particles)
        
        
        cdef double E

        
        
        cdef Photon ph
        for i in range(particles.PHOTONS.size()):
            
            E = particles.PHOTONS.back()
            particles.PHOTONS.pop_back()
            if E < photonCUTOFF: continue
            
            
            ph = Photon._newISOTROPIC(E, self.x, self.y, self.z, self.current_region,  self.genPTR)
            self.nSECONDARY += 1
            self.secondary.append(ph)
        

        for i in range(particles.ELECTRONS.size()):
            E = particles.ELECTRONS.back()
            particles.ELECTRONS.pop_back()
            if E < CUTOFFel: continue
            el = Electron._newISOTROPIC(E, self.x, self.y, self.z, self.current_region, self.genPTR)
            self.nSECONDARY += 1
            self.secondary.append(el)
        
        
        
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
        #         #self.E = E0_el*self.k
                


        #         self.E = self.k*Eel0_eV
                
        #         self.k = INCOH.tau*INCOH.k
                
                
                
                
        #         #E_new = E0_el*self.k
        #         E_new = self.k*Eel0_eV
                
        #         cos = (self.E - E_new*INCOH.cos)/sqrt(self.E**2 + E_new**2 - 2*self.E*E_new*INCOH.cos)
                
                
        #         IF _DEBUGincoh: print(f"CHANGING DIRECTION: cos = {INCOH.cos} ")
        #         self.throwAZIMUTH()
                
            
        #         IF _DEBUGincoh: print(f"passed azimuth 1")
                
                
                
                
                
                
                
        #         E_new = (self.E - E_new)
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
        #         if E_new < CUTOFFel:
        #             self.rotateTHETA(INCOH.cos)
        #             return
                
                
        #         IF _DEBUGincoh: print("--- creating secondary ")
                
        #         p = Electron._new(E_new, self.x, self.y, self.z,
        #                                       self.eyx, self.eyy,self.eyz,
        #                                       self.ezx, self.ezy,self.ezz, 
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

            
            if True if (u1 + u2)*self.genPTR.get_next_float() < u1 else False:
                eps  = .5 + (.5 - 1/self.k)*(2*self.genPTR.get_next_float() - 1)**(1/3)
                phi1, phi2 = (<PP> self.pairproduction).getPhis(eps, self.k)

                if self.genPTR.get_next_float() <= phi1/phiHalf_1:
                    break

            else:
                eps = 1/self.k + (.5 - 1/self.k)*2*self.genPTR.get_next_float()
                phi1, phi2 = (<PP> self.pairproduction).getPhis(eps, self.k)

                if self.genPTR.get_next_float() <= phi2/phiHalf_2:
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
        cdef double Eplus  =  E - Eminus - 2*Eel0_eV
        
        
        cdef double beta_p = sqrt(Eplus * (Eplus + 2*Eel0_eV))   / (Eplus + Eel0_eV)
        cdef double beta_m = sqrt(Eminus * (Eminus + 2*Eel0_eV)) / (Eminus + Eel0_eV)
        
        cdef double r = 2*self.genPTR.get_next_float() - 1
        cdef double cos_p = (r + beta_p)/(r*beta_p + 1)

        r = 2*self.genPTR.get_next_float() - 1
        cdef double cos_m = (r + beta_m)/(r*beta_m + 1)
        
        #cdef double phi_p = 2*pi*self.genPTR.get_next_float()
       # cdef double phi_m = 2*pi*self.genPTR.get_next_float()
        
        
        
        
        self.throwAZIMUTH()
        cdef Electron p
        #Eminus *= 1e6
        if Eminus > CUTOFFel:
            IF _DEBUG: print(" ._pairproduction))")

        
        
            p = Electron._new(Eminus, self.x, self.y, self.z,
                                          self.eyx, self.eyy,self.eyz,
                                          self.ezx, self.ezy,self.ezz, 
                                          self.current_region)
            
            #p.throwAZIMUTH()
            p.rotateTHETA(cos_m)
            self.nSECONDARY += 1
            self.secondary.append(p)
        
        cdef Positron pp
        if Eplus > CUTOFFel:
            pp = Positron._new(Eplus, self.x, self.y, self.z,
                                          -self.eyx, -self.eyy,-self.eyz,
                                          self.ezx, self.ezy,self.ezz, 
                                          self.current_region)
            
            #p.throwAZIMUTH()
            pp.rotateTHETA(cos_p)
            self.nSECONDARY += 1
            self.secondary.append(pp)
            
            
        
        IF _DEBUG: print(" ._pairproduction))")
        return

        
        


    cdef void _photoelectric(Photon self):
        IF not _PH: return
        IF _DEBUG: print("(( .photoelectric ")
        #self.N_photo += 1
        
        self.E = self.k * Eel0_eV
        #cdef int i = self.find_index()
        cdef PARTICLES particles
        #self.curent_material.molecule.PHELionize(self.k*Eel0_eV, &particles)
        (<Mol> self.current_molecule).PHELionize(self.find_index(), self.E,  self.genPTR, &particles)
        
        
        # first energy value of electrons is the first ejected electron
        cdef double E

        
        
        cdef Photon ph
        for i in range(particles.PHOTONS.size()):
            
            E = particles.PHOTONS.back()
            particles.PHOTONS.pop_back()
            if E < photonCUTOFF: continue
            
            
            ph = Photon._newISOTROPIC(E, self.x, self.y, self.z, self.current_region,  self.genPTR)
            self.nSECONDARY += 1
            self.secondary.append(ph)
        
        

            
        E = particles.ELECTRONS.back()
        
        
        
        
        if E < CUTOFFel: 
            # E_el = binding_energy, if this guy is bellow cutoff, all electrons produced by the relaxation will be bellow cut off
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
            
            r = self.genPTR.get_next_float()
            v = 2*A/(A2**2 - r) * (2*r + A2*sqrt(r))

            if self.genPTR.get_next_float()*g0 < (2 - v)*(1/(A + v)  + C):
                
                
                el = Electron._new(E, self.x, self.y, self.z,
                                   self.eyx, self.eyy, self.eyz,
                                   self.ezx, self.ezy, self.ezz,
                                   self.current_region)
                
                
                el.genPTR = self.genPTR
                el.throwAZIMUTH()
                
                el.rotateTHETA(1-v)
                break
 
        
            
        self.secondary.append(el)
        self.nSECONDARY += 1
        
        
        for i in range(particles.ELECTRONS.size()):
            E = particles.ELECTRONS.back()
            particles.ELECTRONS.pop_back()
            if E < CUTOFFel: continue
            el = Electron._newISOTROPIC(E, self.x, self.y, self.z, self.current_region, self.genPTR)
            self.nSECONDARY += 1
            self.secondary.append(el)
        
        
        
        
        

    cdef void _tripletproduction(Photon self):
        IF not _TP: return
        #self.N_trip += 1
        
        cdef PARTICLES particles
        (<Mol> self.current_molecule).ionize(self.genPTR, &particles)
        
        cdef Photon ph
        for i in range(particles.PHOTONS.size()):
            
            E = particles.PHOTONS.back()
            particles.PHOTONS.pop_back()
            if E < photonCUTOFF: continue
            
            
            ph = Photon._newISOTROPIC(E, 
                                      self.x, self.y, self.z, 
                                      self.current_region,
                                      self.genPTR)
            self.nSECONDARY += 1
            self.secondary.append(ph)
        
        
        
        cdef Electron el
        for i in range(particles.ELECTRONS.size()):
            E = particles.ELECTRONS.back()
            particles.ELECTRONS.pop_back()
            if E < electronCUTOFF: continue
            el = Electron._newISOTROPIC(E, 
                                        self.x, self.y, self.z, 
                                        self.current_region, 
                                        self.genPTR)
            self.nSECONDARY += 1
            self.secondary.append(el)

            
        
        
        self._pairproduction()
        
        
        
        
        

            
        
        
        







    cdef inline int find_index(self):
        cdef int i;
        frexp(self.E, &i);

        
        
        #cdef int i = get_exp(E)
        
        # if LIMS[i, 2] is 0:
        #     raise RuntimeError("OUT OF BOUNDS")
        cdef int k = LIMS[i, 2]
        if k is 1:
            return LIMS[i, 0]
        
        if k == 2:
            i = LIMS[i, 0]
            if self.E <= EAX[i + 1]: return i
            return i + 1
        
        if k is 3:
            i = LIMS[i, 0]
            if self.E <= EAX[i + 1]: return i
            if self.E <= EAX[i + 2]: return i + 1
            return i + 2
        
        if k is 4:
            i = LIMS[i, 0]
            if self.E <= EAX[i + 1]: return i
            if self.E <= EAX[i + 2]: return i + 1
            if self.E <= EAX[i + 2]: return i + 2
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
            
            if self.E is xMID: #found the value
                return MID
            
            if self.E < xMID: # discard right side
                END = MID - 1 # do not include mid
                continue
            
            START = MID + 1
        return END 










    ####################################################################################
    ########                           PYTHON                                   ########
    ########                          INTERFACE                                 ########
    ####################################################################################
    
    
    
    # @staticmethod #thin wrapper for python acess
    # def new(Volume space, 
    #         Volume current_region,
    #         E     = 6.,
    #         pos   = Vector(0., 0., 0.),
    #         theta = 0.,
    #         phi   = 0.,
    #         ex    = Vector(1., 0., 0.), 
    #         ey    = Vector(0., 1., 0.), 
    #         ez    = Vector(0., 0., 1.),
    #         simulate_secondary = False):
        
        
    #     return Photon._new(space, current_region, E, pos, 
    #                        theta, phi, 
    #                        ex, ey, ez, 
    #                        simulate_secondary)
    
    #thin wrapper for python access

    
    
    
    def __repr__(self):
        return f"<Photon: pos = {self.x},{self.y},{self.z} , ez = {self.ez}, E = {self.k*Eel0_eV} eV>"

    def __str__(self):
        string = f"""
        Photon:
            pos = {self.pos} ;
            ex = {self.ex} 
            ey = {self.ey} 
            ez = {self.ez} ;
            E = {self.k*Eel0_eV}eV ;
            
        Number of Interactions:
            coh: {self.N_coh}
            incoh: {self.N_incoh}
            pair: {self.N_pair}
            trip: {self.N_trip}
            photo: {self.N_photo}
            
        Material (rho = {(<Mat> self.current_material).density}): 
            {(<Mat> self.current_material).molecule.formula}
        
        Current Region:
            {(<V> self.current_region)}
            
        
        
        """
        return string

    def getTrack(self):
        return self.X, self.Y, self.Z
    
    def getZZ(self):
        print(self.ZZ)
        
    def getEnergy(self):
        return self.energy
    
    def printTrack(self):
        print("x", "y", "z")
        for x, y, z, E in zip(self.X, self.Y, self.Z, self.energy):
            print(f"{x}                 {y}                 {z}                 {E}")
    
    def plotTrack(self, **kwargs):
        import mayavi.mlab as mlab
        mlab.plot3d(self.X, self.Y, self.Z, **kwargs)
        mlab.show()







cdef mixmax_engine gen = mixmax_engine(0,0,0,123);

cdef class PhotonProbe(Photon):
    
    
    def __init__(self, *args, **kwargs):
        raise RuntimeError("use .init")
    
    @staticmethod
    def init(double E, object current_region):
        # self = PhotonProbe._new(E,
        #                         0, 0, 0,
        #                         0, 1, 0,
        #                         0, 0, 1,
        #                         current_region)
        
        self = <PhotonProbe>PhotonProbe.__new__(PhotonProbe)
        
        self.x = 0
        self.y = 0
        self.z = 0
        self.k = E/Eel0_eV
        self.eyx = 0
        self.eyy = 1
        self.eyz = 0
        self.ezx = 0
        self.ezy = 0
        self.ezz = 1
        self.E = E
        self.current_region = current_region
        
        self.setGEN()
        self.update_references()
        self.nSECONDARY = 0
        self.secondary = deque()
        return self

    cdef void setGEN(self):
        self.genPTR = &gen
    
    
    def reset(self):
        self.eyx = 0
        self.eyy = 1
        self.eyz = 0
        
        self.ezx = 0
        self.ezy = 0
        self.ezz = 1

    def test_compton(self, N):
        k = self.k
        THETAS = []
        ENERGY = []
        for _ in range(N):
            self._incoherent()
            
            ENERGY.append(self.k*Eel0_eV)
            self.k = k
            
            theta = self.ezz
            
            THETAS.append(theta)
            
            self.reset()
        import numpy as np
        THETAS = np.array(THETAS)
        ENERGY = np.array(ENERGY)
        THETAS = np.arccos(THETAS)
        THETAS = THETAS*180/pi
        
        return THETAS, ENERGY


    def test_coherent(self, N):
        THETAS = []
        for _ in range(N):
            self._coherent()
            
            
            theta = self.ezz
            
            THETAS.append(theta)
            
            self.reset()
        import numpy as np
        THETAS = np.array(THETAS)
        THETAS = np.arccos(THETAS)
        THETAS = THETAS*180/pi
        
        return THETAS