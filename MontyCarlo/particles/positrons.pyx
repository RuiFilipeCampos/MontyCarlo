# profile = False
# cython: annotate=False
# distutils: language = c++ 

#    ____       U  ___ u   ____                     _____      ____        U  ___ u   _   _     
#  U|  _"\ u     \/"_ \/  / __"| u       ___       |_ " _|  U |  _"\ u      \/"_ \/  | \ |"|    
#  \| |_) |/     | | | | <\___ \/       |_"_|        | |     \| |_) |/      | | | | <|  \| |>   
#   |  __/   .-,_| |_| |  u___) |        | |        /| |\     |  _ <    .-,_| |_| | U| |\  |u   
#   |_|       \_)-\___/   |____/>>     U/| |\u     u |_|U     |_| \_\    \_)-\___/   |_| \_|    
#   ||>>_          \\      )(  (__) .-,_|___|_,-.  _// \\_    //   \\_        \\     ||   \\,-. 
#  (__)__)        (__)    (__)       \_)-' '-(_/  (__) (__)  (__)  (__)      (__)    (_")  (_/  

print("Importing `.particles.positrons`")





#Error messages (to be moved to its own module)
errorMSG1 = "Exhausted allowed number of iterations for rejection sampling."





# Conditional Compilation for Debugging. 
DEF _DEBUG_BASIC = False
DEF _SIGNAL_INTERACTION = False
DEF RECORD = True





# Internal Imports
from ..materials import database as db
from ..settings import __photonCUTOFF__
from ..settings import __electronCUTOFF__

from .particle cimport Particle
from .photons  cimport Photon
from ..geometry.main cimport Volume
from ..tools.vectors cimport Vector
from ..materials.materials cimport Material
from ..materials.materials cimport Material
from ..materials.electron.main cimport Brem
from ..materials.electron.main cimport Inelastic
from ..materials.electron.main cimport Elastic
from ..materials.electron.main cimport DIST
from ..external.mixmax_interface cimport mixmax_engine
from ..materials.cppRelaxAPI cimport PARTICLES

from .._init  import eax
from .._init  cimport EAX
from .._init  cimport LIMS

#External Imports
import numpy as np
from collections import deque

from libcpp.vector cimport vector
from libc.math cimport isnan
from libc.math cimport sin 
from libc.math cimport cos
from libc.math cimport log
from libc.math cimport sqrt
from libc.math cimport pi
from libc.math cimport exp
from libc.math cimport fmin
from libc.math cimport fmax
from libc.math cimport acos
from libc.math cimport pow
from libc.stdlib cimport rand # Deprecated
from libc.stdlib cimport RAND_MAX # Deprecated
from libc.stdlib cimport srand # Deprecated
cimport cython

cdef extern from "<math.h>" nogil:
    double frexp(double x, int* exponent)





# CONSTANTS AND GLOBALS
cdef extern from "math.h":
    float INFINITY

cdef double twoPI = 2*pi
cdef double  ELECTRON_REST_MASS      = 0.51099895000e6
cdef double CUTOFF = __photonCUTOFF__
cdef double photonCUTOFF = __photonCUTOFF__
cdef double MIN_CUT_OFF = min(__photonCUTOFF__, __electronCUTOFF__)
cdef double E0_el = db.E0_electron*1e-3
cdef double CUT_OFF = __electronCUTOFF__
cdef double ELECTRON_REST_ENERGY = 0.51099895000*1e6 #eV
cdef double  _2ELECTRON_REST_ENERGY    = 2 *ELECTRON_REST_ENERGY
cdef double[::1] LOGeax = np.log(eax)
cdef double[::1] diffLOGeax = np.diff(np.array(LOGeax))


cdef inline double g(double v, double gamma):
    return -(gamma+1)**2*v + (gamma**2 + 4*gamma +1) - 1/v





ctypedef Volume V





@cython.boundscheck(False)
@cython.wraparound(False) 
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef class Positron(Particle):
    """Simulates positrons using a class II condensed history scheme.

    Catastrophic Interactions:
        Elastic Scattering
        Bremsstrahlung Emission
        Impact Ionization and/or Excitation.

    """





    @staticmethod
    cdef Positron _new(STATE& state):
        cdef Positron self
        self = <Positron>Positron.__new__(Positron)
        self.state = state
        return self
    
    @staticmethod
    cdef Positron _newISOTROPIC(STATE& state):
        """Create a new `Positron` from the provided `state`. Returns it pointing in a random
        direction but without a polarization vector (no azimuth).
        """

        cdef Positron self
        self = <Positron>Positron.__new__(Positron)
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
        self.index = <int> (10*self.z)
        return self    

    cdef double ENERGY(self):
        return self.state.E





    ####################################################################################
    ########                           RUN                                      ########
    ########                           RUN                                      ########
    ####################################################################################

    cdef void _run(self, mixmax_engine *genPTR):
        """Perform the simulation. Puts all the pieces together.
        """

        IF _DEBUG_BASIC: print("_DEBUG_BASIC: `positron._run(genPTR)`")

        # This needs to be set even if particle is bellow threshold.
        self.secondary  = deque()
        self.nSECONDARY = 0

        if self.state.E < CUT_OFF:
            (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
            return
        
        self.state.genPTR = genPTR
        self.update_references()
        
        IF RECORD: self.record()





        cdef double r           # Random value between 0 and self.imfp_max
        cdef double tau         # Size of first segment.
        cdef double tau2        # Size of second segment.
        cdef double S_soft      # Soft/Restricted stopping power.
        cdef bint delta = False # Indicates if current iteration ends in a delta interaction.

        while True:

            # To be removed:
            if self.state.pos.x**2 + self.state.pos.y**2 + self.state.pos.z**2 > 10_000**2:
                return

            # Choose the total path to be travelled.
            self.s = -log(self.state.genPTR.get_next_float() )/self.imfp_max

            if self.s > self.s_max:
                self.s = self.s_max
                delta = True





            # FIRST PART OF TRAJECTORY
            tau = self.s*self.state.genPTR.get_next_float()
            self.sample_w(self.s)
            S_soft = self.w/self.s
            self.state.L = tau

            if (<V> self.state.current_region).move(self.state, 0):
                #self.state.E -= (tau - self.state.L)*self.w/self.s
                #if self.state.E < CUT_OFF:
                    #(<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
                    #return
                
                self.update_references()
                continue


            self.state.E -= tau*S_soft
            self.do_hinge()

            if self.state.E < CUT_OFF:
                (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
                (<Volume> self.state.current_region).exit()
                return

            (<V> self.state.current_region).depositLOCAL(self.state.pos, self.w)


            IF RECORD: self.record()
            
 
            
            tau = self.s - tau
            self.state.L = tau

            if (<V> self.state.current_region).move(self.state, 0):
                #self.state.E -= (tau - self.state.L)*S_soft
                #if self.state.E < CUT_OFF:
                    #(<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
                    #(<Volume> self.state.current_region).exit()
                    #return

                self.update_references()
                continue
            self.state.E -= S_soft*tau

            IF RECORD: self.record()
            
            if self.state.E < CUT_OFF:
                (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
                (<Volume> self.state.current_region).exit()
                return
            
            
            if delta:
                self.update_imfp()
                delta = False
                continue
            
            
            self.update_imfp_cumul()
            
            r = self.state.genPTR.get_next_float()*self.imfp_max
            
            if r < self.IMFP_CUMUL.C0:
                self._anihilation()

                (<Volume> self.state.current_region).exit()
                return
            
            elif r < self.IMFP_CUMUL.C1: 
                self._inelastic()
                if self.state.E < CUT_OFF:
                    (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)

                    (<Volume> self.state.current_region).exit()
                    return
                self.update_imfp()
                
                
            elif r < self.IMFP_CUMUL.C2:
                self._brem()
                if self.state.E < CUT_OFF:
                    (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
                    (<Volume> self.state.current_region).exit()
                    return
                self.update_imfp()
                
            elif r < self.IMFP_CUMUL.C3: 
                self._elastic()





    ####################################################################################
    ########                           UPDATE                                   ########
    ########                           METHODS                                  ########
    ####################################################################################

    cdef void update_references(self) :
        """
        Updates all references. Called when:
            - there is a region crossing;
            - particle is being run for the first time.
        """
        
        #getting material from current region
        self.current_material = (<Volume> self.state.current_region).material

        self.GOS = self.current_material.positron.inelastic.gosMOLECULE
        
        self.positron = self.current_material.positron

        #these references are used by their corresponding _fooInteraction method
        self.elastic    = self.positron.elastic
        self.inelastic  = self.positron.inelastic
        self.brem       = self.positron.brem
        self.anih       = self.positron.anih

        #since region crossing has ocurred, update the inverse mean free paths
        self.update_imfp()



    cdef inline void update_imfp_cumul(Positron self):
        """
        Updates inverse mean free paths. Constructs cumul.
        Called when there is a region crossing or energy of photon has changed.
        """

        cdef int i = self.find_index()  
 
        self.IMFP_CUMUL.C0 = self.anih.imfpA[i] + self.anih.imfpB[i]*self.state.E
        self.IMFP_CUMUL.C1 = self.IMFP_CUMUL.C0 + self.inelastic.imfpA[i] + self.inelastic.imfpB[i]*self.state.E
        self.IMFP_CUMUL.C2 = self.IMFP_CUMUL.C1 + self.brem.imfpA[i] +    self.brem.imfpB[i]   *self.state.E
        self.IMFP_CUMUL.C3 = self.IMFP_CUMUL.C2 + self.elastic.imfpA[i] + self.elastic.imfpB[i] * self.state.E
        
        #self.positron.IMFP(self.state.E)
        self.s_max = 4/self.IMFP_CUMUL.C3
        cdef double Emin = self.state.E - self.positron.find_wmax(self.s_max, self.state.E)
        
        
        if Emin is 0:
            self.imfp_max = self.IMFP_CUMUL.C3
            return
        
        i = self.positron.find_index(Emin)
        
        
        self.imfp_max = fmax(self.IMFP_CUMUL.C3, 
                             self.positron.imfpA[i] + self.positron.imfpB[i] * Emin
                             )
        
        
        
    cdef inline void _delta(self):
        pass
        #imfp_delta = self.imfp_max - self.IMFP_CUMUL.C3
        
        
        
        self.move(  -log(1e-2 + (1 - 1e-2)*self.state.genPTR.get_next_float()  )/( self.imfp_max - self.IMFP_CUMUL.C3)  )
        # if self.z == INFINITY or self.z == -INFINITY:
        #     print(self.imfp_max, self.IMFP_CUMUL.C3)
        
        
        
    cdef inline void update_imfp(Positron self) :
        
        cdef double Emin = self.state.E - self.positron.find_wmax(self.s_max, self.state.E)
        cdef int i, j
        
        if Emin is 0:
            j = self.find_index()
            self.imfp_max = self.positron.imfpA[j] + self.positron.imfpB[j] * self.state.E
            self.s_max = 4/self.imfp_max
            return
        
        
        
        i = self.positron.find_index(Emin)
        j = self.find_index()
        #cdef double imfpH = 
   
        self.s_max = self.positron.imfpA[j] + self.positron.imfpB[j] * self.state.E
        self.imfp_max = fmax(self.s_max,
                             self.positron.imfpA[i] + self.positron.imfpB[i] * Emin
                             )

        self.s_max = 4/self.s_max





    ####################################################################################
    ########                           RANDOM                                   ########
    ########                           SAMPLING                                  ########
    ####################################################################################
    
    
    
    cdef inline void sample_w(self, double tau):
        
        
        #self.SP     = self.inelastic.softSP._eval(self.state.E) + self.brem.softSP._eval(self.state.E)
        #self.STRAGG = self.inelastic.softSTRAGG._eval(self.state.E) + self.brem.softSTRAGG._eval(self.state.E)

        #self.avgW   = self.s * (self.inelastic.softSP._eval(self.state.E)     + self.brem.softSP._eval(self.state.E)    )
        #self.varW   = self.s * (self.inelastic.softSTRAGG._eval(self.state.E) + self.brem.softSTRAGG._eval(self.state.E))
        
        cdef int i = self.find_index()

        cdef double SP = ( self.positron.softSPA[i]     + self.positron.softSPB[i]*self.state.E)
        cdef double STRAGG = ( self.positron.softSTRAGGA[i] + self.positron.softSTRAGGB[i]*self.state.E)


        self.avgW = self.s * SP * (1 - .5*self.positron.softSPB[i]*self.s)
        self.varW = self.s * STRAGG  - self.s*self.s*(.5*self.positron.softSTRAGGB[i]*SP + STRAGG*self.positron.softSPB[i])

        cdef double sigma = sqrt(self.varW)
        cdef double avgW2 = pow(self.avgW, 2)
        
        if avgW2 > 9*self.varW:

            self.w = 1.015387*sigma*self.current_material.positron.gauss._sample() + self.avgW #CONFIRM
            return
        
        if avgW2 > 3*self.varW:
            #                     6*3**.5
            #self.w = self.avgW + (10.392304845413264 * sqrt(self.varW))*self.state.genPTR.get_next_float()
            self.w = self.avgW - sqrt(3)*sigma + 2*sqrt(3)*sigma*self.state.genPTR.get_next_float()
            return
        
        #cdef doubla a = (3*varW - avgW**2)/(3*varW + 3*avgW**2)
        #if self.state.genPTR.get_next_float() < (3*varW - avgW2)/(3*varW + 3*avgW2):
            
        #if 3*self.state.genPTR.get_next_float() < (3*self.varW - avgW2)/(self.varW + avgW2):
        if 3*self.state.genPTR.get_next_float()*(self.varW + avgW2) < (3*self.varW - avgW2):
            self.w = 0.
            return
        #self.w = self.state.genPTR.get_next_float()*(3*varW + 3*avgW**2)/2/avgW
        #self.w = 1.5*self.state.genPTR.get_next_float()*(self.varW + avgW2)/self.avgW
        self.w = 1.5*self.state.genPTR.get_next_float()*(self.varW + avgW2)/self.avgW
        return
    
    
    cdef inline void do_hinge(self):
        """
        Condensed history hinge.
        """
      # print("HINGE")
    
        #self.imfp0 = self.elastic.imfp0._eval(self.state.E)
        #cdef double rc = 1 - self.elastic.imfp._eval(self.state.E)/imfp0
        #cdef double T1, T2
       # self.T1, self.T2 = self.elastic.sampler.T(self.state.E,   1 - self.elastic.imfp._eval(self.state.E)/self.imfp0           )
        
        #print(rc, imfp0, T1, T2)
        
        #cdef double sIMFP1 = 2*T1/imfp0
        #cdef double sIMFP2 = 6*(T1 - T2)/imfp0
        
        
        
        cdef int i = self.find_index()
        cdef double sIMFP1, sIMFP2
        sIMFP1 = self.elastic.sIMFP1A[i] + self.elastic.sIMFP1B[i]*self.state.E
        sIMFP2 = self.elastic.sIMFP2A[i] + self.elastic.sIMFP2B[i]*self.state.E
        
        #REMINDER : sIMFP interpolations already include the minus sign for the exp()
        cdef double mu1 = .5*(1 - exp(self.s *sIMFP1))
        #cdef double mu2 = mu1 - (1 - exp(self.s*    6*(T1 - T2)/imfp0    ))/6
        
        
        cdef double den = (1 - 2*mu1)
        cdef double mu0 = (2*mu1 - 3*(mu1 - 0.16666666666666666*(1 - exp(self.s*sIMFP2))))/den
        
        cdef double cos
        if self.state.genPTR.get_next_float() < den + mu0: mu0 = self.state.genPTR.get_next_float()*mu0
        else:                   mu0 = mu0 + self.state.genPTR.get_next_float()* (1 - mu0)
        
        
        #cos = 1-2*mu0
        
        #if cos > 1: cos = 1
        
        #self.change_direction(cos , twoPI*self.state.genPTR.get_next_float())
        self.throwAZIMUTH()
        self.rotateTHETA(1-2*mu0)
        
    
    
    
    cdef inline void _elastic(self):
        cdef int i = self.find_index()
        #cdef DIST dist
   
        #LOGeax
        if self.state.genPTR.get_next_float()*(LOGeax[i+1] - LOGeax[i]) < ( LOGeax[i+1] - log(self.state.E) )  :
            
            #dist = self.elastic.DISTRIBUTIONS[i]
            self.mu = (<DIST>(self.elastic.DISTRIBUTIONS[i])).sample(self.state.genPTR)
        else:
            #dist = self.elastic.DISTRIBUTIONS[i + 1]
            #self.mu = dist.sample()
            self.mu = (<DIST>(self.elastic.DISTRIBUTIONS[i+1])).sample(self.state.genPTR)
        
        
        
        #self.mu = self.elastic.sample(self.state.E)
        
        #self.mu = (1-2*self.mu)
        

            
        self.throwAZIMUTH()
        self.rotateTHETA(1-2*self.mu)
        #self.change_direction(self.mu, twoPI*self.state.genPTR.get_next_float())
        
        
    cdef inline void _brem(self):
        IF _SIGNAL_INTERACTION: print("BREM")
        cdef double k, theta
        
        k, theta = self.brem.sampler.full_sample(self.state.E, self.state.genPTR) # <- dont forget secondary particle
        
        
        cdef double dE = k*self.state.E
        
        self.state.E -= dE

        if dE < photonCUTOFF:
            (<V> self.state.current_region).depositLOCAL(self.state.pos, dE)

            return

              
        cdef Photon p = Photon._new(self.state)

        p.state.E = dE


        p.throwAZIMUTH()
        p.rotateTHETA(cos(theta))
        self.nSECONDARY += 1
        self.secondary.append(p)
        
        
        
    cdef inline void _inelastic(self):
        IF _SIGNAL_INTERACTION: print("INELASTIC")
        #cdef double [::1]
        
        #self.GOS.sample(self.state.genPTR, self.find_index(), self.state.E, &particles)
        cdef int j = self.find_index()
        cdef int i
        cdef double cumul = 0
        cdef double r = self.state.genPTR.get_next_float()*self.IMFP_CUMUL.C1
        for i in range(0, self.inelastic.lenposs, 6):
            cumul += self.inelastic.arr[j, i] + self.state.E* self.inelastic.arr[j, i + 1]
            if r < cumul:
                break
            
            
        cdef double Wk = self.inelastic.arr[j, i + 2]
        cdef double Uk = self.inelastic.arr[j, i + 3]
        
        #j = i % 3
        
        #### sample secondary
        cdef Electron el
        


        cdef int _
        cdef double p2E, p2d, Qs, p2Q, a, kc
        cdef double b1, b2, b3, b4, gamma, gp1
        if Uk < self.state.E < Wk:
            if self.state.E < Uk: raise RuntimeError("this hsould not be happening")
            pass
            self.state.E -= Uk 
            # this is equivelent to assuming that the electron has given all its energy, an ejected electron with the same direction
            # 
            # if self.state.E - Uk > CUT_OFF:
            #     self.state.E -= Uk
            #     el = Electron._new(Wk-Uk, self.x, self.y, self.z,
            #                        self.eyx, self.eyy, self.eyz,
            #                        self.ezx, self.ezy, self.ezz,
            #                        self.current_region)
                

            #     self.secondary.append(el)
            #     self.nSECONDARY += 1
        
        elif i/6 % 3 == 0:
            #kc = Wk/self.state.E
            #a = (self.state.E/(self.state.E + ELECTRON_REST_ENERGY))**2
            #p2E = 1+ 5*a*kc/2
            #p2d = (1 - 2*kc)
            #p2Q = p2d/(5*a*kc)
            
            #for _ in range(100_000):
            #    r = p2E*self.state.genPTR.get_next_float()
            #    
            #    if r < 1: r = kc / (1 - r*p2d)
            #    else: r = kc + (r - 1)*p2Q
            #    
            #    if r < kc: continue
            #    if r > .5: continue
            
             #   else: P = (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(k*(1-k))))
                    
            #    if self.state.genPTR.get_next_float()*(1+5*a*r*r) < r*r*  (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(r*(1-r))  )):
            #        break
            #else: 
            #    print(">>>>>> rejection sampling took too mucch")
            #    import time
            #    print("kc", kc)
             #   print(p2E)
            #    print("a", a)
             #   print("E", self.state.E)
             #   print("Wk", Wk)
             #   print("Uk", Uk)
             #   time.sleep(1000)
                
            kc = Wk/self.state.E
            
            
            
            
            gamma = (self.state.E/ELECTRON_REST_MASS + 1)**2
            
            a = (self.state.E/(self.state.E + ELECTRON_REST_MASS))**2
            gp1 = (gamma + 1)**2
            b1 = a*(2*gp1 - 1 )/(gamma**2 - 1 )
            b2 = a*(3*gp1 + 1)/gp1
            b3 = a*2*gamma*(gamma - 1)/gp1
            b4 = a*(gamma - 1)**2 / gp1
            
            
            while 1:
                r = kc / ( 1 - self.state.genPTR.get_next_float()*(1 - kc) )
                if self.state.genPTR.get_next_float() < r*r*(1/(r*r)  - b1/r + b2 - b3*r + b4*r*r): break
                
                
            Wk = r*self.state.E
            
            
            
            self.throwAZIMUTH()
            
            if Wk - Uk > CUT_OFF:
                el = Electron._new(self.state)
                el.state.E = Wk-Uk
                el.invert_axis()
                
                Qs = (self.state.E/ELECTRON_REST_ENERGY + 1)
                
                el.rotateTHETA(sqrt(       (Wk)*(self.state.E + _2ELECTRON_REST_ENERGY)/self.state.E/(Wk + _2ELECTRON_REST_ENERGY)      )     )
                self.secondary.append(el)
                self.nSECONDARY += 1
            else:
                (<V> self.state.current_region).depositLOCAL(self.state.pos, Wk - Uk)

            self.rotateTHETA(sqrt((self.state.E - Wk)*(self.state.E + _2ELECTRON_REST_ENERGY)/self.state.E/(self.state.E - Wk + _2ELECTRON_REST_ENERGY)))
            
            self.state.E -= Wk
            # sample close
        elif i/6 % 3 == 1: # sample L far
            #self.state.E -= Wk
            p2E = self.state.E * (self.state.E + _2ELECTRON_REST_ENERGY)
            p2d = (self.state.E - Wk) * ((self.state.E - Wk) + _2ELECTRON_REST_ENERGY)   #self.p2(E - self.Wk)
            
            
            Qs = sqrt( sqrt( p2E )  - p2d  + ELECTRON_REST_ENERGY*ELECTRON_REST_ENERGY  ) - ELECTRON_REST_ENERGY
            
            
            
            
            Qs = Qs/(1 + Qs / _2ELECTRON_REST_ENERGY)
    
             
       
            Qs = Qs / ((   Qs/Wk *(1 + Wk/_2ELECTRON_REST_ENERGY )      )**(self.state.genPTR.get_next_float()) - Qs/_2ELECTRON_REST_ENERGY)
    
    
            p2Q = Qs * (Qs + _2ELECTRON_REST_ENERGY)
    
            
            
   
            
            
            self.throwAZIMUTH()
            
            if Wk - Uk > CUT_OFF:
                
                el = Electron._new(self.state)
                el.state.E = Wk-Uk
                el.invert_axis()
                
                Qs = (self.state.E/ELECTRON_REST_ENERGY + 1)
                
                el.rotateTHETA(Wk*Wk / ((Qs*Qs-1)/(Qs*Qs) / p2Q * (1 + (p2Q - Wk*Wk)/(2*Wk*(self.state.E + _2ELECTRON_REST_ENERGY)))**2))
                self.secondary.append(el)
                self.nSECONDARY += 1
            else:
                (<V> self.state.current_region).depositLOCAL(self.state.pos, Wk - Uk)

            
            self.rotateTHETA(.5*(p2E + p2d - p2Q)/sqrt(p2E*p2d))
            self.state.E -= Wk
            
            #particles.ELECTRONS.push_back(E - self.Wk)
            #particles.ELECTRONS.push_back(.5*(p2E + p2d - p2Q)/sqrt(p2E*p2d))
            
            
            
            
            
    
            
            #particles.ELECTRONS.push_back(Wk - Uk)
            #particles.ELECTRONS.push_back(Wk*Wk / ((Qs*Qs-1)/(Qs*Qs) / p2Q * (1 + (p2Q - Wk*Wk)/(2*Wk*(E + _2ELECTRON_REST_ENERGY)))**2))
            
        else:
            self.throwAZIMUTH()

            
            if Wk -  Uk > CUT_OFF:
                el = Electron._new(self.state)

                el.state.E = Wk-Uk
                el.invert_dire()
                #el.rotateTHETA(-1)
                self.secondary.append(el)
                self.nSECONDARY += 1
            else:
                (<V> self.state.current_region).depositLOCAL(self.state.pos, Wk - Uk)

            #self.rotateTHETA(1)
                
                
            self.state.E -= Wk
                
            
            #sample transverse
            
        
 
       # self.state.E -= self.inelastic.arr[j, i + 2]
      #  el = Electron._newISOTROPIC(self.inelastic.arr[j, i + 2], self.x, self.y, self.z, self.current_region, self.state.genPTR)
      #  self.nSECONDARY += 1
      #  self.secondary.append(el)
        
      
       # print(Uk)
        if Uk < MIN_CUT_OFF: 
            (<V> self.state.current_region).depositLOCAL(self.state.pos, Uk)
            return 
        
        cdef double Etot = Uk
        cdef PARTICLES particles
        (<_Atom> self.inelastic.arr_atoms[<int>self.inelastic.arr[j, i + 4]] ).run(<int>self.inelastic.arr[j, i + 5], 
                                                                       &particles,
                                                                       self.state.genPTR
                                                                       )


        
         
        cdef Photon ph
        for i in range(particles.PHOTONS.size()):
            E = particles.PHOTONS.back()
            particles.PHOTONS.pop_back()
            if E < photonCUTOFF: 
                continue
            Etot -= E
            ph = Photon._newISOTROPIC(self.state)
            ph.state.E = E
            self.nSECONDARY += 1
            self.secondary.append(ph)
        
        for i in range(particles.ELECTRONS.size()):
            E = particles.ELECTRONS.back()
            particles.ELECTRONS.pop_back()
            if E < CUT_OFF: 
                continue
            Etot -= E
 
            el = Electron._newISOTROPIC(self.state)
            el.state.E = E
            self.nSECONDARY += 1
            self.secondary.append(el)

        (<V> self.state.current_region).depositLOCAL(self.state.pos, Etot)
            
    # @cython.boundscheck(False)
    # @cython.wraparound(False) 
    # @cython.initializedcheck(False)
    # @cython.cdivision(True) 
    # cdef inline void sampleLFAR(self, double Wk):
        
    #            # cdef double dE = E - self.Wk
    #     #return sqrt( (self.momentum(E) - self.momentum(dE))**2  + ELECTRON_REST_ENERGY**2  ) \
    #           # - ELECTRON_REST_ENERGY
              
    #     cdef double p2E = E * (E + _2ELECTRON_REST_ENERGY)
    #     cdef double p2d = (E - Wk) * ((E - Wk) + _2ELECTRON_REST_ENERGY)   #self.p2(E - self.Wk)
        
        
    #     cdef double Qs = sqrt( sqrt( p2E )  - p2d  + ELECTRON_REST_ENERGY*ELECTRON_REST_ENERGY  ) - ELECTRON_REST_ENERGY
        
        
        
        
    #    # if Qm / _2ELECTRON_REST_ENERGY < 0:
    #     #    print(E, Qm, _2ELECTRON_REST_ENERGY, Qm / _2ELECTRON_REST_ENERGY, self.Wk)
        
    #     Qs = Qs/(1 + Qs / _2ELECTRON_REST_ENERGY)

         
    #     #print(_2ELECTRON_REST_ENERGY)
        

    #    # cdef double A = Qs/self.Wk *(1 + self.Wk/_2ELECTRON_REST_ENERGY )
    #     #cdef double B = Qs/_2ELECTRON_REST_ENERGY
    #     Qs = Qs / ((   Qs/self.Wk *(1 + self.Wk/_2ELECTRON_REST_ENERGY )      )**(genPTR.get_next_float()) - Qs/_2ELECTRON_REST_ENERGY)


    #     cdef double p2Q = Qs * (Qs + _2ELECTRON_REST_ENERGY)

    #   #  cdef double cos = p2E + p2d - p2Q
    #     #cos = .5*cos/(p2E*p2d)**.5
    #     particles.ELECTRONS.push_back(E - self.Wk)
    #     particles.ELECTRONS.push_back(.5*(p2E + p2d - p2Q)/sqrt(p2E*p2d))
        
        
        
        
        
    #     #cdef double Wk2 = self.Wk*self.Wk
    #     #cdef double X = (E/ELECTRON_REST_ENERGY + 1)**2
    #     Qs = (E/ELECTRON_REST_ENERGY + 1)
    #     #A = Wk2 / ((X-1)/X) / p2Q

    #     #B = (p2Q - Wk2)/(2*self.Wk*(E + _2ELECTRON_REST_ENERGY))
        
    #     particles.ELECTRONS.push_back(Wk - Uk)
    #     particles.ELECTRONS.push_back(Wk*Wk / ((Qs*Qs-1)/(Qs*Qs) / p2Q * (1 + (p2Q - Wk*Wk)/(2*Wk*(E + _2ELECTRON_REST_ENERGY)))**2))
        
    #     #cdef double cos_sec = A * (1 + B) **2
        #return (newE, cos, Esec, cos_sec)

        
        
        
        # SECONDARY ELECTRON
        # cos2 = (W)*(E + _2ELECTRON_REST_ENERGY)/E/(W + _2ELECTRON_REST_ENERGY)
        # particles.ELECTRONS.push_back(self.Wk - self.Uk)
        # particles.ELECTRONS.push_back(sqrt(cos2))

        
        # self.state.E = particles.ELECTRONS[0]
        
        # if self.state.E == 0:

        #     el = Electron._new(particles.ELECTRONS[2], 
        #                         self.x, self.y, self.z,
        #                         self.eyx, self.eyy, self.eyz,
        #                         self.ezx, self.ezy, self.ezz,
        #                         self.current_region)
            
        #     #self.rotateTHETA(particles.ELECTRONS[1])
        #     #el.rotateTHETA(particles.ELECTRONS[3])
            
            
        #     self.nSECONDARY += 1
        #     self.secondary.append(el)
        #     return
        
        # #cdef double cos = particles.ELECTRONS[1]
        
        # self.throwAZIMUTH()
        
        # if particles.ELECTRONS[2] < CUT_OFF:
        #     self.rotateTHETA(particles.ELECTRONS[1])
        #     return 
        

        # el = Electron._new(particles.ELECTRONS[2], 
        #                     self.x, self.y, self.z,
        #                     -self.eyx, -self.eyy, -self.eyz,
        #                     self.ezx, self.ezy, self.ezz,
        #                     self.current_region)
        
        # self.rotateTHETA(particles.ELECTRONS[1])
        # el.rotateTHETA(particles.ELECTRONS[3])
        
        
        # self.nSECONDARY += 1
        # self.secondary.append(el)
    
        
        
        
        
        
            
     #   print("INELASTIC")
     
     
     
     
     
     
     
        #cdef double cos, Esec, cos_sec
        
        ### RUNNING GOS MODEL
        #self.GOS.update(self.state.E, False)
        
        #self.state.E, self.cos, self.state.Esec, self.cos_sec = self.GOS.sample(self.current_region, self.x, self.y, self.z)
        
        # if self.cos>1:
        #     #print("cos > 1")
        #     self.cos = 1

        #self.throwAZIMUTH()
       # self.rotateTHETA(self.cos)
        #self.change_direction(self.cos, twoPI*self.state.genPTR.get_next_float())

        ## RUNNING pyRelax
        
        #atom.ionize(atom.SHELLS[I])
        
        
        
        
        
        
        
        #self.secondary.extend(self.GOS.secondary)
        #self.nSECONDARY += self.GOS.nSECONDARY


        
        #if self.state.Esec < CUT_OFF:
        #    return
        
        
        
        # if self.cos_sec > 1: 
        #     self.cos_sec = 1
            
            
        
        
        #cdef Electron p = Electron._new(self.state.Esec, self.x, self.y, self.z,
         #                             self.eyx, self.eyy,self.eyz,
         #                             self.ezx, self.ezy,self.ezz, 
          #                            self.current_region)
        #p.genPTR = self.state.genPTR
        #p.throwAZIMUTH()
        #p.rotateTHETA(self.cos_sec)
        #self.nSECONDARY += 1
        #self.secondary.append(p)
                
        
        
        
        
        
        # axis = ez0
        # ey = ey0.rotateAngle(axis, twoPI*self.state.genPTR.get_next_float())
        
        # axis = ey
        # ez = ez0.rotateCos(axis, self.cos_sec)  
            
            
            
            
            
            
            
            
            
        # self.secondary.append(
        #     Electron._new(
        #        self.current_region,
        #        self.state.Esec, 
        #        self.pos, 
        #        ey,
        #        ez, 100)
        #     )
        
        
        

        
 

        
    
    def getindex(self):
        return self.find_index()
    
    def setE(self, double E):
        self.state.E = E
    
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
            MID = START + (END - START)//2 
            
            xMID = EAX[MID]
            
            if self.state.E is xMID: #found the value
                return MID
            
            if self.state.E < xMID: # discard right side
                END = MID - 1 # do not include mid
                continue
            
            START = MID + 1
        return END 
                    
            
            

    cdef void _anihilation(self):
        cdef double gamma = self.state.E/0.511e6 + 1
        cdef double f_min = 1/(gamma + 1 + (gamma**2 - 1)**.5)
        cdef double v, r
        while 1:
            v = f_min * ((1-f_min)/f_min)**self.state.genPTR.get_next_float()
            if self.state.genPTR.get_next_float()*g(f_min, gamma) <= g(v, gamma):
                break
            
            
            

            
        v = min(v, 1-v)
        cdef double Eplus  = (1 - v)*(self.state.E + 2*0.511e6)
        cdef double Eminus = v*(self.state.E + 2*0.511e6)
        

        self.throwAZIMUTH()
        cdef Photon p = Photon._new(self.state)
        p.state.E = Eplus
        
        self.secondary.append(p)
        self.nSECONDARY += 1
        
        r = 1/sqrt(gamma**2 - 1)
        p.rotateTHETA(r*(gamma + 1 - 1/v))
        
        p = Photon._new(self.state)
        p.invert_axis()
        p.state.E = Eminus

        p.rotateTHETA(r*(gamma + 1 - 1/(1-v)))
        
        self.secondary.append(p)
        self.nSECONDARY += 1





cdef mixmax_engine GEN # space to store a generator for the python_hooks.Positron

class python_hooks:
    class Positron(Positron):
        """
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
            if attribute == "E":                return          (<Positron> self).state.E
            if attribute == "current_material": return  <MAT> ( (<Positron> self).current_material )
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
            if   attribute == 'E':                 (<Positron> self).state.E = value
            elif attribute == "current_material":  (<Positron> self).current_material = <void*> value
            elif attribute == "current_region":    (<Positron> self).state.current_region   = <void*> value
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

        def update_references(self):    (<Positron> self).update_references()
        def update_imfp(self):          (<Positron> self).update_imfp()
        def update_imfp_cumul(self):    (<Positron> self).update_imfp_cumul()

        def sample_w(self, double tau): (<Positron> self).sample_w(tau)
        def do_hinge(self):             (<Positron> self).do_hinge()

        def _anihilation(self):         (<Positron> self)._anihilation()
        def _elastic(self):             (<Positron> self)._elastic()
        def _brem(self):                (<Positron> self)._brem()
        def _inelastic(self):           (<Positron> self)._inelastic()
        def _delta(self):               (<Positron> self)._delta()




        def find_index(self): return (<Positron> self).find_index()

        def __repr__(self):
            return "<python_hook.Photon>"

        def __str__(self):
            return "RETURN DEBUG INFO"

