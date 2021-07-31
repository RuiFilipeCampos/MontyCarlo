# cython: profile=False
# cython: annotate=True
# distutils: language = c++ 
# distutils: extra_compile_args = -std=c++11

print("Importing .particles.electrons")

DEF _DEBUG_BASIC = False
DEF _SIGNAL_INTERACTION = False


DEF RECORD = True
from ..materials.cppRelaxAPI cimport PARTICLES


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



from libc.math cimport isnan


from .particle cimport STATE




#Error messages (to be moved to its own module)
errorMSG1 = "Exhausted allowed number of iterations for rejection sampling."


## PYTHON IMPORTS
#Local Imports
from ..materials import database as db
from ..settings import __photonCUTOFF__, __electronCUTOFF__

#external imports
from collections import deque
import numpy as np

## CYTHON IMPORTS
#Local Imports
from .particle cimport Particle
from ..geometry.main cimport Volume
from ..tools.vectors cimport Vector
from .photons cimport Photon
from libcpp.vector cimport vector



from ..materials.materials cimport Material




print("Size of mater", sizeof(Material))
cdef extern from "<math.h>" nogil:
    double frexp(double x, int* exponent)


# cdef int get_exp(double x):
#     cdef int exp;
#     frexp(x, &exp);
#     return exp;


from ..materials.materials cimport Material

from ..materials.electron.main cimport Brem, Inelastic, Elastic, DIST



from .._init  import eax
from .._init  cimport EAX
from .._init  cimport LIMS


cdef double[::1] LOGeax = np.log(eax)
cdef double[::1] diffLOGeax = np.diff(np.array(LOGeax))


# cdef double[:] eax = _eax

# 2695



#from ..materials.electron.main cimport LIMS 



# print(len(eax))
# import time
# time.sleep(100000)
# External Imports
from libc.math cimport sin, cos, log, sqrt, pi, exp, fmin, fmax, acos, pow

cdef double twoPI = 2*pi




cimport cython



cdef double DMAX = .01

# ---- RANDOM SAMPLER --- #######################



from libc.stdlib cimport rand, RAND_MAX, srand

from ..types cimport double3



from .._random.interface cimport mixmax_engine


cdef extern from "math.h":
    float INFINITY

#from numpy.random import rand as urand
###################################################################

# CONSTANTS AND GLOBALS
cdef double CUTOFF = __photonCUTOFF__
cdef double photonCUTOFF = __photonCUTOFF__
cdef double MIN_CUT_OFF = min(__photonCUTOFF__, __electronCUTOFF__)

cdef double E0_el = db.E0_electron*1e-3


ctypedef Volume V
    
    
cdef double CUT_OFF = __electronCUTOFF__
#@cython.cdivision(True)
cdef double ELECTRON_REST_ENERGY = 0.51099895000*1e6 #eV
cdef double  _2ELECTRON_REST_ENERGY    = 2 *ELECTRON_REST_ENERGY
@cython.boundscheck(False)
@cython.wraparound(False) 
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef class Electron(Particle):
 
    # @staticmethod
    # def new(Volume current_region, double E):
    
    #     return Electron._newISOTROPIC( current_region, E, Vector(0, 0, 0))
    
    
    
    cdef double ENERGY(self):
        return self.state.E
    



    @staticmethod
    cdef Electron _new(STATE& state):
        """
        A new Electron is created with its state set equal to the
        provided STATE reference.
        """


        cdef Electron self
        self = <Electron>Electron.__new__(Electron)
        self.state = state
        return self
    
    @staticmethod
    cdef Electron _newISOTROPIC(STATE& state):
        """
        A particle with a random direction in the unit sphere is generated.
        The direction is randomly sampled here and the axis is assumed to be
        thrown in the next interaction.

        The remaining state is copied from the provided STATE reference.
        """

        cdef Electron self
        self = <Electron>Electron.__new__(Electron)
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
        return self    



    ####################################################################################
    ########                           RUN                                      ########
    ########                           RUN                                      ########
    ####################################################################################

    cdef void _run(Electron self, mixmax_engine* genPTR):
        IF _DEBUG_BASIC: print("> ELECTRON")
        self.secondary = deque()
        self.MU = deque()
        self.nSECONDARY = 0
        cdef bint deposit_in_2 = True
        cdef double E0
        cdef double r
        if self.state.E < CUT_OFF:
            (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
            (<Volume> self.state.current_region).exit()
            return

        self.state.genPTR = genPTR
        
        self.update_references()

        IF RECORD: self.record()
        
        #cdef double r
        cdef double  tau, S_soft
        cdef bint delta = False
        cdef double tau2 
        cdef double3 dep_pos
        while True:
            if self.state.pos.x**2 + self.state.pos.y**2 + self.state.pos.z**2 > 10_000**2:
                return

            self.s = -log(self.state.genPTR.get_next_float())/self.imfp_max

            if self.s > self.s_max:
                self.s = self.s_max
                delta = True
            
            if self.s_max > 1000:
                print((<V> self.state.current_region))
         
            

            
            # FIRST PART OF TRAJECTORY
            
            #to_deposit = self.s*self.state.genPTR.get_next_float()
            self.sample_w(self.s)

            #if to_deposit < tau:
            #    dep_pos.x = self.state.pos.x + to_deposit*self.state.dire.x
            ##    dep_pos.y = self.state.pos.y + to_deposit*self.state.dire.y
            #    dep_pos.z = self.state.pos.z + to_deposit*self.state.dire.z
            #    (<V> self.state.current_region).depositLOCAL(dep_pos, self.w)
            #    deposit_in_2 = False


            S_soft = self.w/self.s
            tau = self.s*self.state.genPTR.get_next_float()
            self.state.L = tau

            # first segment, pos0
            if (<V> self.state.current_region).move(self.state, S_soft):
                self.state.E -= (tau - self.state.L)*S_soft
                if self.state.E < CUT_OFF:
                    (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
                    return
                (<V> self.state.current_region).depositLOCAL(self.state.pos, (tau - self.state.L)*S_soft)
                self.update_references()
                continue

            #(<V> self.state.current_region).depositRANDOM(self.state, S_soft*tau, tau)

            self.state.E -= S_soft*tau
            if self.state.E <  CUT_OFF:
               # if deposit_in_2:
                (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
                (<Volume> self.state.current_region).exit()
                return
            



            self.do_hinge()


            #(<V> self.state.current_region).depositLOCAL(self.state.pos, self.w)
            


            IF RECORD: self.record()
   

           # if deposit_in_2:
           #     to_deposit = to_deposit - tau
          #     dep_pos.x = self.state.pos.x + to_deposit*self.state.dire.x
          #      dep_pos.y = self.state.pos.y + to_deposit*self.state.dire.y
          #      dep_pos.z = self.state.pos.z + to_deposit*self.state.dire.z
           #     (<V> self.state.current_region).depositLOCAL(dep_pos, self.w)
           # else: deposit_in_2 = True


            tau = self.s - tau
            self.state.L = tau

            

            

            if (<V> self.state.current_region).move(self.state, S_soft):
                self.state.E -= (tau - self.state.L)*S_soft
                if self.state.E < CUT_OFF:
                   # (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E + (tau - self.state.L)*S_soft)
                    (<Volume> self.state.current_region).exit()
                    return
                

                self.update_references()
                continue
            self.state.E -= S_soft*tau

          #  (<V> self.state.current_region).depositRANDOM(self.state, S_soft*tau, tau)

            if self.state.E < CUT_OFF:
                (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
                (<Volume> self.state.current_region).exit()
                return
            

            IF RECORD: self.record()
                        
            #if self.state.E < CUT_OFF:
            #    (<V> self.state.current_region).depositLOCAL(self.state.pos, self.state.E)
            #    (<Volume> self.state.current_region).exit()
            #    return

            if delta:
                self.update_imfp()
                delta = False
                continue
     
            self.update_imfp_cumul()
            
            r = self.state.genPTR.get_next_float()*self.imfp_max
            
            
            
            if r < self.IMFP_CUMUL.C1: 
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
                continue


            
            
            
            
            
 
            

    










    ####################################################################################
    ########                           UPDATE                                   ########
    ########                           METHODS                                  ########
    ####################################################################################

    cdef void update_references(self):
        """
        Updates all references. Called when there is a region crossing.
        """

        #getting material from current region
        self.current_material = (<Volume> self.state.current_region).material
        self.GOS = self.current_material.electron.inelastic.gosMOLECULE
        self.electron = self.current_material.electron

        #these references are used by their corresponding _fooInteraction method
        self.elastic    = self.electron.elastic
        self.inelastic  = self.electron.inelastic
        self.brem       = self.electron.brem

        #since region crossing has ocurred, update the inverse mean free paths
        self.update_imfp()



    cdef inline void update_imfp_cumul(Electron self):
        """
        Updates hard inverse mean free paths (hIMFP). 
        Constructs corresponding cumulutative table.
        Sets s_max and corresponding maximum energy loss value w_max.
        Sets upper boundary for the total hIFMP (self.imfp_max).
        Residual probability results in a 'delta interaction' in which nothing at all happens.


        Called when a new interaction is to be chosen.
        """

        cdef int i = self.find_index()  
 
        
        self.IMFP_CUMUL.C1 = self.inelastic.imfpA[i] + self.inelastic.imfpB[i]*self.state.E
        self.IMFP_CUMUL.C2 = self.IMFP_CUMUL.C1 + self.brem.imfpA[i] +  self.brem.imfpB[i]   *self.state.E
        self.IMFP_CUMUL.C3 = self.IMFP_CUMUL.C2 + self.elastic.imfpA[i] + self.elastic.imfpB[i] * self.state.E
        #self.electron.IMFP(self.state.E)
        self.s_max = 4/self.IMFP_CUMUL.C3 * (.5+.5*self.state.genPTR.get_next_float())
        self.wmax = self.electron.find_wmax(self.s_max, self.state.E)
        cdef double Emin = self.state.E - self.wmax
     #   print("wmax = ", wmax, "E = ", self.state.E)
        
        
        if Emin is 0:
            self.imfp_max = self.IMFP_CUMUL.C3
            return
        
        i = self.electron.find_index(Emin)
        
        
        self.imfp_max = fmax(self.IMFP_CUMUL.C3, 
                             self.electron.imfpA[i] + self.electron.imfpB[i] * Emin
                             )
        
        
        
    cdef inline void _delta(self):
        """
        TO BE REMOVED.
        """
        pass

        
    cdef inline void update_imfp(Electron self):
        """
        Updates total hIMFP. Called when a new proposal displacement is needed
        but there is no need to calculate the cumulutative table since no interaction
        is going to occur. 
        """

      #  print("update_imfp")

        cdef int i, j

        j = self.find_index()
        cdef double imfp = self.electron.imfpA[j] + self.electron.imfpB[j] * self.state.E
        self.s_max = 4/imfp * (.5+.5*self.state.genPTR.get_next_float())


        self.wmax =  self.electron.find_wmax(self.s_max, self.state.E)
        cdef double Emin = self.state.E - self.wmax
        
        
        if Emin is 0:
            print("Emin is 0 in:")
            print((<V> self.state.current_region))
            self.imfp_max = imfp
            #j = self.find_index()
            #self.imfp_max = self.electron.imfpA[j] + self.electron.imfpB[j] * self.state.E
            #self.imfp_max *= 1.1
            #self.s_max = 4/self.imfp_max
            return
        
        i = self.electron.find_index(Emin)
        
        
        self.imfp_max = fmax(imfp, 
                             self.electron.imfpA[i] + self.electron.imfpB[i] * Emin
                             )



        
        #i = self.electron.find_index(Emin)
        #j = self.find_index()
        ##cdef double imfpH = 
   #
        #self.s_max = self.electron.imfpA[j] + self.electron.imfpB[j] * self.state.E
        #self.imfp_max = fmax(self.s_max,
        #                     self.electron.imfpA[i] + self.electron.imfpB[i] * Emin
        #                     )
#
        #self.s_max = 4/self.imfp_max

    def UT_find_wmax(self, double s, double E):
        print(f"""
Given a value of s = {s}, find the maximum energy loss in accordance to the soft stopping
power. 
            """)

        self.state.E = E
        self.update_references()



        cdef int i0 = self.electron.find_index(self.state.E)
        print(f"Found index i0 = {i0} such that {EAX[i0]} < {self.state.E} < {EAX[i0 + 1]}")

        cdef wmax = self.electron.find_wmax(s, self.state.E)
        print(f"Found value of wmax = {wmax}eV")
        cdef double Emin = self.state.E - wmax
        print(f"Minimum energy value is {Emin}eV")
        cdef int i1 = self.electron.find_index(self.state.E - wmax)
        print(f"Found index i1 = {i1} such that {EAX[i1]} < {Emin} < {EAX[i1 + 1]}")


        print("Building array of sSP values in the given range.")

        import numpy as np
        sSP = self.electron.softSP

        from scipy.interpolate import CubicSpline
        sSP = CubicSpline(eax, 1/sSP)


        new_s = sSP.integrate(self.state.E, Emin)
        print(f"Integration found value of s = {new_s}cm")
        print(f"Provided value of s = {s}cm")
        
        

    ####################################################################################
    ########                           RANDOM                                   ########
    ########                           SAMPLING                                  ########
    ####################################################################################
    
    
    
    cdef inline void sample_w(self, double tau):
        """
        Sample energy loss along the proposed displacement.
        An artificial distribution is used, constructed in such a way that
        the first and second moments of the true energy loss distribution
        are reproduced. 
        """
        
        #self.SP     = self.inelastic.softSP._eval(self.state.E) + self.brem.softSP._eval(self.state.E)
        #self.STRAGG = self.inelastic.softSTRAGG._eval(self.state.E) + self.brem.softSTRAGG._eval(self.state.E)

        #self.avgW   = self.s * (self.inelastic.softSP._eval(self.state.E)     + self.brem.softSP._eval(self.state.E)    )
        #self.varW   = self.s * (self.inelastic.softSTRAGG._eval(self.state.E) + self.brem.softSTRAGG._eval(self.state.E))
        
        cdef int i = self.find_index()

        cdef double SP = ( self.electron.softSPA[i]     + self.electron.softSPB[i]*self.state.E)
        cdef double STRAGG = ( self.electron.softSTRAGGA[i] + self.electron.softSTRAGGB[i]*self.state.E)


        self.avgW = tau * SP * (1. - .5*self.electron.softSPB[i]*tau)
        self.varW = tau * STRAGG - tau*tau*(.5*self.electron.softSTRAGGB[i]*SP + STRAGG*self.electron.softSPB[i])



        cdef double avgW2 = self.avgW*self.avgW
        cdef double sigma = sqrt(self.varW)
        if avgW2 > 9.*self.varW:
            #1.015387*
            self.w = sigma*self.current_material.electron.gauss._sample() + self.avgW #CONFIRM
         #   print("sampled w = ", self.w, "FROM TRUNCATED GAUSSIAN")

            return
        
        if avgW2 > 3.*self.varW:
            #                     6*3**.5
            #self.w = self.avgW + (10.392304845413264 * sqrt(self.varW))*self.state.genPTR.get_next_float()
            self.w = self.avgW - sqrt(3.)*sigma + 2.*sqrt(3.)*sigma*self.state.genPTR.get_next_float()
         #   print("sampled w = ", self.w, "FROM SECOND CASE")

            return
        
        #cdef doubla a = (3*varW - avgW**2)/(3*varW + 3*avgW**2)
        #if self.state.genPTR.get_next_float() < (3*varW - avgW2)/(3*varW + 3*avgW2):
            
        #if 3*self.state.genPTR.get_next_float() < (3*self.varW - avgW2)/(self.varW + avgW2):
        if 3.*self.state.genPTR.get_next_float()*(self.varW + avgW2) < (3.*self.varW - avgW2):
            self.w = 0.
           # print("sampled w = ", self.w)
            return
        #self.w = self.state.genPTR.get_next_float()*(3*varW + 3*avgW**2)/2/avgW
        #self.w = 1.5*self.state.genPTR.get_next_float()*(self.varW + avgW2)/self.avgW
        self.w = 1.5*self.state.genPTR.get_next_float()*(self.varW + avgW2)/self.avgW
       # print("sampled w = ", self.w, "FROM THIRD  CASE")

        return
    




    def UNITTESTING_ENERGYLOSS(self, int N, E = 1e6):
        print(f"""

            """)
        self.state.E = E
        self.update_references()

        print("MAXIMUM ENERGY LOSS ALLOWED: ", self.wmax)
        print("MAXIMUM STEP LENGHT: ", self.s_max)
        print("MAX IMFP_h: ", self.imfp_max)

        import numpy as np
        SAMPLE_s = np.zeros(N)
        SAMPLE_w = np.zeros(N)

        self.s = -log(self.state.genPTR.get_next_float())/self.imfp_max

        if self.s > self.s_max:
            self.s = self.s_max

        print(f"SAMPLED S: {self.s}")


        cdef int i
        for i in range(N):
            self.sample_w(self.s)
            SAMPLE_w[i] = self.w

        cdef double avg = np.mean(SAMPLE_w)
        cdef double var = np.std(SAMPLE_w)**2
        i = self.find_index()
        print("EXPECTED MEAN (Ss*s): ", self.s*( self.electron.softSPA[i]     + self.electron.softSPB[i]*self.state.E)  )
        print("SAMPLED MEAN:", avg)


        print("EXPECTED VAR(stragg*s)", self.s*( self.electron.softSTRAGGA[i] + self.electron.softSTRAGGB[i]*self.state.E))
        print("SAMPLED VAR:", var)
        return SAMPLE_w









    
    def UNITTESTING_HINGE(self, int N, E = 1e6):
        print(f"""
The sampling of the angular deflection of the hinge is made from an artificial distribution
that was taylor made for reproducing the first and second moments of the actual soft scattering
distribution.

ENERGY: {E}eV

            """)
        self.state.E = E
        self.update_references()

        self.s = -log(self.state.genPTR.get_next_float())/self.imfp_max

        if self.s > self.s_max:
            self.s = self.s_max

        print(f"s = {self.s}")
        cdef int i = self.find_index()


        cdef double sIMFP1 = self.elastic.sIMFP1A[i] + self.elastic.sIMFP1B[i]*self.state.E
        cdef double sIMFP2 = self.elastic.sIMFP2A[i] + self.elastic.sIMFP2B[i]*self.state.E


        print(f"sIMFP1 = {sIMFP1}, sIMFP2 = {sIMFP2}")
        cdef double mu1 = .5*(1. - exp(self.s *sIMFP1))
        cdef double mu2 = mu1 -  (1. - exp(self.s*sIMFP2))/6.
        print("VALUES TO BE REPRODUCED: ")
        print(f"<mu> = {mu1} | <mu**2> = {mu2}")


        import numpy as np
        SAMPLE_MU = np.zeros(N)
        SAMPLE_S = np.zeros(N)
        for i in range(N):



            self.do_hinge()
            SAMPLE_MU[i] = self.mu
            SAMPLE_S[i] = self.s

        avg = np.mean(SAMPLE_MU)
        avg2 = np.mean(SAMPLE_MU**2)
        print(f"SAMPLE MEAN: {avg} | error = {100*(avg - mu1)/mu1}%")
        print(f"SAMPLE SEC MOMENT: {avg2} | error = {100*(avg2 - mu2)/mu2}%")

        return {"mu":SAMPLE_MU, "s":SAMPLE_S}



    cdef inline void do_hinge(self):
        """
        Sample angular defflection of the hinge. 
        Sampling is done from an artificial distribution constructed
        in such a way that it reproduces the first and second moments
        of the true distribution.
        """
      # print("HINGE")
    
        #self.imfp0 = self.elastic.imfp0._eval(self.state.E)
        #cdef double rc = 1 - self.elastic.imfp._eval(self.state.E)/imfp0
        #cdef double T1, T2
       # self.T1, self.T2 = self.elastic.sampler.T(self.state.E,   1 - self.elastic.imfp._eval(self.state.E)/self.imfp0           )
        
        #print(rc, imfp0, T1, T2)
        
        #cdef double sIMFP1 = 2*T1/imfp0
        #cdef double sIMFP2 = 6*(T1 - T2)/imfp0
        
      #  input(f"STARTING HINGE | E = {self.state.E} | s = {self.s}")
        
        cdef int i = self.find_index()
       # input(f"found index {i}, corresponding to {EAX[i]} < {self.state.E} < {EAX[i+1]}")

        cdef long double sIMFP1, sIMFP2
        sIMFP1 = self.elastic.sIMFP1A[i] + self.elastic.sIMFP1B[i]*self.state.E
        sIMFP2 = self.elastic.sIMFP2A[i] + self.elastic.sIMFP2B[i]*self.state.E
       # input(f"sIMFP1 = {sIMFP1} | sIMFP2 = {sIMFP2}")




        #REMINDER : sIMFP interpolations already include the minus sign for the exp()
       # input("mu1 = .5*(1 - exp(self.s *sIMFP1))")
        cdef long double mu1 = .5*(1. - exp(self.s *sIMFP1))
       # print(mu1)
        #cdef double mu2 = mu1 - (1 - exp(self.s*    6*(T1 - T2)/imfp0    ))/6
       # print("mu1 = ", mu1)
      #  input("den = (1 - 2*mu1)")
        cdef long double den = (1. - 2.*mu1)
       # print(den)

       # input(f"mu2 = mu1 - 1/6 * (1 - exp({self.s*sIMFP2}))")
        cdef long double mu2 = mu1 -  (1. - exp(self.s*sIMFP2))/6.
       # print(mu2)

       # print("mu2 = ", mu2)

       # input("b = (2*mu1 - 3*mu2)/den")
        cdef long double b = (2.*mu1 - 3.*mu2)/den
      #  print(b)


       #"" print("b = ", b)

      #  cdef double mu0 = (2*mu1 - 3*(mu1 - 0.16666666666666666*(1 - exp(self.s*sIMFP2))))/den
        
        if self.state.genPTR.get_next_float() < den + b:
           # input("CASE 1")
            self.mu = self.state.genPTR.get_next_float()*b
          #  print(f"mu = {self.mu}")
        else:    
          #  input("CASE 2")
            self.mu = b + self.state.genPTR.get_next_float()* (1. - b)
           # print(f"mu  = {self.mu}")
        
        #cos = 1-2*mu0
        
        #if cos > 1: cos = 1
        
        #self.change_direction(cos , twoPI*self.state.genPTR.get_next_float())
        self.throwAZIMUTH()
      #  print("HINGE: ", mu0, self.state.E)
        self.rotateTHETA(1.-2.*self.mu)
        
    
    def UT_elastic(self, int N, double E):
        self.state.E = E
        self.update_references()

        cdef int i = self.find_index()
        #cdef DIST dist
   
        print("CUT OFF VALUES OF MU:")
        print(  (<DIST>(self.elastic.DISTRIBUTIONS[i])).mu_c  )

        print(  (<DIST>(self.elastic.DISTRIBUTIONS[i+1])).mu_c ) 



        import numpy as np
        sample = np.zeros(N)

        for i in range(N):
            self._elastic()
            sample[i] = self.mu
        
        return sample
    
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
        
     #   print("ELASTIC: ", self.mu, self.state.E)

        #self.mu = self.elastic.sample(self.state.E)
        
        #self.mu = (1-2*self.mu)
        

            
        self.throwAZIMUTH()
        self.rotateTHETA(1-2*self.mu)
        #self.change_direction(self.mu, twoPI*self.state.genPTR.get_next_float())
        
        
    cdef inline void _brem(self):
        """Simulate the bremstrahlung interaction.
        """
        
        ## Debugging Flags
        IF _SIGNAL_INTERACTION: print("BREM")
            
            
        ## Implementation
        cdef double k       # Fraction of energy that the electron has lost.
        cdef double theta   # Angular deflection with respect to the electrons direction, in which the secondary photon will be emitted. 
        
        k, theta = self.brem.sampler.full_sample(self.state.E, self.state.genPTR)
        
        cdef double dE      # Energy that the electron has lost. 
        
        dE = k*self.state.E
        self.state.E -= dE

        if dE < photonCUTOFF:
            (<V> self.state.current_region).depositLOCAL(self.state.pos, dE)
            return
    
        cdef Photon p       # The secondary photon.
        
        p = Photon._new(self.state)
        
        p.state.E = dE
        p.throwAZIMUTH()
        p.rotateTHETA(cos(theta))
        self.nSECONDARY += 1
        self.secondary.append(p)
        

        
 
        
    cdef inline void _inelastic(self):
        """
        INELASTIC COLLISION
            GOS MODEL

            CHOOSE SHELL AMONG ALL SHELLS IN THE MOLECULE
            CHOOSE INTERACTION REGIME (FAR LONGITUDINAL, FAR TRANSVERSE, CLOSE)
            SAMPLE INTERACTION REGIME -> EFFECT ON PROJECTILE AND EJECTED ELECTRON
            SIMULATE RELAXATION EFFECTS -> CALL pyRELAX
        """


        IF _SIGNAL_INTERACTION: print("INELASTIC")
        #cdef double [::1]
        
        #self.GOS.sample(self.state.genPTR, self.find_index(), self.state.E, &particles)

      #  input("STARTING _inelastic")
      #  input(f"current energy: {self.state.E}eV")


        cdef int j = self.find_index()
      #  input(f"found index {j}")

        cdef int i
        cdef double cumul = 0
        cdef double r = self.state.genPTR.get_next_float()*self.IMFP_CUMUL.C1
       # input(f"threw r = {r} in {self.IMFP_CUMUL.C1}")


        for i in range(0, self.inelastic.lenposs, 6):
            cumul += self.inelastic.arr[j, i] + self.state.E* self.inelastic.arr[j, i + 1]
          #  input(f"{r} < {cumul}?")
            if r < cumul:
                #print("yes")
                break
            
            
        cdef double Wk = self.inelastic.arr[j, i + 2]
        cdef double Uk = self.inelastic.arr[j, i + 3]
        #input(f"CHOSEN SHELL: Wk = {Wk}eV | Uk = {Uk} eV")

        #j = i % 3
        
        #### sample secondary
        cdef Electron el
        


        cdef int _
        cdef double p2E, p2d, Qs, p2Q, a, kc
        
        if Uk < self.state.E < Wk:
            if self.state.E < Uk: raise RuntimeError("this hsould not be happening")
            #input(f"CASE 1, JUST REDUCING ENERGY BY Wk-Uk")

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
        
        elif i/6 % 3 == 0: # hard close
           # input(f"SAMPLING HARD CLOSE COLLISION")
            kc = Wk/self.state.E
            a = (self.state.E/(self.state.E + ELECTRON_REST_ENERGY))**2
            p2E = 1+ 5*a*kc/2
            p2d = (1 - 2*kc)
            p2Q = p2d/(5*a*kc)
            
            for _ in range(100_000): 
                r = p2E*self.state.genPTR.get_next_float()
                
                if r < 1: r = kc / (1 - r*p2d)
                else: r = kc + (r - 1)*p2Q
                
                if r < kc: continue
                if r > .5: continue
            
             #   else: P = (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(k*(1-k))))
                    
                if self.state.genPTR.get_next_float()*(1+5*a*r*r) < r*r*  (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(r*(1-r))  )):
                    break
            else: 
                print(">>>>>> rejection sampling took too mucch")
                import time
                print("kc", kc)
                print(p2E)
                print("a", a)
                print("E", self.state.E)
                print("Wk", Wk)
                print("Uk", Uk)
                time.sleep(1000)
                
                
                
                
            Wk = r*self.state.E
         #   input(f"SAMPLED ENERGY LOSS: {Wk}/{self.state.E}")
            
            
            self.throwAZIMUTH()
            
            if Wk - Uk > CUT_OFF:
                el = Electron._new(self.state)

                el.invert_axis()
                el.state.E = Wk-Uk

                
                Qs = (self.state.E/ELECTRON_REST_ENERGY + 1)
                
                el.rotateTHETA(sqrt(       (Wk)*(self.state.E + _2ELECTRON_REST_ENERGY)/self.state.E/(Wk + _2ELECTRON_REST_ENERGY)      )     )
                self.secondary.append(el)
                self.nSECONDARY += 1
            else:
                (<V> self.state.current_region).depositLOCAL(self.state.pos, Wk - Uk)
            self.rotateTHETA(sqrt((self.state.E - Wk)*(self.state.E + _2ELECTRON_REST_ENERGY)/self.state.E/(self.state.E - Wk + _2ELECTRON_REST_ENERGY)))
           # input(f"ROTATING PROJECTILE BY COS = {sqrt((self.state.E - Wk)*(self.state.E + _2ELECTRON_REST_ENERGY)/self.state.E/(self.state.E - Wk + _2ELECTRON_REST_ENERGY))}")

            self.state.E -= Wk
            # sample close
        elif i/6 % 3 == 1: # sample L far
            #self.state.E -= Wk
            #input(f"SAMPLING FAR LONGITUDINAL")
            p2E = self.state.E * (self.state.E + _2ELECTRON_REST_ENERGY)
            #input(f"p2E = {p2E}")

            p2d = (self.state.E - Wk) * ((self.state.E - Wk) + _2ELECTRON_REST_ENERGY)   #self.p2(E - self.Wk)
            #input(f"p2d = {p2d}")

            
            Qs = sqrt( (sqrt( p2E )  - sqrt(p2d))**2  + ELECTRON_REST_ENERGY*ELECTRON_REST_ENERGY  ) - ELECTRON_REST_ENERGY
            #input(f"Qs = {Qs}")

            
            
            
            Qs = Qs/(1 + Qs / _2ELECTRON_REST_ENERGY)
            #input(f"Qs = {Qs}")
             
       
            Qs = Qs / ((   Qs/Wk *(1 + Wk/_2ELECTRON_REST_ENERGY )      )**(self.state.genPTR.get_next_float()) - Qs/_2ELECTRON_REST_ENERGY)
            #input(f"Qs = {Qs}")


         #   print("IF THIS IS NOT ZERO, MUST CHANGE TO DEPOSIT THIS ENERGY", Qs - Wk - Uk)


            p2Q = Qs * (Qs + _2ELECTRON_REST_ENERGY)
    
            
            
   
            
            
            self.throwAZIMUTH()
            
            if Wk - Uk > CUT_OFF:
                
                el = Electron._new(self.state)
                el.invert_axis()
                el.state.E = Wk-Uk
                
                Qs = (self.state.E/ELECTRON_REST_ENERGY + 1)
                
                el.rotateTHETA(Wk*Wk / ((Qs*Qs-1)/(Qs*Qs) / p2Q * (1 + (p2Q - Wk*Wk)/(2*Wk*(self.state.E + _2ELECTRON_REST_ENERGY)))**2))
                self.secondary.append(el)
                self.nSECONDARY += 1
            else:
                (<V> self.state.current_region).depositLOCAL(self.state.pos,Wk - Uk)

            
            self.rotateTHETA(.5*(p2E + p2d - p2Q)/sqrt(p2E*p2d))
           # input(f"ROTATING PROJECTILE BY {.5*(p2E + p2d - p2Q)/sqrt(p2E*p2d)}")

            self.state.E -= Wk
            
            #particles.ELECTRONS.push_back(E - self.Wk)
            #particles.ELECTRONS.push_back(.5*(p2E + p2d - p2Q)/sqrt(p2E*p2d))
            
            
            
            
            
    
            
            #particles.ELECTRONS.push_back(Wk - Uk)
            #particles.ELECTRONS.push_back(Wk*Wk / ((Qs*Qs-1)/(Qs*Qs) / p2Q * (1 + (p2Q - Wk*Wk)/(2*Wk*(E + _2ELECTRON_REST_ENERGY)))**2))
            
        else:
            #input(f"TRANSVERSE DISTANT - NO DEFFLECTION")
            self.throwAZIMUTH()

            
            if Wk -  Uk > CUT_OFF:
                print("WARNING: Hard transverse distant interaction")
                print("Above threshold particle is being depoisted locally.")
                print("This is expected behaviour in the current version.")
                print("")
                (<V> self.state.current_region).depositLOCAL(self.state.pos,Wk - Uk)
                #el = Electron._new(self.state)
                #el.state.E = Wk-Uk
                #el.invert_dire()

                #el.rotateTHETA(-1)
                #self.secondary.append(el)
                #self.nSECONDARY += 1
            else:
                (<V> self.state.current_region).depositLOCAL(self.state.pos,Wk - Uk)

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

        (<V> self.state.current_region).depositLOCAL(self.state.pos,Etot)






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
            MID = START + (END - START)//2 #prevents overflow somehow 
            
            xMID = EAX[MID]
            
            if self.state.E is xMID: #found the value
                return MID
            
            if self.state.E < xMID: # discard right side
                END = MID - 1 # do not include mid
                continue
            
            START = MID + 1
        return END 
                    
            
