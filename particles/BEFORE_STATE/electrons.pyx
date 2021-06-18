# cython: profile=False
# cython: annotate=True
# distutils: language = c++ 
# distutils: extra_compile_args = -std=c++11

print(">>>>>   IMPORTING ELECTRONS")



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





# ---- RANDOM SAMPLER --- #######################



from libc.stdlib cimport rand, RAND_MAX, srand




from ..random.mixmax.interface cimport mixmax_engine


cdef extern from "math.h":
    float INFINITY

#from numpy.random import rand as urand
###################################################################

# CONSTANTS AND GLOBALS
cdef double CUTOFF = __photonCUTOFF__
cdef double photonCUTOFF = __photonCUTOFF__
cdef double MIN_CUT_OFF = min(__photonCUTOFF__, __electronCUTOFF__)

cdef double E0_el = db.E0_electron*1e-3


    
    
    
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
        return self.E
    
    
    @staticmethod
    cdef Electron _new(double E,
              double x, double y, double z,
              double eyx, double eyy, double eyz,
              double ezx, double ezy, double ezz, 
              Volume current_region):
        
        cdef Electron self
        self = <Electron>Electron.__new__(Electron)
        self.current_region = current_region
        self.x = x
        self.y = y
        self.z = z
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
    cdef Electron _newISOTROPIC(double E,
              double x, double y, double z,
              Volume current_region,
               mixmax_engine* genPTR):
        
        cdef Electron self
        self = <Electron>Electron.__new__(Electron)
        self.current_region = current_region
        self.x = x
        self.y = y
        self.z = z
        self.E = E
        
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
#        self.throwAZIMUTH()
        self.index = <int> (10*self.z)

        return self    
    
    
    
    
    # @staticmethod #custom constructor for speed
    # cdef Electron _newISOTROPIC(Volume current_region, double E, Vector pos):
        
    #     cdef Electron self
    #     self = <Electron>Electron.__new__(Electron)
    #     self.pos   = pos
    #     self.current_region = current_region
        
    #     axis = ez0
    #     ey = ey0.rotateAngle(axis, 2*pi*rand())
            
    #     axis = ey
    #     ez = ez0.rotateAngle(axis, pi*rand())        
        
        
    #     self.ey = ey
    #     self.ez = ez
    #     self.E = E
    #     self.s_max = 100
    #     return self







    ####################################################################################
    ########                           RUN                                      ########
    ########                           RUN                                      ########
    ####################################################################################

    cdef void _run(Electron self, mixmax_engine* genPTR):
        IF _DEBUG_BASIC: print("> ELECTRON")
        self.secondary = deque()

        self.nSECONDARY = 0
        if self.E < CUT_OFF:
            return
        
        self.genPTR = genPTR
        
        
        self.update_references()
        
        #s#elf.record()
        IF RECORD: self.record()
        
   #     if# self.z == INFINITY or self.z == -INFINITY: print("------------------------------------------- INTIALIZED LIKE THIS")
        cdef double r
        cdef double  tau, S_soft
        cdef bint delta= False
        cdef double tau2 
        self.current_region.p = self
        while True:
            #     raise RuntimeError("ELECTRON: imfp_max = 0")
            
            
            
            
            self.s = -log(1e-6 + (1 - 1e-6)*self.genPTR.get_next_float() )/self.imfp_max
         #   if self.s == INFINITY or self.s == -INFINITY: print(">>><<<<<<<<<< ifni")
            
            
            
            # if self.s == 0:
                
                
            #     raise RuntimeError("ELECTRON: self.s = 0")
            #if self.s - self.s_max > :
            
            if self.s > self.s_max:
                self.s = self.s_max
                delta = True
            
         
            self.sample_w()
 
            
            # if self.z == INFINITY or self.z == -INFINITY:
            #     self.Z.pop_back()
            #     self.X.pop_back()
            #     self.Y.pop_back()
                
            #     print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 11111111")
            #     print(SS)
            #     print("z", self.z)
            #     print("ez", self.ezx, self.ezy, self.ezz)
                
            #     print("tau", tau)
            #     print("s", self.s)
            #     print("tau2", tau2)
            #     print("imfp_max", self.imfp_max)
            #     print("E", self.E)
                
                
                
            #     print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            
            
            tau = self.s*self.genPTR.get_next_float()
            
            
            # self.x += self.ezx*tau
            # self.y += self.ezy*tau
            # self.z += self.ezz*tau
            
            
            # if self.z == INFINITY or self.z == -INFINITY:
            #     print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 2222222222222")
                
            #     print("z", self.z)
            #     print("ez", self.ezx, self.ezy, self.ezz)
                
            #     print("tau", tau)
            #     print("s", self.s)
            #     print("tau2", tau2)
            #     print("imfp_max", self.imfp_max)
            #     print("E", self.E)
                
                
                
            #     print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            
            #self.move(tau)
            #self.record()
            self.state.L = tau2
            if self.current_region.move(self.state):
                self.E -= self.w/self.s *tau
                
                if self.E < CUT_OFF:
                    return
                
                self.update_references()
                continue
            
            IF RECORD: self.record()
            
            
            S_soft = self.w/self.s
            self.E -= S_soft*tau
            
            if self.E < CUT_OFF:
                    return
                
            self.do_hinge()
            
            
            
            tau2 = self.s - tau
            
 
            
            
            # self.x += self.ezx*tau2
            # self.y += self.ezy*tau2
            # self.z += self.ezz*tau2
            
            
            # if self.z == INFINITY or self.z == -INFINITY:
                
            #     print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 3333333333333333")
                
            #     print("z", self.z)
            #     print("ez", self.ezx, self.ezy, self.ezz)
                
            #     print("tau", tau)
            #     print("s", self.s)
            #     print("tau2", tau2)
            #     print("imfp_max", self.imfp_max)
            #     print("E", self.E)
                
                
                
            #     print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            #self.move(tau2)
            #self.record()
            self.state.L = tau2
            if (<Volume> self.state.current_region).move(self.state):
                
            
                self.E -= S_soft*tau2
                if self.E < CUT_OFF:
                    return
                self.update_references()
                continue
            IF RECORD: self.record()
            
            
            self.E -= tau2*S_soft
            
            if self.E < CUT_OFF:
                return
            
            
            
            if delta:
               # self._delta()
                self.update_imfp()
                delta = False
                continue
            
            
            self.update_imfp_cumul()
            
            r = self.genPTR.get_next_float()*self.imfp_max
            
            
            
            if r < self.IMFP_CUMUL.C1: 
                self._inelastic()
                if self.E < CUT_OFF:
                    return
                self.update_imfp()
                
                
            elif r < self.IMFP_CUMUL.C2:
                self._brem()
                if self.E < CUT_OFF:
                    return
                self.update_imfp()
                
            elif r < self.IMFP_CUMUL.C3: 
                self._elastic()
                
           # else: self._delta() 
            
            
            
            
            
 
            

    










    ####################################################################################
    ########                           UPDATE                                   ########
    ########                           METHODS                                  ########
    ####################################################################################

    cdef void update_references(self) :
        """
        Updates all references. Called when there is a region crossing.
        """
        
        #getting material from current region
        self.current_material = self.current_region.material
        

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
        Updates inverse mean free paths. Constructs cumul.
        Called when there is a region crossing or energy of photon has changed.
        """

        cdef int i = self.find_index()  
 
        
        self.IMFP_CUMUL.C1 = self.inelastic.imfpA[i] + self.inelastic.imfpB[i]*self.E
        self.IMFP_CUMUL.C2 = self.IMFP_CUMUL.C1 + self.brem.imfpA[i] +    self.brem.imfpB[i]   *self.E
        self.IMFP_CUMUL.C3 = self.IMFP_CUMUL.C2 + self.elastic.imfpA[i] + self.elastic.imfpB[i] * self.E
        
        #self.electron.IMFP(self.E)
        self.s_max = 4/self.IMFP_CUMUL.C3
        cdef double Emin = self.E - self.electron.find_wmax(self.s_max, self.E)
        
        
        if Emin is 0:
            self.imfp_max = self.IMFP_CUMUL.C3
            return
        
        i = self.electron.find_index(Emin)
        
        
        self.imfp_max = fmax(self.IMFP_CUMUL.C3, 
                             self.electron.imfpA[i] + self.electron.imfpB[i] * Emin
                             )
        
        
        
    cdef inline void _delta(self):
        pass
        #imfp_delta = self.imfp_max - self.IMFP_CUMUL.C3
        
        
        
        self.move(  -log(1e-2 + (1 - 1e-2)*self.genPTR.get_next_float()  )/( self.imfp_max - self.IMFP_CUMUL.C3)  )
        # if self.z == INFINITY or self.z == -INFINITY:
        #     print(self.imfp_max, self.IMFP_CUMUL.C3)
        
        
        
    cdef inline void update_imfp(Electron self) :
        
        cdef double Emin = self.E - self.electron.find_wmax(self.s_max, self.E)
        cdef int i, j
        
        if Emin is 0:
            j = self.find_index()
            self.imfp_max = self.electron.imfpA[j] + self.electron.imfpB[j] * self.E
            self.s_max = 4/self.imfp_max
            return
        
        
        
        i = self.electron.find_index(Emin)
        j = self.find_index()
        #cdef double imfpH = 
   
        self.s_max = self.electron.imfpA[j] + self.electron.imfpB[j] * self.E
        self.imfp_max = fmax(self.s_max,
                             self.electron.imfpA[i] + self.electron.imfpB[i] * Emin
                             )

        self.s_max = 4/self.s_max

        

        
        

    ####################################################################################
    ########                           RANDOM                                   ########
    ########                           SAMPLING                                  ########
    ####################################################################################
    
    
    
    cdef inline void sample_w(self):
        
        
        #self.SP     = self.inelastic.softSP._eval(self.E) + self.brem.softSP._eval(self.E)
        #self.STRAGG = self.inelastic.softSTRAGG._eval(self.E) + self.brem.softSTRAGG._eval(self.E)

        #self.avgW   = self.s * (self.inelastic.softSP._eval(self.E)     + self.brem.softSP._eval(self.E)    )
        #self.varW   = self.s * (self.inelastic.softSTRAGG._eval(self.E) + self.brem.softSTRAGG._eval(self.E))
        
        cdef int i = self.find_index()

        self.avgW = self.s * ( self.electron.softSPA[i]     + self.electron.softSPB[i]*self.E)
        self.varW = self.s * ( self.electron.softSTRAGGA[i] + self.electron.softSTRAGGB[i]*self.E)

        cdef double avgW2 = pow(self.avgW, 2)
        
        if avgW2 > 9*self.varW:

            self.w = sqrt(self.varW)*self.current_material.electron.gauss._sample() + self.avgW #CONFIRM
            return
        
        if avgW2 > 3*self.varW:
            #                     6*3**.5
            self.w = self.avgW + (10.392304845413264 * sqrt(self.varW))*self.genPTR.get_next_float()
            return
        
        #cdef doubla a = (3*varW - avgW**2)/(3*varW + 3*avgW**2)
        #if self.genPTR.get_next_float() < (3*varW - avgW2)/(3*varW + 3*avgW2):
            
        #if 3*self.genPTR.get_next_float() < (3*self.varW - avgW2)/(self.varW + avgW2):
        if 3*self.genPTR.get_next_float()*(self.varW + avgW2) < (3*self.varW - avgW2):
            self.w = 0.
            return
        #self.w = self.genPTR.get_next_float()*(3*varW + 3*avgW**2)/2/avgW
        #self.w = 1.5*self.genPTR.get_next_float()*(self.varW + avgW2)/self.avgW
        self.w = 1.5*self.genPTR.get_next_float()*(self.varW + avgW2)/self.avgW
        return
    
    
    cdef inline void do_hinge(self):
        """
        Condensed history hinge.
        """
      # print("HINGE")
    
        #self.imfp0 = self.elastic.imfp0._eval(self.E)
        #cdef double rc = 1 - self.elastic.imfp._eval(self.E)/imfp0
        #cdef double T1, T2
       # self.T1, self.T2 = self.elastic.sampler.T(self.E,   1 - self.elastic.imfp._eval(self.E)/self.imfp0           )
        
        #print(rc, imfp0, T1, T2)
        
        #cdef double sIMFP1 = 2*T1/imfp0
        #cdef double sIMFP2 = 6*(T1 - T2)/imfp0
        
        
        
        cdef int i = self.find_index()
        cdef double sIMFP1, sIMFP2
        sIMFP1 = self.elastic.sIMFP1A[i] + self.elastic.sIMFP1B[i]*self.E
        sIMFP2 = self.elastic.sIMFP2A[i] + self.elastic.sIMFP2B[i]*self.E
        
        #REMINDER : sIMFP interpolations already include the minus sign for the exp()
        cdef double mu1 = .5*(1 - exp(self.s *sIMFP1))
        #cdef double mu2 = mu1 - (1 - exp(self.s*    6*(T1 - T2)/imfp0    ))/6
        
        
        cdef double den = (1 - 2*mu1)
        cdef double mu0 = (2*mu1 - 3*(mu1 - 0.16666666666666666*(1 - exp(self.s*sIMFP2))))/den
        
        cdef double cos
        if self.genPTR.get_next_float() < den + mu0: mu0 = self.genPTR.get_next_float()*mu0
        else:                   mu0 = mu0 + self.genPTR.get_next_float()* (1 - mu0)
        
        
        #cos = 1-2*mu0
        
        #if cos > 1: cos = 1
        
        #self.change_direction(cos , twoPI*self.genPTR.get_next_float())
        self.throwAZIMUTH()
        self.rotateTHETA(1-2*mu0)
        
    
    
    
    cdef inline void _elastic(self):
        cdef int i = self.find_index()
        #cdef DIST dist
   
        #LOGeax
        if self.genPTR.get_next_float()*(LOGeax[i+1] - LOGeax[i]) < ( LOGeax[i+1] - log(self.E) )  :
            
            #dist = self.elastic.DISTRIBUTIONS[i]
            self.mu = (<DIST>(self.elastic.DISTRIBUTIONS[i])).sample(self.genPTR)
        else:
            #dist = self.elastic.DISTRIBUTIONS[i + 1]
            #self.mu = dist.sample()
            self.mu = (<DIST>(self.elastic.DISTRIBUTIONS[i+1])).sample(self.genPTR)
        
        
        
        #self.mu = self.elastic.sample(self.E)
        
        #self.mu = (1-2*self.mu)
        

            
        self.throwAZIMUTH()
        self.rotateTHETA(1-2*self.mu)
        #self.change_direction(self.mu, twoPI*self.genPTR.get_next_float())
        
        
    cdef inline void _brem(self):
        IF _SIGNAL_INTERACTION: print("BREM")
        cdef double k, theta
        
        k, theta = self.brem.sampler.full_sample(self.E, self.genPTR) # <- dont forget secondary particle
        
        
        cdef double dE = k*self.E
        
        self.E -= dE
        
        
        
        if dE < photonCUTOFF:
            return
        
        
        # axis = ez0ss
        # ey = ey0.rotateAngle(axis, twoPI*self.genPTR.get_next_float())
        
        # axis = ey
        # ez = ez0.rotateAngle(axis, theta)
              
        cdef Photon p = Photon._new(dE, self.x, self.y, self.z,
                                      self.eyx, self.eyy,self.eyz,
                                      self.ezx, self.ezy,self.ezz, 
                                      self.current_region)
        p.genPTR = self.genPTR
        p.throwAZIMUTH()
        p.rotateTHETA(cos(theta))
        self.nSECONDARY += 1
        self.secondary.append(p)
        
        # self.secondary.append(
        #     Photon._new(
        #                self.current_region,
        #                dE, 
        #                self.pos,
        #                ey,
        #                ez,
        #                )
        #     )
        
        
        #self.change_direction(cos(theta), 2*pi*self.genPTR.get_next_float())
        
        #self.E -= k*self.E
        #self.update_imfp()
        
        
        
    # @cython.boundscheck(False)
    # @cython.wraparound(False) 
    # @cython.initializedcheck(False)
    # @cython.cdivision(True) 
    # cdef void sampleCLOSE(self, double E, mixmax_engine *genPTR, PARTICLES *particles):
        

    

        
        
    #     # EFFECT ON TARGET
    #     cdef double W = r*E
    #     cdef double cos2 = (E - W)*(E + _2ELECTRON_REST_ENERGY)/E/(E - W + _2ELECTRON_REST_ENERGY)
    #     particles.ELECTRONS.push_back(E - W)
    #     particles.ELECTRONS.push_back(sqrt(cos2))
        
    #     # SECONDARY ELECTRON
    #     cos2 = (W)*(E + _2ELECTRON_REST_ENERGY)/E/(W + _2ELECTRON_REST_ENERGY)
    #     particles.ELECTRONS.push_back(W)
    #     particles.ELECTRONS.push_back(sqrt(cos2))
        
        
        
    cdef inline void _inelastic(self):
        IF _SIGNAL_INTERACTION: print("INELASTIC")
        #cdef double [::1]
        
        #self.GOS.sample(self.genPTR, self.find_index(), self.E, &particles)
        cdef int j = self.find_index()
        cdef int i
        cdef double cumul = 0
        cdef double r = self.genPTR.get_next_float()*self.IMFP_CUMUL.C1
        for i in range(0, self.inelastic.lenposs, 6):
            cumul += self.inelastic.arr[j, i] + self.E* self.inelastic.arr[j, i + 1]
            if r < cumul:
                break
            
            
        cdef double Wk = self.inelastic.arr[j, i + 2]
        cdef double Uk = self.inelastic.arr[j, i + 3]
        
        #j = i % 3
        
        #### sample secondary
        cdef Electron el
        


        cdef int _
        cdef double p2E, p2d, Qs, p2Q, a, kc
        
        if Uk < self.E < Wk:
            if self.E < Uk: raise RuntimeError("this hsould not be happening")
            pass
            self.E -= Uk 
            # this is equivelent to assuming that the electron has given all its energy, an ejected electron with the same direction
            # 
            # if self.E - Uk > CUT_OFF:
            #     self.E -= Uk
            #     el = Electron._new(Wk-Uk, self.x, self.y, self.z,
            #                        self.eyx, self.eyy, self.eyz,
            #                        self.ezx, self.ezy, self.ezz,
            #                        self.current_region)
                

            #     self.secondary.append(el)
            #     self.nSECONDARY += 1
        
        elif i/6 % 3 == 0:
            kc = Wk/self.E
            a = (self.E/(self.E + ELECTRON_REST_ENERGY))**2
            p2E = 1+ 5*a*kc/2
            p2d = (1 - 2*kc)
            p2Q = p2d/(5*a*kc)
            
            for _ in range(100_000):
                r = p2E*self.genPTR.get_next_float()
                
                if r < 1: r = kc / (1 - r*p2d)
                else: r = kc + (r - 1)*p2Q
                
                if r < kc: continue
                if r > .5: continue
            
             #   else: P = (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(k*(1-k))))
                    
                if self.genPTR.get_next_float()*(1+5*a*r*r) < r*r*  (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(r*(1-r))  )):
                    break
            else: 
                print(">>>>>> rejection sampling took too mucch")
                import time
                print("kc", kc)
                print(p2E)
                print("a", a)
                print("E", self.E)
                print("Wk", Wk)
                print("Uk", Uk)
                time.sleep(1000)
                
                
                
                
            Wk = r*self.E
            
            
            
            self.throwAZIMUTH()
            
            if Wk - Uk > CUT_OFF:
                el = Electron._new(Wk-Uk, self.x, self.y, self.z,
                                   -self.eyx, -self.eyy, -self.eyz,
                                   self.ezx, self.ezy, self.ezz,
                                   self.current_region)
                
                Qs = (self.E/ELECTRON_REST_ENERGY + 1)
                
                el.rotateTHETA(sqrt(       (Wk)*(self.E + _2ELECTRON_REST_ENERGY)/self.E/(Wk + _2ELECTRON_REST_ENERGY)      )     )
                self.secondary.append(el)
                self.nSECONDARY += 1
        
            self.rotateTHETA(sqrt((self.E - Wk)*(self.E + _2ELECTRON_REST_ENERGY)/self.E/(self.E - Wk + _2ELECTRON_REST_ENERGY)))
            
            self.E -= Wk
            # sample close
        elif i/6 % 3 == 1: # sample L far
            #self.E -= Wk
            p2E = self.E * (self.E + _2ELECTRON_REST_ENERGY)
            p2d = (self.E - Wk) * ((self.E - Wk) + _2ELECTRON_REST_ENERGY)   #self.p2(E - self.Wk)
            
            
            Qs = sqrt( sqrt( p2E )  - p2d  + ELECTRON_REST_ENERGY*ELECTRON_REST_ENERGY  ) - ELECTRON_REST_ENERGY
            
            
            
            
            Qs = Qs/(1 + Qs / _2ELECTRON_REST_ENERGY)
    
             
       
            Qs = Qs / ((   Qs/Wk *(1 + Wk/_2ELECTRON_REST_ENERGY )      )**(self.genPTR.get_next_float()) - Qs/_2ELECTRON_REST_ENERGY)
    
    
            p2Q = Qs * (Qs + _2ELECTRON_REST_ENERGY)
    
            
            
   
            
            
            self.throwAZIMUTH()
            
            if Wk - Uk > CUT_OFF:
                
                el = Electron._new(Wk-Uk, self.x, self.y, self.z,
                                   -self.eyx, -self.eyy, -self.eyz,
                                   self.ezx, self.ezy, self.ezz,
                                   self.current_region)
                
                Qs = (self.E/ELECTRON_REST_ENERGY + 1)
                
                el.rotateTHETA(Wk*Wk / ((Qs*Qs-1)/(Qs*Qs) / p2Q * (1 + (p2Q - Wk*Wk)/(2*Wk*(self.E + _2ELECTRON_REST_ENERGY)))**2))
                self.secondary.append(el)
                self.nSECONDARY += 1
            
            self.rotateTHETA(.5*(p2E + p2d - p2Q)/sqrt(p2E*p2d))
            self.E -= Wk
            
            #particles.ELECTRONS.push_back(E - self.Wk)
            #particles.ELECTRONS.push_back(.5*(p2E + p2d - p2Q)/sqrt(p2E*p2d))
            
            
            
            
            
    
            
            #particles.ELECTRONS.push_back(Wk - Uk)
            #particles.ELECTRONS.push_back(Wk*Wk / ((Qs*Qs-1)/(Qs*Qs) / p2Q * (1 + (p2Q - Wk*Wk)/(2*Wk*(E + _2ELECTRON_REST_ENERGY)))**2))
            
        else:
            self.throwAZIMUTH()

            
            if Wk -  Uk > CUT_OFF:
                el = Electron._new(Wk-Uk, self.x, self.y, self.z,
                                   self.eyx, self.eyy, self.eyz,
                                   -self.ezx, -self.ezy, -self.ezz,
                                   self.current_region)
                #el.rotateTHETA(-1)
                self.secondary.append(el)
                self.nSECONDARY += 1
            #self.rotateTHETA(1)
                
                
            self.E -= Wk
                
            
            #sample transverse
            
        
 
       # self.E -= self.inelastic.arr[j, i + 2]
      #  el = Electron._newISOTROPIC(self.inelastic.arr[j, i + 2], self.x, self.y, self.z, self.current_region, self.genPTR)
      #  self.nSECONDARY += 1
      #  self.secondary.append(el)
        
      
       # print(Uk)
        if Uk < MIN_CUT_OFF: return 
        
        cdef PARTICLES particles
        (<_Atom> self.inelastic.arr_atoms[<int>self.inelastic.arr[j, i + 4]] ).run(<int>self.inelastic.arr[j, i + 5], 
                                                                       &particles,
                                                                       self.genPTR
                                                                       )


        
         
        cdef Photon ph
        for i in range(particles.PHOTONS.size()):
            E = particles.PHOTONS.back()
            particles.PHOTONS.pop_back()
            if E < photonCUTOFF: 
                continue
            
            ph = Photon._newISOTROPIC(E, self.x, self.y, self.z, self.current_region,  self.genPTR)
            self.nSECONDARY += 1
            self.secondary.append(ph)
        
        for i in range(particles.ELECTRONS.size()):
            E = particles.ELECTRONS.back()
            particles.ELECTRONS.pop_back()
            if E < CUT_OFF: 
                continue
 
            el = Electron._newISOTROPIC(E, self.x, self.y, self.z, self.current_region, self.genPTR)
            self.nSECONDARY += 1
            self.secondary.append(el)
            
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

        
        # self.E = particles.ELECTRONS[0]
        
        # if self.E == 0:

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
        #self.GOS.update(self.E, False)
        
        #self.E, self.cos, self.Esec, self.cos_sec = self.GOS.sample(self.current_region, self.x, self.y, self.z)
        
        # if self.cos>1:
        #     #print("cos > 1")
        #     self.cos = 1

        #self.throwAZIMUTH()
       # self.rotateTHETA(self.cos)
        #self.change_direction(self.cos, twoPI*self.genPTR.get_next_float())

        ## RUNNING pyRelax
        
        #atom.ionize(atom.SHELLS[I])
        
        
        
        
        
        
        
        #self.secondary.extend(self.GOS.secondary)
        #self.nSECONDARY += self.GOS.nSECONDARY


        
        #if self.Esec < CUT_OFF:
        #    return
        
        
        
        # if self.cos_sec > 1: 
        #     self.cos_sec = 1
            
            
        
        
        #cdef Electron p = Electron._new(self.Esec, self.x, self.y, self.z,
         #                             self.eyx, self.eyy,self.eyz,
         #                             self.ezx, self.ezy,self.ezz, 
          #                            self.current_region)
        #p.genPTR = self.genPTR
        #p.throwAZIMUTH()
        #p.rotateTHETA(self.cos_sec)
        #self.nSECONDARY += 1
        #self.secondary.append(p)
                
        
        
        
        
        
        # axis = ez0
        # ey = ey0.rotateAngle(axis, twoPI*self.genPTR.get_next_float())
        
        # axis = ey
        # ez = ez0.rotateCos(axis, self.cos_sec)  
            
            
            
            
            
            
            
            
            
        # self.secondary.append(
        #     Electron._new(
        #        self.current_region,
        #        self.Esec, 
        #        self.pos, 
        #        ey,
        #        ez, 100)
        #     )
        
        
        

        
 

        
    
    def getindex(self):
        return self.find_index()
    
    def setE(self, double E):
        self.E = E
    
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
                    
            
            

    # cdef int _anihilation(Electron self) except -1:
    #     cdef double gamma = E/0.511e6 + 1
    #     cdef double f_min = 1/(gamma + 1 + (gamma**2 - 1)**.5)
    #     cdef double v, r
    #     while 1:
    #         r = rand()
    #         v = f_min * ((1-f_min)/f_min)**r
    #         r = rand()
    #         if r*g(f_min, gamma) >= g(v, gamma):
    #             break
            
    #     v = min(v, 1-v)
    #     cdef double Eplus  = (1 - v)*(self.E + 2*0.511e6)
    #     cdef double Eminus = v*(self.E + 2*0.511e6)
        

        



    ####################################################################################
    ########                           PYTHON                                   ########
    ########                          INTERFACE                                 ########
    ####################################################################################
    
    # @staticmethod #thin wrapper for python acess
    # def new(Volume space, Volume current_region,
    #              E     = 6.,
    #              pos   = Vector(0., 0., 0.),
    #              theta = 0.,
    #              phi   = 0.,
    #              ex    = Vector(1., 0., 0.), 
    #              ey    = Vector(0., 1., 0.), 
    #              ez    = Vector(0., 0., 1.),
    #              simulate_secondary = False,
    #              s_max = 100):
        
        
    #     return Electron._new(space, current_region, E, pos, 
    #                          theta, phi, 
    #                          ex, ey, ez, 
    #                          simulate_secondary,
    #                          s_max)    
    

    
    


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
