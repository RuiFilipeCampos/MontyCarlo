# distutils: language = c++
print(">>>> IMPORTING main.pyx")




class MAP(dict):
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError


# #### PREPARING ENERGY AXIS AND GRID
from ...settings import __montecarlo__
__path__ = __montecarlo__/'materials'/'electron'


# PATH = __path__/'elastic'
# PATH = str(PATH)


# #mu = np.load(path + "/muGRID.npy")

# LEeax = np.load(PATH + "/LEeax.npy")
# HEeax = np.load(PATH + "/HEeax.npy")
# _eax =  np.append(LEeax, HEeax[1:])


# eax = _eax



from ..._init import eax
from ..._init  cimport LIMS
from ..._init  cimport EAX








from ..._random.interface cimport mixmax_engine

from collections import deque

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import cumtrapz


from .. import database as db


from ...tools.interpol1 cimport InvRationalInterpolation, LinLinInterpolation, FastCubicSpline, hLinLinInterpolation
#from libcpp.vector cimport vector
from numpy import logspace
from ...tools cimport search, CubicInverseTransform
from ...tools import  CubicInverseTransform as CIT
from numpy cimport ndarray
from libc.math cimport fmax, fmin, log





# frexp gets exponent of double number = y*2**(exp)
cdef extern from "<math.h>" nogil:
    double frexp(double x, int* exponent)


# cdef int get_exp(double x):
#     cdef int exp;
#     frexp(x, &exp);
#     return exp;


cdef double E



# hashed =  np.array([get_exp(E) for E in eax], dtype = int)
# indexes = np.arange(0, len(hashed))
# Imax = int(max(hashed))
# lims = []


# for i in range(Imax + 1): #every possible value of the hash, index = hash
#     selected = indexes[hashed == i]
#     n = len(selected)
#     if n == 0: # no values in this range
#         lims.append(np.array([0, 0, 0], dtype = int))
#         continue
#     if n == 1: # only one value in this range -> 1 interval -> use selected index 
#         lims.append(np.array([0, 0, 1], dtype = int))
#         continue
#     # at least two intervals
#     lims.append(np.array([selected[0], selected[-1], n], dtype = int))

# LIMS = np.array(lims)





# hashed =  np.array([get_exp(E) for E in eax], dtype = int)
# indexes = np.arange(0, len(hashed))
# Imax = int(max(hashed))



# lims = [np.array([0, 0, 0], dtype = int)]


# cdef int i 
# for i in range(Imax + 1): #every possible value of the hash, index = hash
#     selected = indexes[hashed == i]
#     n = len(selected)
#     if n == 0: # either out of bounds or no values in this range
#         #if no values in this range, interpolate using last interval
#         n_last = lims[-1][2]
#         if n_last == 0: #out of bounds
#             lims.append(np.array([0, 0, 0], dtype = int))
#             continue
#         i_last = lims[-1][1]
#         lims.append(np.array([i_last, i_last, 1], dtype = int))
#         continue

#     lims.append(np.array([selected[0], selected[-1] , n], dtype = int))

# LIMS = np.array(lims[1:], dtype = int)

def makeLinLin(x, y):
    m = np.diff(y)/np.diff(x)
    
    #y  = m*x - m*x[i] + y[i]
    #m*(x - x[i]) + y[i]
    return - m*x[:-1] + y[:-1], m

getLinLin = makeLinLin



#getDCS = lambda Z: db.EEDL(Z)[(9, 8, 0, 0, 9, 22)]

#from numpy import array





















from numba import njit
















cimport cython


import scipy
import numpy as np
np.seterr(all='raise')

cdef double const = 2*(1.015387)**2
cdef double C = 2*(1.015387)**2
from numpy import exp
from numba import njit

#@njit
def G1(w, avg = 0, std = 1):
    dw2 = (w - avg)**2
    std2 = std**2
    if not dw2 < 9*std2:
        return 0
    
    
    return exp(-dw2/ std2 / const   )

def G(w, avg = 0, std = 1):
    cdef double dw2 = (w - avg)**2
    
    if dw2 > 9*std**2:
        return 0
    cdef double std2 = 2*(1.015387*std)**2
    
    return exp(-dw2/ std2  )



def rebuildElectron(this):
    cdef Electron self
    self = <Electron> Electron.__new__(Electron)

    self.inelastic = this.inelastic
    self.brem =  this.brem

    self.softSP = this.softSP
    self.softSPA = this.softSPA
    self.softSPB = this.softSPB
    
    self.softSTRAGG = this.softSTRAGG
    self.softSTRAGGA = this.softSTRAGGA
    self.softSTRAGGB = this.softSTRAGGB


    self.elastic = this.elastic

    self.imfpA = this.imfpA
    self.imfpB = this.imfpB
    self.Itable = this.Itable


    self.integral = this.integral
    self.invI = this.invI
    self.gauss = this.gauss
    return self



@cython.auto_pickle(True)
cdef class Electron:
    def __reduce__(self):
        this = MAP()
        this.inelastic = self.inelastic
        this.brem =  self.brem

        from numpy import array 
        this.softSP = array(self.softSP)
        this.softSPA = array(self.softSPA)
        this.softSPB = array(self.softSPB)
        
        this.softSTRAGG = array(self.softSTRAGG)
        this.softSTRAGGA = array(self.softSTRAGGA)
        this.softSTRAGGB = array(self.softSTRAGGB)


        this.elastic = self.elastic

        this.imfpA = array(self.imfpA)
        this.imfpB = array(self.imfpB)
        this.Itable = array(self.Itable)


        this.integral = self.integral
        this.invI = self.invI
        this.gauss = self.gauss
        return rebuildElectron, (this,)




    def __init__(Electron self, object formula):
        #self.elastic = Elastic(formula, density)
        
        #print("    MAKING INELASTIC")
        formula.log.add_header("\t &emsp MAKING INELASTIC")
        self.inelastic = Inelastic(formula, formula.density)

        formula.log.add_header("\t &emsp MAKING BREM")
        self.brem = Brem(formula)
        
        
        cdef ndarray fullSP = self.brem.fullSP + self.inelastic.fullSP

        cdef ndarray fullSTRAGG = self.brem.fullSTRAGG + self.inelastic.fullSTRAGG


        #print("    MAKING ELASTIC")
        #SPcol = lambda x: self.inelastic.fullSP._eval(x)
        #SPrad = lambda x: self.brem.fullSP._eval(x)
        
        formula.log.add_header("\t &emsp MAKING ELASTIC")

        formula.log.add_header("MAKING CONDENSED HISTORY")  
        import numpy as np
        self.softSP = self.inelastic.softSP + self.brem.softSP
        self.softSPA, self.softSPB = getLinLin(eax, self.softSP)
        
        self.softSTRAGG = self.inelastic.softSTRAGG + self.brem.softSTRAGG
        self.softSTRAGGA, self.softSTRAGGB = getLinLin(eax, self.softSTRAGG)


        self.elastic = Elastic(formula, fullSP, -np.array(self.inelastic.sIMFP1), -1e-3*np.array(self.inelastic.sIMFP2))
        
        
        #eax = logspace(-3, 3, 2500)*1e6
        #imfp = [self.inelastic.imfp(E) + 
        #        self.brem.imfp(E) +
        #        self.elastic.imfp(E) for E in eax]
        
        #self.imfp = hLinLinInterpolation(eax, imfp)
        
        #self.imfp = self.inelastic.imfp + self.brem.imfp + self.elastic.imfp
        
        
        import numpy as np
        arr = np.array
        self.imfpA = arr(self.inelastic.imfpA) + arr(self.brem.imfpA) + arr(self.elastic.imfpA)
        self.imfpB = arr(self.inelastic.imfpB) + arr(self.brem.imfpB) + arr(self.elastic.imfpB)
        #self.imfpA, self.imfpB = getLinLin(eax, self.imfp)
        del arr
        
        

        


        fig = formula.log.new_plot()

        formula.log.add_to_plot(fig, eax*1e-6, fullSP*1e-6/formula.density, label = "FULL SP")

        formula.log.add_to_plot(fig, eax*1e-6, self.softSP*1e-6/formula.density, label = "Soft Stopping Power")
        formula.log.finish_plot(fig, title = "Soft SP's", xlabel = "E(MeV)", ylabel = "SP (MeV cm^2 /g)", logscale = True)
        
        fig = formula.log.new_plot()
        formula.log.add_to_plot(fig, eax, fullSTRAGG/formula.density, label = "Full straggling")
        formula.log.add_to_plot(fig, eax, self.softSTRAGG/formula.density, label = "Soft straggling")
        formula.log.finish_plot(fig, title = "STRAHH", xlabel = "E(eV)", ylabel = "SP (eV**2 cm^2 /g)", logscale = True)

        
        
        #Emin = np.log10(50)
        #Emax = np.log10(99e6)


 
        #eax = np.logspace(Emin, Emax, 1000)
        
        #eax = logspace(-3, 3, 2500)*1e6
        
        #SPrad = [SPrad(E) for E in eax]
        #SPcol = [SPcol(E) for E in eax]
        
        #SPrad, SPcol = map(np.array, (SPrad, SPcol) )
        import numpy as np


        integrand = 1/self.softSP
        integrand = CubicSpline(eax, integrand)
        integral = integrand.antiderivative()
        self.Itable = integral(eax)
        self.integral = FastCubicSpline(integral.x, integral(integral.x))
        self.invI = LinLinInterpolation(self.Itable, eax)

        fig = formula.log.new_plot()

        formula.log.add_to_plot(fig, eax, self.Itable, label = "Soft Stopping Power")
        #formula.log.add_to_plot(fig, eax, self.softSTRAGG, label = "Soft straggling")
        formula.log.finish_plot(fig, title = "Soft SP's", xlabel = "E(eV)", ylabel = "SP (eV/cm)", logscale = True)


        
        
        print(f"    > creating gauss (to do: move outside of class)")
        self.gauss = CIT.fromCallable(G, -3, 3, num = 1500)


        import numpy as np
        SAMPLES = deque(self.gauss._sample() for i in range(500_000))
        SAMPLES = np.array(SAMPLES)
        formula.log.add_attribute("gaussian avg ", np.mean(SAMPLES))
        formula.log.add_attribute("gaussian std ", np.std(SAMPLES))

        

        
        #Wcc = formula.Wcc
        
        #print(f"    FINISHING UP INELASTIC (introducing cut off | Wcc = {Wcc})")
        #self.inelastic.cut(Wcc)
        

        
    @property
    def ELASTIC(self):
        return self.elastic
    
    def plotSP(self):
        import matplotlib.pyplot as plt
        import numpy as np
        arr = np.array
        # coherent
        
        
        
     #   plt.plot(eax[:-1], totalCS[0] + totalCS[1]*eax[:-1])

        plt.plot(eax[:-1], arr(self.softSPA) + eax[:-1]*arr( self.softSPB))



        plt.legend()
        plt.xscale("log")
        plt.yscale("log")
        plt.show()
        
    def plot(self):
        import matplotlib.pyplot as plt
        import numpy as np
        arr = np.array
        # coherent
        
        
        
     #   plt.plot(eax[:-1], totalCS[0] + totalCS[1]*eax[:-1])

        plt.plot(eax[:-1], arr(self.inelastic.imfpA) + eax[:-1]*arr( self.inelastic.imfpB), label = "inel")
        plt.plot(eax[:-1], arr(self.elastic.imfpA) + eax[:-1]*arr( self.elastic.imfpB), label = "incoh")
        plt.plot(eax[:-1], arr(self.brem.imfpA) + eax[:-1]*arr( self.brem.imfpB), label = "pair")

        plt.legend()
        plt.xscale("log")
        plt.yscale("log")
        plt.show()
        
        
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    @cython.cdivision(True)
    cdef int find_index(self, double E):
        cdef int i;
        frexp(E, &i);
        
        
        #cdef int i = get_exp(E)
        
        # if LIMS[i, 2] == 0:
        #     raise RuntimeError(f">>>>>  OUT OF BOUNDS E = {E} eV")
        cdef int k = LIMS[i, 2]
        if k == 1:
            return LIMS[i, 0]
        
        if k == 2:
            i = LIMS[i, 0]
            if E <= EAX[i + 1]: return i
            return i + 1
        
        if k == 3:
            i = LIMS[i, 0]
            if E <= EAX[i + 1]: return i
            if E <= EAX[i + 2]: return i + 1
            return i + 2
        
        if k == 4:
            i = LIMS[i, 0]
            if E <= EAX[i + 1]: return i
            if E <= EAX[i + 2]: return i + 1
            if E <= EAX[i + 2]: return i + 2
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
            
            if E == xMID: #found the value
                return MID
            
            if E < xMID: # discard right side
                END = MID - 1 #do not include mid
                continue
            
            START = MID + 1
        return END 
                    
            
            
        
        
    @cython.boundscheck(False)
    @cython.initializedcheck(False)
    @cython.cdivision(True)
    cdef double find_wmax(self, double smax, double E0):
        cdef double CONST, e
        CONST =  - smax + self.integral._eval(E0)
        e = self.invI._eval(CONST)
        if e <= 0:
            print("in find_wmax:", f"wmax = {e} | E = {E0} | smax = {smax}")
            return 0

        return E0 - e

    @property
    def inel(self): return self.inelastic
    
    
    def getbrem(self):
        return self.brem
    
    def EXTRACT_ALL(self):
        dic = {}
        #inelastic
        dic["el"] = self.elastic.EXTRACT_ALL()
        dic["brem"] = self.brem.EXTRACT_ALL()
        
        return dic
        




def rebuildInelastic(this):
    cdef Inelastic self
    self = <Inelastic> Inelastic.__new__(Inelastic)
    self.gosMOLECULE = this.gosMOLECULE

    self.fullSP     = this.fullSP     #FULL[:, 1]
    self.fullSTRAGG = this.fullSTRAGG #FULL[:, 2]
    self.softSP     = this.softSP     #SOFT[:, 1]
    self.softSTRAGG = this.softSTRAGG #SOFT[:, 2]  
    self.imfpA      = this.imfpA      #interp[0]
    self.imfpB      = this.imfpB      #interp[1]

    try:
        self.arr = this.arr
        self.lenposs = this.lenposs
        self.arr_atoms = this.arr_atoms
    except AttributeError:
        pass

    return self


from . import GOSmake_softIMFP
from . cimport GOSfinal
from . cimport GOS
from . import makeGOS
@cython.auto_pickle(True)
cdef class Inelastic:
    def __reduce__(self):
        this = MAP()

        this.gosMOLECULE = self.gosMOLECULE

        from numpy import array
        this.fullSP     = array(self.fullSP    ) #FULL[:, 1]
        this.fullSTRAGG = array(self.fullSTRAGG) #FULL[:, 2]
        this.softSP     = array(self.softSP    ) #SOFT[:, 1]
        this.softSTRAGG = array(self.softSTRAGG) #SOFT[:, 2]  
        this.imfpA      = array(self.imfpA     ) #interp[0]
        this.imfpB      = array(self.imfpB     ) #interp[1]

        try:
            this.arr = array(self.arr)
            this.lenposs = self.lenposs
            this.arr_atoms = array(self.arr_atoms)
        except AttributeError:
            pass

        return rebuildInelastic, (this,)


    def __init__(self, formula, double density):
        
        
        print("        DATA PROCESSING IN PYTHON")
        molecule, cb, delta = makeGOS.pyGOS(formula, density)
        print(molecule)
        print("        GENERATING GOS MODEL IN CYTHON")
        self.GOSmodel = GOS.CMolecule(molecule, cb, delta)
        
        print("        GENERATING FULL SP")
        # these values are not supposed to change after renormalization
        # I'll be collecting data here, and use that to renormalize
        FULL = [np.array(self.get(E, True)) for E in eax]
        FULL = np.array(FULL)
        
        self.gosMOLECULE = GOSfinal.gosMolecule(self.GOSmodel, formula)
        
        
        
        
        
       # interp = np.array(self.gosMOLECULE.totalCS)
       # interp = interp*formula.N
       # self.imfpA = interp[0]
      #  self.imfpB = interp[1]

        
        
        
        print("        CUTTING SHELLS FROM MODEL")
        formula.log.add_attribute("Wcc = ", formula.Wcc)
        #cutting the GOS model should be made in the new scheme
        self.GOSmodel.cut(formula.Wcc)
        shells_for_softIMFP = []
        empty = self.gosMOLECULE.cut(formula.Wcc, shells_for_softIMFP)
        
        print("        GENERATING hard sp")

        HARD = [np.array(self.get(E, False)) for E in eax]
        HARD = np.array(HARD)
        
        print("        GENERATING SOFT sp")

        SOFT = FULL - HARD
        
        
        self.fullSP = FULL[:, 1]
        self.fullSTRAGG = FULL[:, 2]

        self.softSP = SOFT[:, 1]
        self.softSTRAGG = SOFT[:, 2]        
        #self.imfp = HARD[:, 0]
        #self.imfpA, self.imfpB = makeLinLin(eax, self.imfp)
        #NUMBER DENSITY IS MISSING
        
        print(self.gosMOLECULE.totalCS)
        
        
        interp = np.array(self.gosMOLECULE.totalCS)
        interp = interp*formula.N
        self.imfpA = interp[0]
        self.imfpB = interp[1]
        #self.imfpA, self.imfpB = interp[0], interp[1]
        
        if len(shells_for_softIMFP) != 0:
            
            self.sIMFP1, self.sIMFP2 = GOSmake_softIMFP.do(formula, shells_for_softIMFP)
        else:
            self.sIMFP1 = np.zeros(len(eax))
            self.sIMFP2 = np.zeros(len(eax))
        if not empty:
            self.arr, self.lenposs, self.arr_atoms = self.gosMOLECULE.makearray()



        
        
        
        #self.fullSP = hLinLinInterpolation(eax, FULL[:, 1])
        #self.softSP = hLinLinInterpolation(eax, SOFT[:, 1])
        
        
        
        
        
        
        #self.softSTRAGG = hLinLinInterpolation(eax, SOFT[:, 2])

        #self.softSTRAGGA, self.softSTRAGGB = makeLinLin(eax, SOFT[:, 2])
  
        #self.imfp = hLinLinInterpolation(eax, HARD[:, 0])

        
        
        
    @property
    def GOS(self): return self.gosMOLECULE
    
    @property
    def GOSmodel(self): return self.GOSmodel
        
    cpdef (double, double, double) get(Inelastic self, double E, bint record):
        self.GOSmodel.update(E, record)
        return self.GOSmodel.imfp, self.GOSmodel.SP, self.GOSmodel.stragg
    
    cpdef double SP(Inelastic self, double E, bint record):
        self.GOSmodel.update(E, record)
        return self.GOSmodel.SP
    
    # cpdef void cut(self, double Wcc):
    #     soft_shells = self.GOSmodel.cut(Wcc)
        
    #     eax = logspace(-3, 3, 2500)*1e6
        
    #     hardSP = [self.SP(e) for e in eax]
    #     hardSP = np.array(hardSP)
        
    #     fullSP = [self.fullSP(e) for e in eax]
    #     fullSP = np.array(fullSP)

    #     self.softSP = LinLinInterpolation(eax, fullSP - hardSP)        








def rebuildBrem(this):
    cdef Brem self
    self = <Brem> Brem.__new__(Brem)
    self.sampler    = this.sampler    # BREM.sampler(molecule)
    self.softSP     = this.softSP     # y*1e6
    self.imfp       = this.imfp       # y
    self.imfpA      = this.imfpA      #
    self.imfpB      = this.imfpB      # getLinLin(x*1e6, y)
    self.softSTRAGG = this.softSTRAGG # y*1e6*1e6 # <--------- ------------------- confirm this, probly correct, just in case
    self.fullSP     = this.fullSP     # y*1e6
    self.fullSTRAGG = this.fullSTRAGG # y*1e6*1e6
    return self








from . import makeBrem
from . cimport BREM
@cython.auto_pickle(True)
cdef class Brem:

    def __reduce__(self):
        this = MAP()
        from numpy import array
        this.sampler    = self.sampler    # BREM.sampler(molecule)
        this.softSP     = array(self.softSP    ) # y*1e6
        this.imfp       = array(self.imfp      ) # y
        this.imfpA      = array(self.imfpA     ) #
        this.imfpB      = array(self.imfpB     ) # getLinLin(x*1e6, y)
        this.softSTRAGG = array(self.softSTRAGG) # y*1e6*1e6 # <--------- ------------------- confirm this, probly correct, just in case
        this.fullSP     = array(self.fullSP    ) # y*1e6
        this.fullSTRAGG = array(self.fullSTRAGG) # y*1e6*1e6
        return rebuildBrem, (this, )




    def __init__(self, formula):
        print("        MAKING BREM SAMPLER")
        molecule = makeBrem.makeX(formula)
        molecule = makeBrem.X(molecule)
        molecule.Wcr = formula.Wcr
        self.sampler = BREM.sampler(molecule)
        
        
        #formula or molecule?
        
        
        
        #eax = logspace(-3, 3, 2500)
        Wcr = formula.Wcr*1e-6
        
        
        print(f"        MAKING SOFT SP | Wcr = {Wcr} MeV")
        x, y = makeBrem.getSP(molecule, Wcr, eax = eax*1e-6)
        self.softSP = y*1e6
        #self.softSP = hLinLinInterpolation(x*1e6, y*1e6)
        
        print("        MAKING HARD IMFP")
        x, y = makeBrem.getIMFP(molecule, Wcr, eax = eax*1e-6)
        #self.imfp = hLinLinInterpolation(x*1e6, y)
        self.imfp = y
        self.imfpA, self.imfpB = getLinLin(x*1e6, y)
        
        print("        MAKING SOFT STRAGG")
        x, y = makeBrem.getSTRAGG(molecule, Wcr, eax = eax*1e-6)
        self.softSTRAGG = y*1e6*1e6 # <--------- ------------------- confirm this, probly correct, just in case
        #self.softSTRAGG = hLinLinInterpolation(x*1e6, y*1e6*1e6)
        
        print("        MAKING FULL SP")
        x, y = makeBrem.getFULLSP(molecule, eax = eax*1e-6)
        self.fullSP = y*1e6

        print("        MAKING FULL STRAGGA")
        x, y = makeBrem.getFULLSTRAGG(molecule, eax = eax*1e-6)
        self.fullSTRAGG = y*1e6*1e6






        #self.fullSP = hLinLinInterpolation(x*1e6, y*1e6)
        
        

    cpdef (double, double, double) get(self, double E):
        cdef int i = self.imfp.getINDEX(E)
        
        return self.imfp.evalbyINDEX(i, E), self.SP.evalbyINDEX(i, E), self.stragg.evalbyINDEX(i, E)
    
    def _get(self, double E):
        return self.get(E)
    
    def EXTRACT_ALL(self):
        dic = {}
        
        
        dic["SP"] = self.SP
        dic["imfp"] = self.imfp,
        dic["stragg"] = self.stragg,
        dic["fullSP"] = self.fullSP,
        dic["get"] = self.get
        
        return dic



#equations for determination of A0
#@njit
def eqn(A0, avg):
    a = 1 + A0
    return A0*(a*np.log(a/A0) - 1) - avg


#@njit
def eqn1(A0, avg):
    a = 1 + A0
    return   np.log(a/A0)*( 2*A0 + 1) - 2



#@njit
def eqn2(A0, avg):
    a = A0 + 1
    return (2*A0*a*np.log(a/A0) - 2*A0 - 1)/a/A0

#equations for determination of A
#@njit
def _mu(A):
    a = 1 + A
    return A*(a*np.log(a/A) - 1)

#@njit
def _mu2(A):
    return A*(1 - 2*_mu(A))

#@njit
def _eqn(A, const1, const2, const3):
    return const1*_mu(A) + const2*_mu2(A) + const3



@cython.auto_pickle(True)
cdef class DIST:
    def __init__(self):
        pass
    
    cdef double sample(self, mixmax_engine *genPTR):
        raise RuntimeError(".sample method called from DIST")
    
@cython.auto_pickle(True)
cdef class dist1(DIST):
    cdef double A, B, avgMU, r0, mu0, T10, T20, dr

    
    def __init__(self, A, B, avgMU, rc):
        self.A = A
        self.B = B
        self.avgMU = avgMU
        self.r0 = (1 - B)*(1 + A)*avgMU/(A + avgMU)
        self.mu0 = self.invCUMUL(self.r0)
        self.T10 = (1-self.B)*self.I1(self.mu0)
        self.T20 = (1-self.B)*self.I2(self.mu0)
        self.rc = rc
        self.dr = 1 - self.rc
        self.mu_c  = self.invCUMUL(self.rc)
        
        self.T1, self.T2 = self.T(self.mu_c, self.rc)

        
        
        
    cdef double sample(self, mixmax_engine *genPTR):
        return self.invCUMUL(self.rc + genPTR.get_next_float()*self.dr)
        
    cdef double invCUMUL(self, double r):
        if r < self.r0:
            return r*self.A/((1-self.B)*(1+self.A) - r)
        if r < self.r0 + self.B:
            return self.avgMU
        return (r - self.B)*self.A/((1-self.B)*(1+self.A) - r + self.B)
    
    
    cdef double I1(self, double mu):
        return self.A*((1+self.A)*log((self.A + mu)/self.A) - (1+self.A)*mu/(self.A + mu)    )
    
    cdef double I2(self, double mu):
        return self.A*((1+self.A)*mu**2/(self.A + mu) - 2*self.I1(mu))
    
    cdef (double, double) T(self, double mu_c, double rc):
        cdef double c = 1 - self.B
        
        if rc < self.r0:
            return (c*self.I1(mu_c), c*self.I2(mu_c))
        cdef double dr
        if rc < self.r0 + self.B:
            dr = rc - self.r0
            return (self.T10 + dr*self.mu0, self.T20 + dr*self.mu0**2)
        
        return c*self.I1(mu_c) + self.B*self.mu0, c*self.I2(mu_c) + self.B*self.mu0**2
        
 
        
@cython.auto_pickle(True)
cdef class dist2(DIST):
    
    cdef double A, B
    
    def __init__(self, double A, double B, double rc):
        self.A = A
        self.B = B
        self.rc = rc
        
        self.mu_c = scipy.optimize.bisect(self.CUMUL, 0, 1)
        self.T1, self.T2 = self.T(self.mu_c)
        
        if self.T1 == 0 or self.T2 == 0:
            print("")
            print("T1", self.T1)
            print("T2", self.T2)
            raise ValueError("dist1, one of T values is zero")
        
        
        
    cdef double CUMULeqn(self, double mu):
        """
        Cumulutative function of the second case distribution.
        Used to determine the mu_c values from the rc values.
        """
        
        if mu <.5:
            return (1 - self.B)*(1 + self.A)*mu/(self.A + mu) - self.rc
        
        return (1 - self.B)*(1 + self.A)*mu/(self.A + mu) + 4*self.B*(mu**2 - mu - .25) - self.rc
    
    cdef double sample(self, mixmax_engine *genPTR):
        """
        Sample from the distribution in the restricted domain (rc, 1).
        Uses composition rule since 0 < B < 1.
        """
        if genPTR.get_next_float() < self.B:
            return self.invCUMUL2(self.rc + genPTR.get_next_float()*(1 - self.rc))
            
        return self.invCUMUL1(self.rc + genPTR.get_next_float()*(1 - self.rc))

    cdef double invCUMUL1(self, double r):
        return self.A*r/(self.A - r + 1)
    
    cdef double invCUMUL2(self, double r):
        return .5*(1 + r**.5)
    
    cdef (double, double) T(self, double mu_c):
        cdef double logA = log(self.A)
        cdef double logAmu = log(self.A + mu_c)
        
        cdef double factor = self.B*(self.A**2 + self.A)
        
        cdef double I1 = (self.A/(self.A + mu_c) + logAmu) - logA + 1
        I1 = I1 * factor
        
        cdef double I2 =  self.A*(2*logA - 1) - self.A**2 /(self.A + mu_c) - 2* self.A*logAmu + mu_c
        I2 = I2 * factor
                          
        
        
        if mu_c < .5:
            return I1, I2
        
        I1 += (2.6666666666666666666666666666666667*mu_c - 2)*mu_c**2*(1-self.B)
        I2 += mu_c**3*(1 - self.B) * (2*mu_c - 1.3333333333333333333333333333333)
        
        return I1, I2
    


def rebuildElastic(this):
    cdef Elastic self
    self = <Elastic> Elastic.__new__(Elastic)
    self.imfp          =  this.imfp          
    self.imfpA         =  this.imfpA        
    self.imfpB         =  this.imfpB         
    self.DISTRIBUTIONS =  this.DISTRIBUTIONS 
    self.imfp          =  this.imfp          
    self.sIMFP1A       =  this.sIMFP1A      
    self.sIMFP1B       =  this.sIMFP1B       
    self.sIMFP2A       =  this.sIMFP2A      
    self.sIMFP2B       =  this.sIMFP2B      
    return self




#from . cimport makeElastic
#from . import makeElastic as mkEl
#from . import makeElastic
cdef class Elastic:


    def __reduce__(self):
        this = MAP()

        from numpy import array

        this.imfp          =  array(self.imfp         ) 
        this.imfpA         =  array(self.imfpA        )
        this.imfpB         =  array(self.imfpB        ) 
        this.DISTRIBUTIONS =  array(self.DISTRIBUTIONS) 
        this.imfp          =  array(self.imfp         ) 
        this.sIMFP1A       =  array(self.sIMFP1A      )
        this.sIMFP1B       =  array(self.sIMFP1B      ) 
        this.sIMFP2A       =  array(self.sIMFP2A      )
        this.sIMFP2B       =  array(self.sIMFP2B      )
        return rebuildElastic, (this,)


    def __init__(self, formula, fullSP, inel_sIMFP1, inel_sIMFP2):
        formula.log.add_header("\t MAKING ELASTIC")
        
        print("        MAKING ELASTIC SAMPLER")
        C1, C2 = formula.C1, formula.C2
        
        path = __path__/'elastic'
        path = str(path)
        self.path = path
        
        mu = np.load(path + "/muGRID.npy")
        
        #LEeax = np.load(path + "/LEeax.npy")
        #HEeax = np.load(path + "/HEeax.npy")
        #eax =  np.append(LEeax, HEeax)
        
        
        allDCS, SIGMA = self.compose(formula)
        SIGMA = SIGMA


        fig = formula.log.new_plot()

        formula.log.add_to_plot(fig, mu, allDCS[0], label = "E = " + str(eax[0]) + "eV")
        formula.log.add_to_plot(fig, mu, allDCS[100], label =  "E = " + str(eax[100]) + "eV" )
        formula.log.add_to_plot(fig, mu, allDCS[199], label =  "E = " + str(eax[199]) + "eV")

        formula.log.finish_plot(fig, title = "ELECTRON.ELASTIC EXAMPLES OF DCS", xlabel = "mu = (1-cos(theta))/2", ylabel = "DCS(cm**2)", logscale = True)

        fig = formula.log.new_plot()

        formula.log.add_to_plot(fig, eax, SIGMA[0], label = "CS0")
        formula.log.add_to_plot(fig, eax, SIGMA[1], label = "CS1")
        formula.log.add_to_plot(fig, eax, SIGMA[2], label = "CS2")

        formula.log.finish_plot(fig, title = "ELECTRON: Elastic CS (from ELESPA)", xlabel = "Energy(eV)", ylabel = "CS(cm^2)", logscale = True)


        
        # CONSTRUCTING MW MODEL
        cdef int i 
        i = np.searchsorted(eax, 100e6, side = "right") + 1 #index of element 100e6
        # print("INDEX", i)
        # print(eax[i])
        # print(eax[200])
        sigma0 = SIGMA[0, 200:]
        sigma1 = SIGMA[1, 200:]
        sigma2 = SIGMA[2, 200:]
        
        avgMU  = .5*sigma1/sigma0
        avgMU2 = avgMU - sigma2/sigma0/6.
        
        # CONSTRUCT A0 TABLE
        

        
        
        A0 = deque()
        cdef double avg
        import scipy
        for avg in avgMU:


            try:
                a0 = scipy.optimize.bisect(eqn, 1e-24, 1, args = (avg,))
                #print(f"A0 = {a0}, <mu> = {avg}")
                A0.append(a0)
                #A0.append(   scipy.optimize.newton(eqn, 0.001, 
                               #   fprime = eqn1, 
                               #   fprime2 = eqn2,
                               #   args = (avg,))    )
            except RuntimeError:
                import matplotlib.pyplot as plt
                _x = np.arange(-10, 10, .0001)
                __eqn = np.vectorize(lambda x: eqn(x, avg))
                plt.plot(_x, __eqn(_x))
                plt.plot(_x, [0]*len(_x))
                plt.show()
                raise RuntimeError("did not converge")
                
        A0 = np.array(A0)
        MWavgMU  = avgMU
        MWavgMU2 = A0*(1-2*MWavgMU)
        
        
        # CONSTRUCT TABLE OF IMFP0, IMFP1, IMFP2
        SIGMA0 = SIGMA[0, :]
        SIGMA1 = SIGMA[1, :]
        SIGMA2 = SIGMA[2, :]
        
        
        imfp0 = SIGMA0*formula.N
        imfp1 = SIGMA1*formula.N
        imfp2 = SIGMA2*formula.N
        
        # CONSTRUCT TABLE OF TOTAL SP
        #SPtot = []
        #cdef double E
        #for E in eax:
        #    SPtot.append(SPcol(E) + SPrad(E))
            
        #SPtot = fullSP
        
        # CONSTRUCT TABLE OF IMFP
        mfp = deque()
        
        vmax = np.vectorize(max)
        vmin = np.vectorize(min)
        
        mfp = vmax(1/imfp0, vmin(C1/imfp1,  C2*eax/fullSP))
        imfp = 1/mfp
        self.imfp = imfp
        self.imfpA, self.imfpB = getLinLin(eax, imfp)



        fig = formula.log.new_plot()

        formula.log.add_to_plot(fig, eax, 1/imfp0 * formula.density, label = "MFP")
        formula.log.add_to_plot(fig, eax, mfp * formula.density, label = "HARD MFP")
        formula.log.add_to_plot(fig, eax, C1/imfp1 *formula.density, label = "C1*MFP1")
        formula.log.add_to_plot(fig, eax, C2*eax/fullSP * formula.density, label = "C2*eax/fullSP")
       # formula.log.add_to_plot(fig, eax, SIGMA[2], label = "CS2")

        formula.log.finish_plot(fig, title = "ELECTRON.ELASTIC MFP (g/cm**2)", xlabel = "Energy(eV)", ylabel = "MFP", logscale = True)
        formula.log.add_paragraph("NOTE: mfp = max(mfp0, min(C1*mfp1,  C2*E/fullSP))")


        

        
        # CONSTRUCT rc, T1 and T2 tables
        rc = 1 - imfp/imfp0

        fig = formula.log.new_plot()
        formula.log.add_to_plot(fig, eax, rc, label = "rc")
        formula.log.finish_plot(fig, title = "ELECTRON.ELASTIC rc", xlabel = "Energy(eV)", ylabel = "rc", logscale = True)





        
        if np.any(rc < 0):
            INDEXES = np.arange(0, len(rc), 1)
            INDEXES = INDEXES[rc < 0]
            
            print(">>>> found the following negative values in rc")
            for ii in INDEXES:
                print(rc[ii])
                rc[ii] = 0.
                
            print(">>> I have assumed that their origin is numerical error and corrected them to 0")
            print(">>> please verify that they are small (~ 1e-16) ")
                    
            
        if np.any(rc > 1):
            raise RuntimeError(">>> VALUES GREATER THAN 1 in rc")
        #cdef int i 
        
        DISTRIBUTIONS = [] #since im already iterating over all dcs, might as well construct the distributions
        
        T1, T2 = deque(), deque()
        cdef DIST dist
        from scipy.integrate import trapz

        formula.log.add_header("ITERATING LOW E DCS'S")
        for i, dcs in enumerate(allDCS):
            prob = dcs/trapz(dcs, mu)#SIGMA0[i]
            
            cumul = cumtrapz(prob, mu, initial = 0)
            cumul, _mu = self.remove_duplicates( cumul, mu)
            dist = sFastCubicSpline(cumul, _mu, rc[i])

            DISTRIBUTIONS.append(dist)
            
            mu_c = dist.mu_c
            t1 = cumtrapz(prob*mu,    mu, initial = 0)
            T1.append( np.interp(mu_c, mu, t1) )
            
            t2 = cumtrapz(prob*mu**2, mu, initial = 0)
            T2.append( np.interp(mu_c, mu, t2) )
            
        cdef int j
        cdef double b, const1, const2, const3
        CASE = deque()
        A, B = deque(), deque()
        cdef double _mwavgMU2, _avgMU2, _avgMU
        cdef double __mu
        for j in range(i+1, len(imfp0)):
            k = j - 200
            _mwavgMU2 = MWavgMU2[k]
            _avgMU2 = avgMU2[k]
            _avgMU = avgMU[k]
            
            if _mwavgMU2 > _avgMU2: #CASE 1
                
                # MAKING DISTRIBUTIONS
                dist = dist1(A0[k],  (_mwavgMU2 - _avgMU2)/(_mwavgMU2 - _avgMU**2), _avgMU, rc[j])
                
                
                # CALCULATING T1 and T2
                #mu_c = dist.invCUMUL(rc[j])
                
                #t1, t2 = dist.T(mu_c, rc[j])
                
                T1.append(dist.T1)
                T2.append(dist.T2)
                DISTRIBUTIONS.append(dist)
                
                
            else: #CASE 2

                print(17/24, "IF THIS IS ZERO SEE LINE 1019")
                # CALCULATING COEFFICIENTS
                const1 = 17/24 - _avgMU2
                const2 = -(5/6 - _avgMU)
                const3 =  (-17/24*_avgMU + 5/6*_avgMU2)
                try:
                    print(j, k)
                    a = scipy.optimize.bisect(_eqn, 1e-60, A0[k], args=(const1, const2, const3))
                    print(a)
                except ValueError:
                   print(const1, const2, const3, A0[k])
                   raise ValueError("see previou traceback")
                
                __mu = _mu(a)

                dist = dist2(a, (_avgMU - __mu)/(5/6 - __mu), rc[j])
                
                T1.append(dist.T1)
                T2.append(dist.T2)
                DISTRIBUTIONS.append(dist)
                
        self.DISTRIBUTIONS = np.array(DISTRIBUTIONS)
        T1 = np.array(T1)
        T2 = np.array(T2)
                
            
        
        formula.log.add_paragraph("T1/T2")
        formula.log.log_table(T1, T2)
        # finally CONSTRUCT MFP1 and MFP2
        sIMFP1 = -2*T1*imfp0
        sIMFP2 = -6*(T1 - T2)*imfp0
        self.imfp = imfp0

        formula.log.add_paragraph("imfp0/imfp1")
        formula.log.log_table(imfp0, imfp1)
  
        
        self.sIMFP1A, self.sIMFP1B = getLinLin(eax, (sIMFP1 + inel_sIMFP1))
        self.sIMFP2A, self.sIMFP2B = getLinLin(eax, (sIMFP2 + inel_sIMFP2))



        from numpy import array as arr

        fig = formula.log.new_plot()


        x = eax[:-1]
        y = arr(self.sIMFP1A) + eax[:-1]*arr( self.sIMFP1B)
        formula.log.add_to_plot(fig, x, -y, label = "sIMFP1")
        formula.log.log_table(x, y)



        y = arr(self.sIMFP2A) + eax[:-1]*arr( self.sIMFP2B)
        formula.log.log_table(x, y)

        formula.log.add_to_plot(fig, x, -y, label = "sIMFP2")

        formula.log.finish_plot(fig, xlabel = "Energy (eV)", ylabel = "sIMFP", logscale = True)






        

        
    #cdef double sample(self, double E):
    @property
    def IMFP0(self): return self.imfp
    def plotSOFT(self):
        import matplotlib.pyplot as plt
        import numpy as np
        arr = np.array
        # coherent
        
        
        
     #   plt.plot(eax[:-1], totalCS[0] + totalCS[1]*eax[:-1])

        plt.plot(eax[:-1], arr(self.sIMFP1A) + eax[:-1]*arr( self.sIMFP1B), label = "sIMFP1")
        plt.plot(eax[:-1], arr(self.sIMFP2A) + eax[:-1]*arr( self.sIMFP2B), label = "sIMFP2")

        plt.legend()
        #plt.xscale("log")
        #plt.yscale("log")
        plt.show()
        
        
    def compose(self, object formula):
        formula.log.add_paragraph("> COMPOSING ELASTIC")
        #cdef ndarray allDCS, SIGMA0, SIGMA1, SIGMA2
        
        allDCS = np.zeros((200, 606))
        #SIGMA0 = zeros(200)
        #SIGMA1 = zeros(200)
        #SIGMA2 = zeros(200)
        
        shape = np.append(np.load(self.path + f"/{1}/HEtransportTCS.npy"),
                             np.load(self.path + f"/{1}/LEtransportTCS.npy")[:, 1:], axis = 1)
        
        
        print(shape.shape)
        print("SUBSTITUTE THIS VALUE ")
        
        SIGMA = np.zeros(shape.shape)
        
        cdef int Z
        cdef double x
        
        for Z, x in formula.items():
            formula.log.add_paragraph(f"Z = {Z}, x = {x}")
            dcs = np.load(self.path + f"/{Z}/DCS.npy")
            
            sigma  = np.append(np.load(self.path + f"/{Z}/LEtransportTCS.npy"),
                             np.load(self.path + f"/{Z}/HEtransportTCS.npy")[:, 1:], axis = 1)
            
            
            #dcs, sigma0, sigma1, sigma2, eax = self.getData(Z)
            
            SIGMA += x*sigma
            
            #SIGMA0 += x*sigma0
            #SIGMA1 += x*sigma1
            #SIGMA2 += x*sigma2
            
            allDCS += x*dcs
            
        return allDCS, SIGMA

    cdef object remove_duplicates(self, ndarray x, ndarray Y):
        cdef ndarray u, c, dup
        u, c = np.unique(x, return_counts=True)
        dup = u[c > 1]
        
        cdef bint keep = True
        cdef list new_y = []
        cdef int i
        cdef double y
        
        for i, y in enumerate(Y):
            
            if x[i] in dup:
                if keep:
                    new_y.append(y)
                    keep = False
                    continue
                else: continue
            
            if keep is False:
                keep = True
            new_y.append(y)
        
        return u, np.array(new_y)        



        
    cdef double sMFP(self, double E):
        pass
        
    
    cpdef double getSample(self, double E):
        return 0
        #return self.sampler.sample(E)
    
    def __repr__(self):
        return "<Elastic >"
    
    def EXTRACT_ALL(self):
        
        dic = {
            "sampler": self.sampler,
               "imfp": self.imfp,
               "imfp0": self.imfp0,
               "imfp1": self.imfp1
               }
        return dic
            
            
        
        
        
        
        
        
        
        
        
# this thing really shouldn't be here, leaving it for now because I need
# to inhearit from DIST
 

def rebuildsFastCubicSpline(this):
    cdef sFastCubicSpline self = <sFastCubicSpline> sFastCubicSpline.__new__(sFastCubicSpline)
    self.N  = this.N  
    self.rc = this.rc
    from numpy import array
    self.LIMS = array(this.LIMS, dtype = int)

    self.mu_c = this.mu_c
    
    self.x = this.x 
    self.c = this.c     
    return self


cimport cython 
from libc.math cimport floor
cdef class sFastCubicSpline(DIST):
    cdef double dx
    cdef int N
    cdef double y
    cdef Py_ssize_t fr
    cdef double[:] x
    cdef double[:, :] c
    cdef int[:, :] LIMS
    cdef double r


    def __reduce__(self):
        this = MAP()
        from numpy import array
 
        this.N  = self.N  
        this.rc = self.rc
        this.LIMS = array(self.LIMS, dtype = int)

        this.mu_c = self.mu_c
        
        this.x = array(self.x )
        this.c =array( self.c)

        return rebuildsFastCubicSpline, (this, )


    def __init__(self, ndarray invcum, ndarray y, double rc):
        self.rc = rc
        
        
        
        
        
        lims = []
        cdef int i
        
        self.N = len(invcum)
        hashed = np.floor(invcum*self.N)
        
        
        
        
        
        
        
        #self.N -= 1
        # NOW DO THE SAME FOR MAIN
        # ALSO, TAKE PRINT OUT OF ELECTRON PARTICLE  
        indexes = np.arange(0, len(hashed))
        Imax = int(max(hashed))
        
        
        lims = [np.array([0, 0, 0], dtype = int)]
        
        
        
        for i in range(Imax + 1): #every possible value of the hash, index = hash
            selected = indexes[hashed == i]
            n = len(selected)
            if n == 0: # either out of bounds or no values in this range
                #if no values in this range, interpolate using last interval
                n_last = lims[-1][2]
                if n_last == 0: #out of bounds
                    lims.append(np.array([0, 0, 0], dtype = int))
                    continue
                i_last = lims[-1][1]
                lims.append(np.array([i_last, i_last, 1], dtype = int))
                continue

            lims.append(np.array([selected[0], selected[-1], n], dtype = int))
        
        self.LIMS = np.array(lims[1:])
        
        
        
        
        
        
        
        
        
        # hashed = np.floor(invcum*len(invcum))
        # indexes = np.arange(0, len(hashed))
         
        # Imax = int(max(indexes))
         
        # lims = []
        # for i in range(Imax + 2):
        #     _indexes = indexes[hashed == i]
             
        #     n = len(_indexes)
             
        #     if n < 2:
        #         lims.append(np.array([0, 0, n], dtype = int))
        #         continue
             
        #     _lims = [_indexes[0], _indexes[n-1], n]
        #     lims.append(np.array(_lims, dtype = int))
            
        # self.lims = np.array(lims)
        #print(np.array(self.lims))
        
        
        cdef object spline = CubicSpline(invcum, y)
        self.mu_c = spline(self.rc)
        
        self.x = spline.x
        assert len(self.x) == len(invcum)
        self.c = spline.c
        #self.N = len(spline.c[0]) -1    
        
        
        
    @cython.boundscheck(True)
    @cython.wraparound(False) 
    @cython.initializedcheck(True)
    @cython.cdivision(True)
    cdef double sample(self, mixmax_engine *genPTR):
        #LOG(f"sampler_rc(rc = {rc})")
        #cdef double r = urand()
        self.r = self.rc + (1 - self.rc)*genPTR.get_next_float()
        cdef int i = <int> floor(self.r*self.N)
        #self.fr = <Py_ssize_t> self.r*self.N
        
        

        
        #cdef int i = get_exp(E)
        
        if self.LIMS[i, 2] == 0:

            

            raise RuntimeError(f">>>>> {self.r < self.rc} OUT OF BOUNDS mu = {self.r}  rc = {self.rc} i = {i}")
        
        
        
        
        if self.LIMS[i, 2] == 1:
            return self.eval(self.LIMS[i, 0])
        
        
        
        
        if self.LIMS[i, 2] == 2:
            i = self.LIMS[i, 0]
            if self.r <= self.x[i+1]: return self.eval(i)
            return self.eval(i + 1)
        
        if self.LIMS[i, 2] == 3:
            i = self.LIMS[i, 0]
            if self.r <= self.x[i+1]: return self.eval(i)
            if self.r <= self.x[i+2]: return self.eval(i + 1)
            return self.eval(i + 2)
        
        if self.LIMS[i, 2] == 4:
            i = self.LIMS[i, 0]
            if self.r <= self.x[i + 1]: return self.eval(i)
            if self.r <= self.x[i + 2]: return self.eval(i + 1)
            if self.r <= self.x[i + 3]: return self.eval(i + 2)
            return self.eval(i + 3)
        
        cdef int START, END, MID
        START = self.LIMS[i, 0]
        END   = self.LIMS[i, 1]
        
        cdef double xMID
        
        #do binary search 
        while START <= END:
            #find middle
            MID = START + (END - START)//2 #prevents overflow somehow 
            
            xMID = self.x[MID]
            
            if self.r == xMID: #found the value
                return self.eval(MID)
            
            if self.r < xMID: # discard right side
                END = MID - 1 #do not include mid
                continue
            
            START = MID + 1
        return self.eval(END)
        
    @cython.boundscheck(True)
    @cython.wraparound(False) 
    @cython.initializedcheck(True)
    @cython.cdivision(True)
    cdef double eval(self, int i ):
       
        self.y = self.c[3, i]
        
        self.dx = self.r - self.x[i]
        self.y += self.c[2, i]*self.dx
        
        self.dx *= self.dx
        self.y += self.c[1, i]*self.dx
        
        self.dx *= self.dx
        self.y += self.c[0, i]*self.dx
        
        
        return self.y 
        
        
        
        
        
        
        
        
        
        
        
        

        
        # if self.lims[self.fr, 2] < 3:
        #     self.fr = self.lims[self.fr, 0]


        #     self.y = self.c[3, self.fr]
            
        #     self.dx = self.r - self.x[self.fr]
        #     self.y += self.c[2, self.fr]*self.dx
            
        #     self.dx *= self.dx
        #     self.y += self.c[1, self.fr]*self.dx
            
        #     self.dx *= self.dx
        #     self.y += self.c[0, self.fr]*self.dx
            
            
        #     return self.y            
            
            
            
        

        # self.fr = <Py_ssize_t> search._sortedArrayDOUBLE(self.x, 
        #                                     self.r, 
        #                                     self.lims[self.fr, 0], 
        #                                     self.lims[self.fr, 1] )
        

            
        # self.y = self.c[3, self.fr]

        # self.dx = self.r - self.x[self.fr]
        # self.y += self.c[2, self.fr]*self.dx
        
        # self.dx *= self.dx
        # self.y += self.c[1, self.fr]*self.dx
        
        # self.dx *= self.dx
        # self.y += self.c[0, self.fr]*self.dx
        
        
        # return self.y 
        
        
        
        



# cdef class DCS:
#     cdef list __container__
#     cdef vector[double] E
#     cdef LinLinInterpolation T
#     cdef int N
    
#     def __init__(self, Z):
#         self.__container__ = []
        
#         cdef dict xE, pE
#         el = self.getDCS(Z)
#         xE, pE = el.Y1, el.Y2
        
#         cdef list E = list(xE.keys())
#         self.N = len(E)
#         cdef double e

            
        

#         self.__container__ = []
#         cdef double[:] x
#         for e in E:
#             self.E.push_back(e)
#             x = xE[e]
            
            
#             p = CubicSpline(x, pE[e])
#             sampler = InvRationalInterpolation(p, min(x), max(x))
#             self.__container__.append(sampler)
        
        
#         TCS = self.getTCS(Z)
#         self.T = LinLinInterpolation(TCS.X, TCS.Y)
        
#     cdef InvRationalInterpolation _eval(self, double E):
#         i = search._sortedListDOUBLE(self.E, E, 0, self.N)
#         return self.__container__[i]
    
#     def eval(self, double E):
#         return self._eval(E)
    

        
    
#     @classmethod
#     def getDCS(cls, Z):
#         return db.EEDL(Z)[(9, 8, 0, 0, 9, 22)]
    
#     @classmethod
#     def getTCS(cls, Z):
#         return db.EEDL(Z)[(9, 8, 0, 0, 0, 0)]
    





# import numpy as np
# from numpy.random import rand



# cdef class Elastic:
#     cdef vector[double] coefs
#     cdef DCS[:] elements
#     cdef int N
#     cdef vector[double] cache
#     cdef double number_dens
#     cdef LinLinInterpolation[:] transport
    
#     def __init__(self, formula):
        
#         a = list(formula.values())
        
#         for x in a:
#             self.coefs.push_back(x)
#             self.cache.push_back(0.)
            

        
        
#         self.elements = array([DCS(Z) for Z in formula])
        
#         transport = []
        
#         for Z in formula:
#             data = db.EEDL(Z)[()]
#             transport.append(LinLinInterpolation(data.X, data.Y))
            
#         self.transport = array(transport)
        
#         self.N = len(self.elements)
#         self.number_dens = 1
        
        
#     cdef double _imfp(self, double E, double S):
#         cdef double result = 0
#         cdef int i
#         cdef double sigma_el, sigma_el1, chosen
        
        
#         for i in range(self.N):
            
#             sigma_el  = self.elements[i].T._eval(E) 
#             sigma_el1 = self.transport[i]._eval(E)
            
#             chosen = min(sigma_el, max(sigma_el1/C1, S/E / C2)) * self.coefs[i]
            
#             self.cache[i] = chosen
#             result += chosen
            
            
#         return result*self.number_density

#     def choose_element(self):
#         return self._choose_element()

    
#     def imfp(self, E):
#         return self._imfp(E)
        
    
#     cpdef _choose_element(self):
#         cdef double r = rand()*self.cache[self.N-1]
#         print(r)
#         cdef int i
#         for i in range(self.N):
#             if self.cache[i] < r:
#                 return self.elements[i]





# water = Elastic({1:2, 8:1})
# print(water.imfp(1))


        #self.sMFP1 = hLinLinInterpolation(eax, sMFP1)
        #self.sMFP2 = hLinLinInterpolation(eax, sMFP2)
        
        
        
        
        
        
                # imfp = []
        # cdef double e, SPtot, val0, val1, mfp0, mfp1
        # for e in eax:
        #     SPtot = SPcol(e) + SPrad(e)
            
        #     mfp0 = 1/self.imfp0._eval(e)
        #     mfp1 = 1/self.imfp1._eval(e)
            
            
        #     val0 = fmin(mfp1*C1, C2*e/SPtot)
        #     val = fmax(mfp0, val0)
        #     imfp.append(1/val)
            
        # imfp = np.array(imfp)




        
        #self.sampler = makeElastic.sampler(mu, eax, allDCS, SIGMA0)
        
        
        
        

        # print("        INTERPOLATING TRANSPORT CS as IMFP's")
        # print(f"        > E_min = {eax[0]}eV  |  E_max{eax[-1]}eV")
        # self.imfp0 = hLinLinInterpolation(eax, SIGMA0*formula.N)
        # self.imfp1 = hLinLinInterpolation(eax, SIGMA1*formula.N)
        
        
        # print("        MAKING HARD IMFP USING (SPrad + SPcol)soft")
        # print("        > eax = logspace(-3, 2, 2500)*1e6")
        # #
        # eax = logspace(-3, 2, 5500)*1e6
        # imfp = []
        
        
        
        
        
        
        # cdef double e, SPtot, val0, val1, mfp0, mfp1
        # for e in eax:
        #     SPtot = SPcol(e) + SPrad(e)
            
        #     mfp0 = 1/self.imfp0._eval(e)
        #     mfp1 = 1/self.imfp1._eval(e)
            
            
        #     val0 = fmin(mfp1*C1, C2*e/SPtot)
        #     val = fmax(mfp0, val0)
        #     imfp.append(1/val)
        
        
        # imfp = np.array(imfp)
        
        # print("        INTERPOLATING HARD IMFP")
        # self.imfp = hLinLinInterpolation(eax, imfp)