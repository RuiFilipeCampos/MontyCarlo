# cython: profile=False
print("Importing `.material.photon.pyx`")



cdef extern from "math.h":
    bint isnan(double x)

#External Imports
from numpy import arctan
import numpy as np

from libc.math cimport log
#Internal Imports
from ...tools.RITA import RationalInterpolation
from .CrossSection cimport CSLOGIC
from ...tools.interpol1 cimport LinLinInterpolation
from ...settings import __montecarlo__

from .. import database as db

path = __montecarlo__/'materials'/'photon'
path = str(path) + "/factor.txt"

# from .coherent.coherent                   import Coherent
# from .incoherent.incoherent               import Incoherent
# from .photoelectric.photoelectric         import Photoelectric
# from .pairproduction.pairproduction       import Pairproduction
# from .tripletproduction.tripletproduction import Tripletproduction




import numpy as np


h = 4.135667696e-15 #eV s
m_e = 9.10938e-31 #kg
c = 299792458 #m/s

C = h/(m_e*c)
CONST = C**2 * 2e-4
data = db.EPDL[0][(7, 93, 0, 0, 0, 941)]


XX = np.array(data.X, dtype = np.longdouble)
NN = len(data.X)
XX = np.array(CONST*XX**2, dtype =np.longdouble)




cdef long double[::1] XXX = XX

cdef extern from "<math.h>" nogil:
    double frexp(double x, int* exponent)


cdef int get_exp(double x):
    cdef int exp;
    frexp(x, &exp);
    return exp;


cimport cython
from scipy.interpolate import CubicSpline

hashed =  np.array([get_exp(x) for x in XX], dtype = np.int32)
indexes = np.arange(0, len(hashed))
Imax = int(max(hashed))



lims = [np.array([0, 0, 0], dtype = np.int32)]


cdef int i 
for i in range(Imax + 1): #every possible value of the hash, index = hash
    selected = indexes[hashed == i]
    n = len(selected)
    if n == 0: # either out of bounds or no values in this range
        #if no values in this range, interpolate using last interval
        n_last = lims[-1][2]
        if n_last == 0: #out of bounds
            lims.append(np.array([0, 0, 0], dtype = np.int32))
            continue
        i_last = lims[-1][1]
        lims.append(np.array([i_last, i_last, 1], dtype = np.int32))
        continue

    lims.append(np.array([selected[0], selected[-1] , n], dtype = np.int32))

cdef int[:, ::1] LIMS = np.array(lims[1:], dtype = np.int32) ### memory view defined in pxd, cdef double[::1] EAX











with open(path, "r") as f:
    text = f.readlines()
    
Z = []
factor = []
for line in text:
    numbers = line.split()
    Z.append(float(numbers[0]))
    factor.append(float(numbers[1]))


newZ = []
newfactor = []
for tup in sorted(zip(Z, factor)):
    newZ.append(tup[0])
    newfactor.append(tup[1])



cdef LinLinInterpolation factors = LinLinInterpolation(np.array(newZ), np.array(newfactor))
    




@cython.auto_pickle(True)
cdef class Photon:
    """
    DATA BASE DOC:
        https://www-nds.iaea.org/epics/DOCUMENTS/ENDL2002.pdf

    REFERENCES:
        PENELOPE
        https://drive.google.com/file/d/1rb_wKkICOyL5UMuG4chuRxBQHuqR_8q1/
        -----------------------------------------------------------------------------------------------
        ??? :: Atomic Form Factors, Incoherent Scattering Functions,and Photon Scattering Cross Sections
        https://drive.google.com/file/d/1hbhDTCn1NGYZIB31K-OW21RsbFEPO9ta/
        ->  by: Hubble et al.
        -----------------------------------------------------------------------------------------------
        photonCS1994 :: Analytical cross sections for Monte Carlo simulation of photon transport.
        https://drive.google.com/file/d/1Rt2DqkhwJINQC1S469adqn5Whz9110ep/view?usp=sharing
        ->  savat et. al
        -----------------------------------------------------------------------------------------------

    """
    def __init__(self, formula, density):
        self.density = density


        formula.log.add_paragraph("     > Note: Photoelectric has been fully compiled by 'Molecule'")

        formula.log.add_paragraph("     > Constructing Coherent")
        self.coherent = Coherent(formula, density)

        formula.log.add_paragraph("     > Constructing Incoherent")
        self.incoherent = Incoherent(formula, density)

        formula.log.add_paragraph("     > Constructing Pairproduction")
        #self.photoelectric     = Photoelectric     (formula, density)
        self.pairproduction = Pairproduction(formula, density)

        formula.log.add_paragraph("     > Constructing Tripletproduction")
        self.tripletproduction = Tripletproduction(formula, density)

    def __repr__(self):
        return "<Material.photon>"
    
    @property
    def COHERENT(self):
        return self.coherent
    
    @property
    def ls(self):
        print("""
.coherent
.incoherent    
.photoelectric  
.pairproduction 
.tripletproduction
""")



from ...tools.main cimport remove_duplicates



def reconstruct_Coherent(imfpA, imfpB, xSPLINE, ySPLINE, X, Y, xLIMS, yLIMS, xMAX, xMIN, xADDER, yADDER):
    cdef Coherent new = Coherent({}, 0, pickle = True)
    new.imfpA = imfpA
    new.imfpB = imfpB
    new.xSPLINE = xSPLINE
    new.ySPLINE =ySPLINE
    new.X = X
    new.Y = Y
    new.xLIMS = xLIMS
    new.yLIMS = yLIMS

    return new


@cython.auto_pickle(True)
cdef class Coherent(CSLOGIC):
    """
    PENELOPE SECTION: https://drive.google.com/file/d/1F-0JUO4Ucf_Z755IlMqJyz7T1pj7UbLz/
    """
    
    def __reduce__(self):
        imfpA = np.array(self.imfpA)
        imfpB = np.array(self.imfpB)
        xSPLINE = np.array(self.xSPLINE)
        ySPLINE = np.array(self.ySPLINE)
        X, Y = np.array(self.X), np.array(self.Y)
        xLIMS, yLIMS = np.array(self.xLIMS, dtype = np.int32), np.array(self.yLIMS, dtype = np.int32)
        tup = (imfpA, imfpB, xSPLINE, ySPLINE, X, Y, xLIMS, yLIMS, self.xMAX, self.xMIN, self.xADDER, self.yADDER)
        to_pickle = list(tup)
        return (reconstruct_Coherent,  tup)
    
    def __init__(self, formula, density, pickle = False):
        
        
        if pickle:
            return
        formula.log.add_paragraph("         INTERPOLATING COHERENT IMFP")
        
        super().__init__((7, 71, 0, 0, 0, 0), formula, density)
        
        h = 4.135667696e-15 #eV s
        m_e = 9.10938e-31 #kg
        c = 299792458 #m/s

        C = h/(m_e*c)
        CONST = C**2 * 2e-4
        
        SPLINES = []
        
        for Z, x in formula.items():
            data = db.EPDL[Z-1][(7, 93, 0, 0, 0, 941)]
            FF2 = x*data.Y**2
            
            X = CONST*data.X**2 #k*(1 - cos)
            SPLINES.append(CubicSpline(X, FF2))
            
        
        Y = np.zeros(NN, dtype = np.longdouble)
        
        for FF2 in SPLINES:
            Y += np.array(FF2(XX), dtype = np.longdouble)
            
        from scipy.integrate import trapz, cumtrapz
        
        

        
        norm = trapz(Y, XX)
        Y = Y/norm
        
        cdef ndarray[double] CUMUL
        
        CUMUL = cumtrapz(Y, XX, initial = 0)
        
        for exclude_from in range(len(CUMUL)):
            if CUMUL[exclude_from] == 1:
                break
            
        CUMUL = CUMUL[:exclude_from + 1]
        myX = XX[:exclude_from + 1]
        
        

        
        CUMUL, myX = remove_duplicates(CUMUL, myX)


        if len(CUMUL) != len(myX):
            formula.log.add_paragraph("         > Warning: line 305")

            print(">>>>>>>>>>>>> you'l need to patch this eventually, line 263 on `material.photon.coherent`")
            
            __XX, __YY = [], []
            for __xx, __yy in zip(myX, CUMUL):
                __XX.append(__xx)
                __YY.append(__yy)
            
            myX = np.array(__XX)
            CUMUL = np.array(__YY)
      #  for x in CUMUL: print(x)
        
        cumul = CubicSpline(myX, CUMUL)
        
        self.ySPLINE  = cumul.c
        xADDER, xLIMS = self.construct_LIMS(myX)
        self.xADDER   = xADDER
        self.xLIMS    = xLIMS
        self.X        = np.array(myX, dtype = float)
        self.xMAX     = max(myX)
        self.xMIN     = min(myX)
        

        dx = np.diff(CUMUL)
        if np.any(dx <= 0):
            print("\n\n > DEBUG INFO BEFORE RAISING `ValueError`")
            print("> The following array should be strictly increasing:")
            print(CUMUL)
            print("> If this has failed, add this case to the `/tests/test_tools.py` script.")
            raise ValueError("`x` must be strictly increasing sequence.")
        del dx

        invcumul = CubicSpline(CUMUL, myX)

        self.xSPLINE = invcumul.c
        yADDER, yLIMS = self.construct_LIMS(CUMUL)
        self.yADDER = yADDER
        self.yLIMS = yLIMS
        self.Y = CUMUL
        
        

        

    
        # hashed =  np.array([get_exp(x) for x in CUMUL], dtype = int)
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
        
        # self.LIMS = np.array(lims[1:], dtype = int) ### memory view defined in pxd, cdef double[::1] EAX
            
        
    def construct_LIMS(self, ARR):
        hashed =  np.array([get_exp(x) for x in ARR], dtype = np.int32)

        adder = min(hashed)
        if adder < 0: adder = abs(adder)
        else: adder = 0
        hashed += adder
        indexes = np.arange(0, len(hashed))
        Imax = int(max(hashed))
        
        
        
        lims = [np.array([0, 0, 0], dtype = np.int32)]
        
        
        cdef int i 
        for i in range(Imax + 1): #every possible value of the hash, index = hash
            selected = indexes[hashed == i]
            n = len(selected)
            if n == 0: # either out of bounds or no values in this range
                #if no values in this range, interpolate using last interval
                n_last = lims[-1][2]
                if n_last == 0: #out of bounds
                    lims.append(np.array([0, 0, 0], dtype = np.int32))
                    continue
                i_last = lims[-1][1]
                lims.append(np.array([i_last, i_last, 1], dtype = np.int32))
                continue
        
            lims.append(np.array([selected[0], selected[-1] , n], dtype = np.int32))

        return adder, np.array(lims[1:], dtype = np.int32) ### memory view defined in pxd, cdef double[::1] EAX

        
        
    cdef int find_index_Y(self, double x):
        cdef int i;
        frexp(x, &i);


        i += self.yADDER
        
        if i < 0: i = 0
        
        cdef int k = self.yLIMS[i, 2]
        if k is 1:
            return self.yLIMS[i, 0]
        
        if k == 2:
            i = self.yLIMS[i, 0]
            if x <= self.Y[i + 1]: return i
            return i + 1
        
        if k is 3:
            i = self.yLIMS[i, 0]
            if x <= self.Y[i + 1]: return i
            if x <= self.Y[i + 2]: return i + 1
            return i + 2
        
        if k is 4:
            i = self.yLIMS[i, 0]
            if x <= self.Y[i + 1]: return i
            if x <= self.Y[i + 2]: return i + 1
            if x <= self.Y[i + 3]: return i + 2
            return i + 3
        
        cdef int START, END, MID
        START = self.yLIMS[i, 0]
        END   = self.yLIMS[i, 1] 
        
        cdef double xMID
        
        #do binary search 
        while START <= END:
            #find middle
            MID = START + (END - START)//2 #prevents overflow somehow 
            
            xMID = self.Y[MID]
            
            if x is xMID: #found the value
                return MID
            
            if x < xMID: # discard right side
                END = MID - 1 # do not include mid
                continue
            
            START = MID + 1
        return END

    cdef int find_index_X(self, double x):

        cdef int i;
        frexp(x, &i);
        
        i += self.xADDER

        

        cdef int k = self.xLIMS[i, 2]
        if k is 1:
            return self.xLIMS[i, 0]
        
        if k == 2:
            i = self.xLIMS[i, 0]
            if x <= self.X[i + 1]: return i
            return i + 1
        
        if k is 3:
            i = self.xLIMS[i, 0]
            if x <= self.X[i + 1]: return i
            if x <= self.X[i + 2]: return i + 1
            return i + 2
        
        if k is 4:
            i = self.xLIMS[i, 0]
            if x <= self.X[i + 1]: return i
            if x <= self.X[i + 2]: return i + 1
            if x <= self.X[i + 3]: return i + 2
            return i + 3
        
        cdef int START, END, MID
        START = self.xLIMS[i, 0]
        END   = self.xLIMS[i, 1] 
        
        cdef double xMID
        
        #do binary search 
        while START <= END:
            #find middle
            MID = START + (END - START)//2 #prevents overflow somehow 
            
            xMID = self.X[MID]
            
            if x is xMID: #found the value
                return MID
            
            if x < xMID: # discard right side
                END = MID - 1 # do not include mid
                continue
            
            START = MID + 1
        return END
        

    def public_evalY(self, x):
        return self.evalY(x)
    
    def public_evalX(self, x):
        return self.evalX(x)
    
    
    @cython.initializedcheck(False)
    @cython.cdivision(True)
    @cython.boundscheck(False)  # Deactivate bounds checking
    @cython.wraparound(False)   # Deactivate negative indexing.
    cdef double evalY(self, double x):
        
        
        
        
        if x < self.xMIN:
            return 0
        
        if x > self.xMAX:
            return 1
        
        cdef int i = self.find_index_X(x)

        cdef double y = 0
        cdef double xi = self.X[i]
        cdef int k
        for k in range(0, 4):
            y += self.ySPLINE[k, i]*(x - xi)**(3-k)
        return y
    

    cdef double evalX(self, double r): 
        #cdef double rc = self.eval_CUMUL(qmax2)

            
 
   
        
        cdef int  i, k
        cdef double y, xi
        try:
            i = self.find_index_Y(r)
            
            if i < 0: print(i)
    
            y = 0
            xi = self.Y[i]
            k
            for k in range(0, 4):
                y += self.xSPLINE[k, i]*(r - xi)**(3-k)
            return y
        except IndexError:
            import time
            print(r, get_exp(r), self.yADDER)
            return 0

        
    
        
      #   print("    GETTING FF RITA INTERPOLATION <- to be improved")
      #   from .coherent import FormFactor
      #   FF = FormFactor.composeFF(formula)
        
      # #   from ...tools.CubicInverseTransform import fromCallable
        
      #   # self.FFnew = fromCallable(FF, 0, 10, num = 2500)
      #   import numpy as np 
      #   for x in np.arange(1. ,1000):
      #       if FF(x) < 1e-16: break
      #   self.FF = RationalInterpolation(FF, 0, x, True)
      #   print(self.FF.invCum(0.5))
    
    def ls(self):
        print(".FF")



def reconstruct_Incoherent(imfpA, imfpB):
    cdef Incoherent new = Incoherent({}, 0, pickle = True)
    new.imfpA = imfpA
    new.imfpB = imfpB
    return new
    

from ...tools.interpol1 cimport CSa
@cython.auto_pickle(True)
cdef class Incoherent(CSLOGIC):
    def __init__(self, formula, density, pickle = False):
        
                
        if pickle:
            return

        formula.log.add_paragraph("         INTERPOLATING Incoherent IMFP")

        super().__init__((7, 72, 0, 0, 0, 0), formula, density)
        
        
        
        
        
       # self.S = CSa(x, y)
        
        
       # print("    GETTING S fit")
       # from .incoherent import IncoherentFormFactor
      #  self.S = IncoherentFormFactor.composeIFF(formula)
    def __reduce__(self):
        imfpA = np.asarray(self.imfpA)
        imfpB = np.asarray(self.imfpB)
        return (reconstruct_Incoherent,( imfpA, imfpB) )
    
    
    def ls(self):
        print(".S")

class Photoelectric(CSLOGIC):
    """
    PENELOPE SECTION: https://drive.google.com/file/d/1GY5ZvvnZSyJDwedNNvB5qg52P8Jxh5Q1/
    """
    def __init__(self, formula, density):
        print("    INTERPOLATING Photoelectric CS")
        super().__init__((7, 73, 0, 0, 0, 0), formula, density)




def reconstruct_Pairproduction(imfpA, imfpB, factor, CONST, alpha, a, Zeq, fC):
    cdef Pairproduction new = Pairproduction({}, 0, pickle = True)
    new.imfpA = imfpA
    new.imfpB = imfpB
    
    new.factor = factor
    new.CONST = CONST
    new.alpha = alpha
    new.a = a
    new.Zeq = Zeq
    new.fC = fC
    return new

@cython.auto_pickle(True)
cdef class Pairproduction(CSLOGIC):
    def __init__(self, formula, density, pickle = False):  
        if pickle:
            return

        formula.log.add_paragraph("         INTERPOLATING Pairproduction IMFP")
        super().__init__((7, 74, 0, 0, 0, 0), formula, density)

        from ..database import EADL

        Zeq = 0
        Am  = 0

        for Z in formula:
            x = formula[Z]
            Aw = EADL[Z-1]['Aw']

            Am  += x*Aw
            Zeq += x*Z*Aw

        self.Zeq = Zeq/Am

        self.alpha = 1/137.035999074  #constante de estrutura fina?
        self.a     = self.alpha*self.Zeq

        self.fC = self.a**2 *( (1+self.a**2)**-1     \
                               + 0.202059            \
                               - 0.03693*self.a**2   \
                               + 0.00835*self.a**4 \
                               - 0.00201*self.a**6 \
                               + 0.00049*self.a**8   \
                               - 0.00012*self.a**10\
                               + 0.00003*self.a**12  )
            
        self.factor = factors._eval(self.Zeq)
        formula.log.add_paragraph(f"         >>>> Zeq = {self.Zeq}   | factor = {self.factor}")

        self.CONST = 4*log(self.factor)

    def __reduce__(self):
        imfpA = np.asarray(self.imfpA)
        imfpB = np.asarray(self.imfpB)
        tup = (imfpA, imfpB, self.factor, self.CONST, self.alpha, self.a, self.Zeq, self.fC)
        return (reconstruct_Pairproduction, tuple(tup))
        
        
        
        
        
    cpdef double F0(self, double k):
        cdef double a = self.a
        cdef double a2 = self.a**2
        cdef double k2 = 2/k


        return    (-1.774 - 12.10*a  + 11.18*a2) * (k2)**.5  \
                + (8.523 + 73.26*a - 44.1*a2 )   * k2        \
                - (13.52 + 121.1*a - 96.41*a2)   * k2**(3/2) \
                + (8.946 + 62.05*a - 63.41 *a2)  * k2**2


  #      return  (-1.774 - 12.10*a  + 11.18*a2) * (k2)**.5    \
   #             + (8.523 +  73.26*a   − 44.41*a2) * (k2)        \
    #            - (13.52  + 121.1*a  − 96.41*a2) * (k2)**(3/2) \
     #           + (8.946  + 62.05*a  − 64.41*a2) * (k2)**2


    cpdef double g1(self, double b):
        return 7/3 - 2*log(1+b**2) - 6*b*arctan(b**-1) \
               -b**2 * (4-4*b*arctan(b**-1) - 3*log(1+b**-2))

    cpdef double g2(self, double b):
        return 11/6 - 2*log(1+b**2) - 3*b*arctan(b**-1) + .5*b**2*(4-4*b*arctan(b**-1) - 3*log(1+b**-2))

    cpdef double g0(self, double k):
        """
        FALTA O FATOR Rmec/h

        """
        return self.CONST-4*self.fC + self.F0(k)

    cpdef double b(self, double eps, double k):
        """
        FALTAM COISAS AQUI  <--- faltam? acho que ja as meti - CONFIRMAR (era o self.factor) <--- esta correto
        """
        return self.factor/(2*k*eps*(1-eps))

    cpdef (double, double) getPhis(self, double eps, double k):
        g0 = self.g0(k)

        b = self.b(eps, k)
        g1 = self.g1(b)
        g2 = self.g2(b)
        
        res1 = g1 + g0
        res2 = g2+g0
        
        return res1 if res1 > 0 else 0, res2 if res2 > 0 else 0



def reconstruct_Tripletproduction(imfpA, imfpB):
    cdef Tripletproduction new = Tripletproduction({}, 0, pickle = True)
    new.imfpA = imfpA
    new.imfpB = imfpB
    return new

@cython.auto_pickle(True)
cdef class Tripletproduction(CSLOGIC):
    def __init__(self, formula, density, pickle = False):
        
        
        
                
        if pickle:
            return
            
        formula.log.add_paragraph("         INTERPOLATING Tripletproduction IMFP")

        super().__init__((7, 75, 0, 0, 0, 0), formula, density)

    def __reduce__(self):
        imfpA = np.asarray(self.imfpA)
        imfpB = np.asarray(self.imfpB)
        return (reconstruct_Tripletproduction, (imfpA, imfpB) )
        


















