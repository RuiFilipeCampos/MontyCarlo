# cython: profile=True

#External Imports
from numpy import arctan
import numpy as np

from libc.math cimport log
#Internal Imports
from ...tools.RITA import RationalInterpolation
from .CrossSection cimport CSLOGIC
from ...tools.interpol1 cimport LinLinInterpolation
from ...settings import __montecarlo__

path = __montecarlo__/'materials'/'photon'
path = str(path) + "/factor.txt"

# from .coherent.coherent                   import Coherent
# from .incoherent.incoherent               import Incoherent
# from .photoelectric.photoelectric         import Photoelectric
# from .pairproduction.pairproduction       import Pairproduction
# from .tripletproduction.tripletproduction import Tripletproduction



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
    





class Photon:
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

        self.coherent          = Coherent          (formula, density)
        self.incoherent        = Incoherent        (formula, density)
        #self.photoelectric     = Photoelectric     (formula, density)
        self.pairproduction    = Pairproduction    (formula, density)
        self.tripletproduction = Tripletproduction (formula, density)
        
    def __repr__(self):
        return "<Material.photon>"
    
    
    @property
    def ls(self):
        print(""".coherent
.incoherent    
.photoelectric  
.pairproduction 
.tripletproduction""")




cdef class Coherent(CSLOGIC):
    """
    PENELOPE SECTION: https://drive.google.com/file/d/1F-0JUO4Ucf_Z755IlMqJyz7T1pj7UbLz/
    """
    def __init__(self, formula, density):
        
        print("    INTERPOLATING COHERENT CS")
        
        super().__init__((7, 71, 0, 0, 0, 0), formula, density)
        
        
        
        print("    GETTING FF RITA INTERPOLATION <- to be improved")
        from .coherent import FormFactor
        self.FF = FormFactor.composeFF(formula)
        self.FF = RationalInterpolation(self.FF, 0, 300**2, True)
    
    def ls(self):
        print(".FF")

cdef class Incoherent(CSLOGIC):
    def __init__(self, formula, density):
        print("    INTERPOLATING Incoherent CS")
        super().__init__((7, 72, 0, 0, 0, 0), formula, density)
        
        print("    GETTING S fit")
        from .incoherent import IncoherentFormFactor
        self.S = IncoherentFormFactor.composeIFF(formula)
    
    def ls(self):
        print(".S")

class Photoelectric(CSLOGIC):
    """
    PENELOPE SECTION: https://drive.google.com/file/d/1GY5ZvvnZSyJDwedNNvB5qg52P8Jxh5Q1/
    """
    def __init__(self, formula, density):
        print("    INTERPOLATING Photoelectric CS")
        super().__init__((7, 73, 0, 0, 0, 0), formula, density)

cdef class Pairproduction(CSLOGIC):
    def __init__(self, formula, density):
        print("    INTERPOLATING Pairproduction CS")
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
        
        print(self.factor)
        self.CONST = 4*log(self.factor)


    def F0(self, double k):
        cdef double a = self.a
        cdef double a2 = self.a**2
        cdef double k2 = 2/k


        return    (-1.774 - 12.10*a  + 11.18*a2) * (k2)**.5  \
                + (8.523 - 73.26*a - 44.1*a2 )   * k2        \
                - (13.52 + 121.1*a - 96.41*a2)   * k2**(3/2) \
                + (8.946 + 62.05*a - 64.41 *a2)  * k2**2


  #      return  (-1.774 - 12.10*a  + 11.18*a2) * (k2)**.5    \
   #             + (8.523 +  73.26*a   − 44.41*a2) * (k2)        \
    #            - (13.52  + 121.1*a  − 96.41*a2) * (k2)**(3/2) \
     #           + (8.946  + 62.05*a  − 64.41*a2) * (k2)**2


    def g1(self, double b):
        return 7/3 - 2*log(1+b**2) - 6*b*arctan(b**-1) \
               -b**2 * (4-4*b*arctan(b**-1) - 3*log(1+b**-2))

    def g2(self, double b):
        return 11/6 - 2*log(1+b**2) - 3*b*arctan(b**-1) + .5*b**2*(4-4*b*arctan(b**-1) - 3*log(1+b**-2))

    def g0(self, double k):
        """
        FALTA O FATOR Rmec/h

        """
        return self.CONST-4*self.fC * self.F0(k)

    def b(self, double eps, double k):
        """
        FALTAM COISAS AQUI  <--- faltam? acho que ja as meti - CONFIRMAR (era o self.factor)
        """
        return self.factor*(2*k*eps*(1-eps))

    def getPhis(self, double eps, double k):
        g0 = self.g0(k)

        b = self.b(eps, k)
        g1 = self.g1(b)
        g2 = self.g2(b)

        return g1+g0, g2+g0




class Tripletproduction(CSLOGIC):
    def __init__(self, formula, density):
        print("    INTERPOLATING Tripletproduction CS")
        super().__init__((7, 75, 0, 0, 0, 0), formula, density)



















