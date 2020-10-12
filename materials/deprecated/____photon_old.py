
#External Imports
from numpy import *
from scipy.interpolate import *
from scipy.integrate import *

#Internal Imports
from ..tools.distributions import *

def getAxis(text):
    """Get axis from list of lines of text."""
    xAxis, yAxis = [], []
    for line in text:
        XY = [float(number) for number in line.split()]
        xAxis += [XY[0]]
        yAxis += [XY[1]]
    xAxis, yAxis = map(array, (xAxis, yAxis))
    return xAxis, yAxis


class Photon:
    def __init__(self, upper_dir, bookmarked_text, full_init = False):

        #COHERENT SCATTERING
        CS = bookmarked_text[(7, 71, 0, 0, 0, 0)]   #integrated cross section
        F  = bookmarked_text[(7, 93, 0, 0, 0, 941)] #form factor
        I  = bookmarked_text[(7, 93, 0, 0, 0, 943)] #imaginary anomalous scattering
        R  = bookmarked_text[(7, 93, 0, 0, 0, 944)] #real anomalous scattering
        E  = bookmarked_text[(7, 71, 0, 0, 7, 10)]  #average energy of the scattered photon
        
        self.coherent = Coherent(upper_dir, CS, F, R, I, E)

        #INCOHERENT SCATTERING
        CS = bookmarked_text[(7, 72, 0, 0, 0, 0)]    #integrated cross section
        S  = bookmarked_text[(7, 93, 0, 0, 0, 942)]
        
        self.incoherent = Incoherent(upper_dir, CS, S)



        #PHOTOELECTRIC
        CS = bookmarked_text[(7, 73, 0, 0, 0, 0)]

        self.photoelectric = Photoelectric(upper_dir, CS)



        #PAIR PRODUCTION
        CS = bookmarked_text[(7, 74, 0, 0, 0, 0)]
                             
        self.pairproduction = Pairproduction(upper_dir, CS)



        #TRIPLET PRODUCTION
        CS = bookmarked_text[(7, 75, 0, 0, 0, 0)]
                             
        self.tripletproduction = Tripletproduction(upper_dir, CS)
                             
        
    def __add__(self, other):
        self.coherent = self.coherent + other.coherent
        self.incoherent = self.incoherent + other.incoherent
        return self
        
    def __mul__(self, other):
        self.coherent = other*self.coherent
        self.incoherent = other*self.incoherent
        return self

    def __radd__(self, other):
        return self.__add__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def final_init(self, density):
        self.coherent.final_init(density)
        self.incoherent.final_init(density)







class TotalCrossSection:
    def __init__(self, upper_dir, CS):
        """
        Converts the text data to xAxis and yAxis arrays.
        Interpolates the CS data to create a standard xAxis grid.
        This standard grid is needed, otherwise we wouldn't be able to perform
        composition operations (2*H+O).
        """
        self.upper_dir = upper_dir
        self.Iflag  = CS[0]
        self.data   = CS[1]

        xAxis, yAxis = getAxis(self.data)

        assert self.Iflag == 2 or self.Iflag == 0
        
        #creating a standard grid (which is also denser)
        m, M = min(xAxis), max(xAxis)
        f = interp1d(xAxis, yAxis, assume_sorted = True)

        dx = (M-m)/600
        self.xAxis = arange(m, M+dx, dx)
        self.yAxis = f(self.xAxis)
        
        
            
    def final_init(self, density):
        Na = 6.0221409E23
        number_density = Na * density/self.upper_dir.A
        self.yAxis = number_density*self.yAxis
        self.imfp = interp1d(self.xAxis, self.yAxis, assume_sorted = True)
        
    def __call__(self, E):
        return self.imfp(E)

    def __add__(self, other):
        self.yAxis = self.yAxis + other.yAxis
        return self

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        self.yAxis = other*self.yAxis
        return self

    def __rmul__(self, other):
        return self.__mul__(other)





















class Coherent:
    """
    UPDATE DOCSTRING****
    Calling this object with an energy value will return an inverse mean free path.
    Available tools:
    - self.pi -> form factor as a UnivariateDistribution (|F(x**2)|**2)
    - self.g  -> Unnormalized angular distribution of thomson scattering.
    """
    
    def __init__(self, upper_dir, CS, F, R, I, E):
        
        self.CS = TotalCrossSection(upper_dir, CS)
        self.F  = FormFactor(F)
        


    def g(self, c):
        """Thomson DCS."""
        return (1 + c**2)/2

    

    def final_init(self, density):
        self.F.final_init()
        self.CS.final_init(density)

    def __add__(self, other):
        self.CS = self.CS + other.CS
        self.F  = self.F + other.F
        return self
    
    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        self.CS = self.CS*other
        self.F  = self.F*other
        return self
    
    def __rmul__(self, other):
        return self.__mul__(other)
        


class FormFactor(UnivariateDistribution):
    """Holds the form factor."""
    
    def __init__(self, F):
        self.Iflag = F[0]

        assert self.Iflag == 0 or self.Iflag == 2
        
        xAxis, yAxis = getAxis(F[1])
        
        f = interp1d(xAxis, yAxis, assume_sorted = True)

        m, M = min(xAxis), max(xAxis)
        dx = (M-m)/2000
        
        self.xAxis = arange(m, M + dx, dx)
        self.yAxis = f(self.xAxis)
        
        self.xAxis, self.yAxis = self.xAxis**2, self.yAxis**2


    def __add__(self, other):
        self.yAxis = self.yAxis + other.yAxis
        return self

    def __mul__(self, other):
        self.yAxis = other*self.yAxis
        return self

    def __rmul__(self, other):
        return self.__mul__(other)

    def __radd__(self, other):
        return self.__add__(other)

    def final_init(self, *args):
        super().__init__(xAxis = self.xAxis, yAxis = self.yAxis)
        #self.interpolation = interp1d(self.xAxis, self.yAxis)

    #def __call__(self, x2):
        #return self.interpolation(x2)
    
        
        
        
    











class Incoherent:
    
    def __init__(self, upper_dir, CS, S):
        self.CS = TotalCrossSection(upper_dir, CS)
        #self.S  = FormFactor(F)

    def __add__(self, other):
        self.CS = self.CS + other.CS
        #self.F  = self.F + other.F
        return self
    
    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        self.CS = self.CS*other
        #self.F  = self.F*other
        return self
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def final_init(self, density):
        self.CS.final_init(density)



class Photoelectric:
    
    def __init__(self, upper_dir, CS):
        self.CS = TotalCrossSection(upper_dir, CS)
        #self.S  = FormFactor(F)

    def __add__(self, other):
        self.CS = self.CS + other.CS
        #self.F  = self.F + other.F
        return self
    
    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        self.CS = self.CS*other
        #self.F  = self.F*other
        return self
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def final_init(self, density):
        self.CS.final_init(density)


class Pairproduction:
    
    def __init__(self, upper_dir, CS):
        self.CS = TotalCrossSection(upper_dir, CS)
        #self.S  = FormFactor(F)

    def __add__(self, other):
        self.CS = self.CS + other.CS
        #self.F  = self.F + other.F
        return self
    
    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        self.CS = self.CS*other
        #self.F  = self.F*other
        return self
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def final_init(self, density):
        self.CS.final_init(density)


class Tripletproduction:
    
    def __init__(self, upper_dir, CS):
        self.CS = TotalCrossSection(upper_dir, CS)
        #self.S  = FormFactor(F)

    def __add__(self, other):
        self.CS = self.CS + other.CS
        #self.F  = self.F + other.F
        return self
    
    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        self.CS = self.CS*other
        #self.F  = self.F*other
        return self
    
    def __rmul__(self, other):
        return self.__mul__(other)

    def final_init(self, density):
        self.CS.final_init(density)



##class Compton:
##    """
##    Compton interaction model.
##    Note that if we ignore the init, there is not material dependence in the code.
##    This is because the KN model assumes the electron is unbound.
##    There is room for improvment here.
##    """
##
##    def __init__(self, E, coefs):
##        
##
##
##
##    def DCS(self, E, g):
##        """
##        Unnormalized differential cross section for a compton interaction.
##        It is the klein nishina formula(without multiplicative constants).
##        """
##        P = self.P(E, g)
##        return P**2 * (P + 1/P - 1+g**2) * sqrt(1-g**2)
##
##
##    def P(self, E, g):
##        """
##        Returns fraction of the energy that the photon loses given that its
##        incident energy is 'E' and its polar deflection is 'theta'.
##        """
##
##        return 1/(1+E/0.511 * (1 - g))
    
##
##
##
##class Photoelectric(Interaction):
##    def __init__(self, E, coefs):
##        super().__init__(E, coefs)
##
##        
##
##class Pairproduction(Interaction):
##    def __init__(self, E, coefs):
##        super().__init__(E, coefs)
