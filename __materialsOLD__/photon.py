
#External Imports
from numpy import * 
from pandas import read_csv
from pynverse import inversefunc
from scipy.interpolate import CubicSpline
from scipy.integrate import *

#Internal Imports
from ..tools.distributions import *

class Photon:
    def __init__(self, path, density, F):
        df = read_csv(path + "photon.txt", sep=" ")
        data = array(df)[:,0:6]
        E = data[:, 0]

        coefs = data[:, 1]*density
        self.rayleigh = Rayleigh(E, coefs, F)
        
        coefs = data[:, 2]*density
        self.compton       = Compton(E, coefs)

        coefs = data[:, 3]*density
        self.photoelectric = Photoelectric(E, coefs)

        coefs = (data[:, 4] + data[:, 5])*density
        self.pairproduction = Pairproduction(E, coefs)
        


##
##    def imfp(self, E):
##        index  = searchsorted(self.E, E, side="left")
##        result = .5*(self.coefs[index-1] + self.coefs[index])
##        return result
    

class Interaction:
    def __init__(self, E, coefs):
        self.coefs = coefs
        self.E = E
        self.imfp = interp1d(E, coefs)
        
    def __call__(self, E):
        return self.imfp(E)


class Rayleigh(Interaction):
    """
    Calling this object with an energy value will return an inverse mean free path.
    Available tools:
    - self.pi -> form factor as a UnivariateDistribution (|F(x**2)|**2)
    - self.g  -> Unnormalized angular distribution of thomson scattering.
    """
    def __init__(self, E, coefs, F):
        super().__init__(E, coefs)

        xAxis = list(map(lambda x: x**2, F[0]))
        self.pi = UnivariateDistribution(xAxis = xAxis, yAxis = F[1])

    def g(self, c):
        return (1 + c**2)/2



class Compton(Interaction):
    """
    Compton interaction model.
    Note that if we ignore the init, there is not material dependence in the code.
    This is because the KN model assumes the electron is unbound.
    There is room for improvment here.
    """

    def __init__(self, E, coefs):
        super().__init__(E, coefs)



    def DCS(self, E, g):
        """
        Unnormalized differential cross section for a compton interaction.
        It is the klein nishina formula(without multiplicative constants).
        """
        P = self.P(E, g)
        return P**2 * (P + 1/P - 1+g**2) * sqrt(1-g**2)


    def P(self, E, g):
        """
        Returns fraction of the energy that the photon loses given that its
        incident energy is 'E' and its polar deflection is 'theta'.
        """

        return 1/(1+E/0.511 * (1 - g))
    



class Photoelectric(Interaction):
    def __init__(self, E, coefs):
        super().__init__(E, coefs)

        

class Pairproduction(Interaction):
    def __init__(self, E, coefs):
        super().__init__(E, coefs)
