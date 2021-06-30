

from ...settings import __montecarlo__


__path__ = __montecarlo__/'materials'/'electron'/'NRC_BREM'

__data__ = str(__path__) + '/nrc_brems_'

from numpy import array, zeros, searchsorted, log, logspace, insert, append
from matplotlib.pyplot import *
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d
from numpy import exp


def getX(Z):
    """
    Read X from NRC_BREM_Z file.
    Returns all X in units of *barn*.
    """
    with open(__data__ + str(Z), "r") as f:
        ds = f.read()

    dsigma = [float(x) for x in ds.split()]

    Xcoll = []
    X = []

    for i, x in enumerate(dsigma):
        if i % 54 == 0:
            Xcoll.append(array(X))
            X = [x]
        else: X.append(x)
            
    Xcoll.append(array(X))
    return array(Xcoll[1:])*1e-3

def makeX(formula):
    Zeq2, ntot = 0, 0
    X = zeros((49, 54))
    
    for Z, x in formula.items():
        ntot += x
        Zeq2 += x*Z**2
        
        X += getX(Z)*x*Z**2
    
    Zeq2 = Zeq2/ntot
    formula.Zeq2 = Zeq2
    formula.dcs = X/Zeq2
    return formula


class X:
    
    with open(__data__ + "grid", "r") as f:
        grid = f.read()
    
    grid = grid.split()[3:]
    E, k = [], []
    
    for i in range(49):
        val = float(grid[i])
        E.append(val)
        
    for i in range(49, 49 + 54):
        val = float(grid[i])
        k.append(val)
        
    k = array(k)
    E = array(E)
    
    
    def __init__(self, formula):
        
        self.Z = formula.Zeq2**.5
        self.ds = formula.dcs
        self.N = formula.N
        
        Z = [Z for Z in formula.keys()]
        n = [n for n in formula.values()]
        n = array(n)
        Z = array(Z)
        p = n/sum(n)
        self.Zeff = sum(p*Z*(Z+1))**.5
    

        

    def __iter__(self):
        yield from self.ds
        
    def __getitem__(self, i):
        return self.ds[i]
    
    
    def __call__(self, E):
        
        # a[i-1] <= v < a[i]
        i = searchsorted(self.E, E, side = "right")
        
        if i == 0:
            E2 = self.E[i]
            logE = log(E)
            log2 = log(E2)
            log1 = 0

            log_diff = log2 - log1

            #pi_1 = ( log2 - logE ) / log_diff
            pi_2 = (logE - log1)/log_diff            
            
                
            return self[0]*pi_2
        
        if E == self.E[i-1]:
            return self[i-1]
        
        dcs1, E1 = self[i-1], self.E[i-1]
        dcs2, E2 = self[i], self.E[i]
        
        assert E1 < E < E2
        
        logE = log(E)
        log2 = log(E2)
        log1 = log(E1)

        log_diff = log2 - log1

        pi_1 = ( log2 - logE ) / log_diff
        pi_2 = (logE - log1)/log_diff
        return pi_1*dcs1 + pi_2*dcs2
   
    def plot_dcs(self, i):
        figure()
        scatter(self.k, self[i], color = "b", s=1)
        plot(self.k, self[i], color = "b", lw = 1)
    
    def plotAll3d(self):
        fig = pyplot.figure()
        ax = Axes3D(fig)
        for i, ds in enumerate(self):
            n = len(ds)
            E = [self.E[i]]*n
            E = array(E)
            plot(log10(E), self.k, ds, color="black", lw=1)
        xlabel("log(E)");ylabel("k")        
   
    
   
    
   
    
   
    
   
def to_beta2(E):
    ELECTRON_REST_MASS = 0.51099895000
    X = (E/ELECTRON_REST_MASS + 1)**2
    return (X-1)/X  
   
    
from scipy.integrate import trapz
from numpy import exp
def F(Z, E):
    t = np.log(1 + 1e6*E/0.51099895000/Z**2)
    
    arg = -0.12359*t + 6.1274e-2*t**2 - 3.1516e-2 * t**3 + 7.7446e-3*t**4 - 1.0595e-3 * t**5 + 7.0568e-5*t**6 - 1.8080e-6*t**7
    return 1 - exp(arg)

def integrate1(molecule, E, k0, kf):

    if kf >= 1:
        kf = 1
    
    dcs = molecule(E)*F(molecule.Z, E)
    dcs = interp1d(molecule.k, dcs)
    
    k = molecule.k
    k = k[k0 < k]
    k = k[k < kf]
    k = insert(k, 0, k0)
    k = append(k, kf)
    
    dcs = dcs(k)
    
    
    return trapz(dcs, k)*E/to_beta2(E)
   

def getFULLSP(molecule, eax = logspace(-8, 3, 444)):
    sp = []
    
    for e in eax:
        spp = integrate1(molecule, e, 0, 1)

        sp.append(spp)
        
    return eax, array(sp)*molecule.N*molecule.Z**2 * 1e-24
    



def getSP(molecule, Wcr, eax = logspace(-8, 3, 444)):
    sp = []
    for e in eax:
        spp = integrate1(molecule, e, 0, Wcr/e)

        sp.append(spp)
        
    return eax, array(sp)*molecule.N*molecule.Z**2 * 1e-24
    
    
    
    
   
def integrate0(molecule, E, k0, kf):
    #print(k0, kf)
    if k0 >= 1:
        return 0
    
    

    #if k0 > 1:
    #    kf = 1
    
    dcs = molecule(E)*F(molecule.Z, E)
    dcs = interp1d(molecule.k, dcs)
    
    k = molecule.k
    if k0 < k[1]:
        k0 = k[1]
    
    
    k = k[k0 < k]
    k = k[k < kf]
    k = insert(k, 0, k0)
    k = append(k, kf)
    
    dcs = dcs(k)
    dcs = dcs/k
    
    return trapz(dcs, k)*to_beta2(E) 
   
    
def getIMFP(molecule, Wcr, eax = logspace(-3, 3, 699)):
    
    s0 = []
    for e in eax:
        
        s00 = integrate0(molecule, e, Wcr/e, 1)

        s0.append(s00)
        
    return eax, array(s0)*molecule.N*molecule.Z**2 * 1e-24




def integrate2(molecule, E, k0, kf):

    if kf >= 1:
        kf = 1
    
    dcs = molecule(E)*F(molecule.Z, E)
    dcs = interp1d(molecule.k, dcs)
    
    k = molecule.k
    k = k[k0 < k]
    k = k[k < kf]
    k = insert(k, 0, k0)
    k = append(k, kf)
    
    dcs = dcs(k)
    
    
    return trapz(dcs*k, k)*E**2/to_beta2(E)


def getSTRAGG(molecule, Wcr, eax = logspace(-3, 3, 699)):
    
    s0 = []
    for e in eax:
        
        s00 = integrate2(molecule, e, 0, Wcr/e)

        s0.append(s00)
        
    return eax, array(s0)*molecule.N*molecule.Z**2 * 1e-24