from numpy import * #array, geomspace, flip, load, searchsorted


#to be improved - this catches all exceptions
from . import electron as el
from . import photon   as ph
#import montecarlo.materials.electron as el
#import montecarlo.materials.photon   as ph


import os
directory = os.path.dirname(__file__)

def get_array(dic, key):
    """
    find closest key and interpolates corresponding data
    DEPRECATED: probably not going to use it now since I have the Distribution class. 
    """
    try:
        return array(dic[key])
    except:
        keys = list(dic.keys())
        i = searchsorted(keys, key, side = "left")
        key0, keyf = keys[i-1], keys[i]
        return .5*(array(dic[key0]) + array(dic[keyf]))




def count_spaces(string):
    """Count the number of spaces at the beggining of a string."""
    return len(string)-len(string.strip(' '))

class Material:
    """Container class for all materials."""
    def __init__(self, mat):
        """
        Initializes material. Reads and funnels data down to Electron and Photon.
        """
        self.name = mat
        self.path = directory + r"/" + mat + r"/" #to open material directory
        filename = mat + ".mat"                   #to open penelope input file


        with open(self.path + filename, "r") as file:
            text = file.readlines()

        text = [line.strip('\n') for line in text]



        ##################################################
        #GENERAL INFORMATION - all lines before first star
        pointer  = 1
        self.title = text[pointer][11:]

        pointer += 1
        self.density = float(text[pointer][16:-7])

        pointer += 1
        self.number_of_elements = int(text[pointer][37:])

        pointer += 1
        molecule = {} #'atomic_number':'number of atoms in molecule'
        for line in text[pointer:]:
            if count_spaces(line) != 3: break
            numbers_in_line = []
            for word in line.split():
                try:
                    numbers_in_line += [float(word.strip(","))]
                except ValueError:
                    pass
            element = numbers_in_line[0]
            number_of_elements = numbers_in_line[1]
            molecule[int(element)] = number_of_elements
            pointer += 1
        self.molecule = molecule

        line = text[pointer]
        self.mean_exitation_energy = float(line[25:-2])

        pointer += 1
        self.number_of_oscilators = [int(x) for x in text[pointer].split() if x.isdigit()][0]

        pointer += 1
        self.number_of_shells = [int(x) for x in text[pointer].split() if x.isdigit()][0]

        self.A = sum(A*freq for A, freq in self.molecule.items())
        self.N = 6.0221409e23*self.density/self.A
        ##################################################








        ##################################################
        bookmarks = []
        for n, line in enumerate(text):
            if line[0:4] == ' ***':
                bookmarks += [n]

        #gather all bookmarked text into a dict
        bookmark_ = bookmarks[0:-1]
        _ookmarks = bookmarks[1:]
        
        self.bookmarked_text = {}
        for i, j in zip(bookmark_, _ookmarks):
            self.bookmarked_text[text[i]] = text[i+1:j]



        #### ----> ELECTRONS <---- ####
        #Electron scaled bremss x-section
        text = self.getData('Electron scaled bremss x-section')

        bookmarks = []
        for n, line in enumerate(text):
            if count_spaces(line) == 1:
                bookmarks += [n]
        
        bookmarks_ = bookmarks[:-1]
        _bookmarks = bookmarks[1:]

        X = {}
        for start, end in zip(bookmarks_, _bookmarks):
            block = text[start:end]
            block = [float(value) for line in block for value in line.split()]
            E   = block[0]
            dcs = block[1:33]
            cs  = block[-1]
            X[E] = (cs, dcs)

        #Bremss angular distribution
        text = self.getData('Bremss angular distribution')
        ang = text
        
        #Stopping powers for electrons and positrons
        text = self.getData('Stopping powers for electrons and positrons')

        E, spC, spR = [], [], []
        for line in text:
            values = [float(value) for value in line.split()]
            E   += [values[0]]
            spC += [values[1]]
            spR += [values[2]]
        E, spC, spR = map(array, [E, spC, spR])
        sp_args = (E, spC, spR)

        #Electron and positron elastic cross sections
        text = self.getData('Electron and positron elastic cross sections')

        E, sigma, sigma1, sigma2 = [], [], [], []
        for line in text:
            values = [float(value) for value in line.split()]
            E       += [values[0]]
            sigma   += [values[1]]
            sigma1  += [values[2]]
            sigma2  += [values[3]]
        E, sigma, sigma1, sigma2 = map(array, [E, sigma, sigma1, sigma2])
        elastic_args = (E, sigma, sigma1, sigma2)

        #Electron elastic differential cross sections - contains integrated CS
        text = self.getData('Electron elastic differential cross sections')
        bookmarks = []
        for n, line in enumerate(text):
            if line[0:3]== ' -1':
                bookmarks += [n]
        bookmarks_ = bookmarks[:-1]
        _bookmarks = bookmarks[1:]

        
        DCS = {}
        for start, end in zip(bookmarks_, _bookmarks):
            block = text[start:end]
            block = [float(value) for line in block for value in line.split()]
            E   = block[1]
            sigma  = block[2]
            sigma1 = block[3]
            sigma2 = block[4]
            dcs = block[5:]
            DCS[E] = (sigma, sigma1, sigma2, dcs)
        #Electron ionisation cross sections
        text = self.getData('Electron ionisation cross sections')
        
        self.electron = el.Electron(X, ang,
                                    sp_args,
                                    elastic_args, DCS,
                                    self.N, self.density)







        


        #### ----> PHOTONS <---- ####
        text = self.getData('Compton and pair-production cross sections')
        text = self.getData('Compton and pair-production cross sections')
        
        
        text = self.getData('Rayleigh scattering')

        Y, X = [], []
        for line in text:
            numbers = line.split()
            X += [float(numbers[0])]
            Y += [float(numbers[1])]
        rayleigh = (X, Y)
        
        self.photon   = ph.Photon(self.path, self.density, rayleigh)
        ##################################################


    
    def getData(self, search):
        for key in self.bookmarked_text:
            if search in key:
                return self.bookmarked_text[key]

    def printKeys(self):
        for key in self.bookmarked_text:
            print(key)

















##class Photon(Material):
##    def __init__(self, path, density):
##        df = read_csv(path + "photon.txt", sep=" ")
##        data = array(df)[:,0:4]*density
##        E = data[:, 0]
##
##        coefs = data[:, 1]
##        self.rayleigh      = Rayleigh(E, coefs)
##
##        coefs = data[:, 2]
##        self.compton       = Compton(E, coefs)
##
##        coefs = data[:, 3]
##        self.photoelectric = Photoelectric(E, coefs)
##
##
##    def imfp(self, E):
##        index  = searchsorted(self.E, E, side="left")
##        result = .5*(self.coefs[index-1] + self.coefs[index])
##        return result
##        
##
##
##class Rayleigh(Photon):
##    def __init__(self, E, coefs):
##        self.coefs = coefs
##        self.E = E
##        pass
##    
##    def DCS(self, theta):
##        return (3/4)*(1 + cos(2*theta))
##
##class Compton(Photon):
##    def __init__(self, E, coefs):
##        self.coefs = coefs
##        self.E = E
##        pass
##    
##    def DCS(self, E, theta):
##        return 1/(1+E/0.511 * (1-cos(theta)))
##
##class Photoelectric(Photon):
##    def __init__(self, E, coefs):
##        self.coefs = coefs
##        self.E = E
##        pass


print("> Imported materials!")














if __name__ == '__main__':
    from matplotlib.pyplot import * 





if __name__=='__main__':
    water = Material("Water")
    ee = water.electron.elastic

    pr = water.photon.rayleigh
    pass
##    water = Material("Water")
##    #brain = Material("Brain")
##
##    def imfp(E):
##        return water.electron.brem.imfp(E)
##    
##    from numpy import *
##    x = 0.1
##    k = geomspace(x, 1, 32)
##    k = 1 + x- k
##    k = flip(k)
##
##    #Y = [imfp(E)/water.density for E in arange(1e3, 100e3, 1000)]
##    #DCS = water.electron.brem.DCS(1e6)
##    #plot(Y)
##    #yscale('log')
##    #show()
##    
##    img = imread("graph.png")
##    figure()
##    xlim(0,1), ylim(0,13)
##    imshow(img, extent=[0,1,0,13], aspect='auto')
##    
##    
##    #k = load("K.npy")
##    
##    
##    
##    #k = arange(32)/32
##    X1 = water.electron.brem.scaledDCS[1e3]
##    X2 = water.electron.brem.scaledDCS[10e3]
##    X3 = water.electron.brem.scaledDCS[100e3]
##
##    X4 = water.electron.brem.scaledDCS[1e6]
##    X5 = water.electron.brem.scaledDCS[10e6]
##    X6 = water.electron.brem.scaledDCS[100e6]
##
##
##    scatter(k, X1, s=1)
##    scatter(k, X2, s=1)
##    scatter(k, X3, s=1)
##    scatter(k, X4, s=1)
##    scatter(k, X5, s=1)
##    scatter(k, X6, s=1)
##
##    
##
##    show()

























##
##class Electron(Material):
##    def __init__(self, sp, brem, elastic):
##        
##        self.brem = Brem(brem)
##        self.elastic = Elastic(elastic)
##
##
##
##
##        self.density = sp[1]
##        df = read_csv(sp[0], sep=' ')
##        nistTable = array(df)[:,0:3]
##
##        sp = nistTable[:, 1] + nistTable[:, 2]
##        self.sp = CubicSpline(nistTable[:,0], sp*self.density)
##        
##        
##        
##
##class Brem(Electron):
##    k = load(directory + "/K.npy")
##
##    def __init__(self, brem):
##        #SCALED BREMSTRAUGHT SECTIONS X(k)
##        self.N = brem[1]
##        text   = brem[0]
##
##        CS = {}
##        scaledDCS = {}
##        for i, line in enumerate(text):
##            if count_spaces(line) == 1: #new energy entry
##                numbers_in_line = [float(number) for number in line.split()]
##                energy = numbers_in_line[0]
##                data = numbers_in_line[1:]
##                scaledDCS[energy] = data 
##            else:
##                numbers_in_line = [float(number) for number in line.split()]
##                scaledDCS[energy] += numbers_in_line
##
##
##        #the last number in each energy entry is the CS for that energy, it must be filtered out
##
##        CS = {energy:value[-1] for energy, value in scaledDCS.items()}
##          
##        self.CS = CS
##
##        energy = list(CS.keys())
##        CS     = array(list(CS.values()))
##
##        self.imfp = CubicSpline(energy, CS*1e-24*self.N)
##
##        
##        self.scaledDCS = {energy:dcs[:-1] for energy, dcs in scaledDCS.items()}
##
####    def imfp(self, E):
####        """Returns inverse of mean free path."""
####        
####        if E in self.CS: return self.CS[E]*self.N
####
####        closest_match = min(self.CS.keys(), key = lambda k: abs(k-E))
####        return 1e-24*self.CS[closest_match]*self.N
##
##
##    def X(self, E):
##        """Return callable - X(k) where k = W/E. and maximum height of graph for sampling."""
##        if E in self.scaledDCS: DCS = self.scaledDCS[E]
##        else:
##            closest_match = min(self.scaledDCS.keys(), key = lambda k: abs(k-E))
##            DCS = self.scaledDCS[closest_match]
##            
##        return (max(DCS), CubicSpline(self.k, DCS))
##
##    
##    def p(self, E):
##        """Return callable. p(theta) where theta is angle of emited photon."""
##        pass
##
##
##
##
##
##
##
##   #     elastic = (self.getData('Electron and positron elastic cross sections'),
##  #                 self.getData('Electron elastic differential cross sections'),
##   #                self.N)
##
##
##
##
##
##
##class Elastic(Electron):
##    def __init__(self, elastic):
##        #Electron and positron elastic cross sections
##        #self.N = elastic[-1]
##
##        textCS     = elastic[0]
##        textDCS    = elastic[1]
##        self.N     = elastic[2]
##        
##        self.E     = []
##        self.sigma = []
##
##        for line in textCS[1:]:
##            list_it   = line.split()
##            to_number = list(map(float, list_it))
##            self.E += [to_number[0]]
##            self.sigma += [to_number[1:4]]
##
##        self.E      = array(self.E)
##        self.sigma  = array(self.sigma)
##        self.L      = 1/(self.sigma*self.N)
##
##
##        C1 = 0.01
##        self.L_hard   = [max(L[0], C1*L[1]) for L in self.L]
##        self.mfp      = CubicSpline(self.E, self.L[:, 0])
##        self.mfp_hard = CubicSpline(self.E, self.L_hard)
##        #################################################
##
##
##        
##        
##
##        
##        #energy = array(list(self.CS.keys()))
##        #cross  = array(list(self.CS.values()))
##
##        #self.imfp = CubicSpline(energy, cross*self.N)
##
##        
##
##            
##
##        DCS    = elastic[1]
##        starts = []
##        for i, line in enumerate(DCS):
##            if line[0:3] == ' -1':
##                starts += [i]
##
##        starts1 = starts[0:-1]
##        starts2 = starts[1:]
##
##        
##        self._DCS = {}
##
##        for start,end in zip(starts1, starts2):
##            joining = ''
##            for i in range(start, end):
##                if i == start: joining += DCS[i][3:]
##                else:          joining += DCS[i]
##
##            joined  = joining.split()
##            numbers = list(map(float, joined))
##            self._DCS[numbers[0]] = numbers[1:]
##
##    def DCS(self, E):
##        """Give energy in eV."""
##        DCS = get_array(self._DCS, E)*self.N
##        N   = len(DCS)
##
##        thetas = arange(0, pi, pi/N)
##        func = CubicSpline(thetas, DCS)
##        return (max(DCS), func)
##
##    def theta_c(self, E):
##        from scipy.optimize import fsolve
##        from scipy.integrate import quad
##        
##        c         = 1/self.mfp_hard(E)/2/pi
##        integrand = self.DCS(E)[1]
##        F = lambda x: c - quad(integrand, x, pi)[0]
##        return fsolve(F,  0.02)
##        
##        
##
##
##    def P(l, x):
##        l += 1
##        coefs     = [0]*l
##        coefs[-1] = 1
##        return legval(x, coefs)
##
##
##
##    def plot(self, E):
##        """Plot DCS closest to E and prints CS closest to E."""
##        DCS = get_array(self._DCS, E)
##        #CS = get_array(self.CS, E)
##        N   = len(_DCS)
##        thetas = arange(0, pi, pi/N)
##
##        plot(thetas, DCS)
##
##
##
##    def show(self, E):
##        DCS = get_array(self.DCS, E)
##        CS  = get_array(self.CS, E)
##
##        print(DCS)
##        print(f"Total Cross Section = {CS}")
##
##
##
##        
##
##
##


##
##class StopingPower(Electron):
##    def __init__(self, sp):
##        self.density = sp[1]
##        df = read_csv(sp[0], sep=' ')
##        self.sp = array(df)[:,0:3]
##
##    def get(self, E):
##        """Give energy in MeV"""
##        
##        i = searchsorted(self.sp[:,0], E, side='left')
##        if self.sp[i, 0] == E: return self.sp[i,1]*self.density
##        return .5*(SP[i-1,1] + SP[i, 1])*self.density
##

















