# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 12:32:06 2020

@author: araujoj
"""

from pandas import read_csv #to be removed
from numpy import *#array,geomspace,flip,load,searchsorted
from numpy.polynomial.legendre import *
from scipy.interpolate import CubicSpline


    


#from scipy.integrate import quad

if __name__ == '__main__':
    from matplotlib.pyplot import * 



import os
directory = os.path.dirname(__file__)

def count_spaces(string):
    return len(string)-len(string.strip(' '))

def get_array(dic, key):
    """ find closest key and interpolates corresponding data"""
    try:
        return array(dic[key])
    except:
        keys = list(dic.keys())
        i = searchsorted(keys, key, side = "left")
        key0, keyf = keys[i-1], keys[i]
        return .5*(array(dic[key0]) + array(dic[keyf]))





class Material:
    def __init__(self, mat):
        self.name = mat
        self.path = directory + r"/" + mat + r"/" #to open material directory
        filename = mat + ".mat"                   #to open penelope input file


        #reading penelope input file @ material_direc/material_file.mat
        with open(self.path + filename, "r") as file:
            text = file.readlines() #outputs list of strings

        text = [line.strip('\n') for line in text]

        #locate information
        bookmarks = []
        for n,line in enumerate(text):
            if line[0:4] == ' ***':
                bookmarks += [n]

        bookmarks1 = bookmarks[0:-1]
        bookmarks2 = bookmarks[1:]

        #general information - all lines before first star
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

        self.data = {}
        for i, j in zip(bookmarks1, bookmarks2):
            self.data[text[i]] = text[i+1:j]

        #filter data, pass electron related data to electron instance
        electronData = {}
        for key in self.data:
            if 'electron' in key or 'Electron' in key:
                electronData[key] = self.data[key]

        #self.electron = Electron(self.path, electronData, self.density, self.N)
        self.photon   = Photon(self.path, self.density)


        sp   = (self.path + "electron.txt", self.density)
        brem = (self.getData('Electron scaled bremss x-section'), self.N)
        elastic = (self.getData('Electron and positron elastic cross sections'),
                   self.getData('Electron elastic differential cross sections'),
                   self.N)
        
        self.electron = Electron(sp, brem, elastic)


    
    def getData(self, search):
        for key in self.data:
            if search in key:
                return self.data[key]

    def printKeys(self):
        for key in self.data: print(key)













class Electron(Material):
    def __init__(self, sp, brem, elastic):
        
        self.brem = Brem(brem)
        self.elastic = Elastic(elastic)




        self.density = sp[1]
        df = read_csv(sp[0], sep=' ')
        nistTable = array(df)[:,0:3]

        sp = nistTable[:, 1] + nistTable[:, 2]
        self.sp = CubicSpline(nistTable[:,0], sp*self.density)
        
        
        

class Brem(Electron):
    k = load(directory + "/K.npy")

    def __init__(self, brem):
        #SCALED BREMSTRAUGHT SECTIONS X(k)
        self.N = brem[1]
        text   = brem[0]

        CS = {}
        scaledDCS = {}
        for i, line in enumerate(text):
            if count_spaces(line) == 1: #new energy entry
                numbers_in_line = [float(number) for number in line.split()]
                energy = numbers_in_line[0]
                data = numbers_in_line[1:]
                scaledDCS[energy] = data 
            else:
                numbers_in_line = [float(number) for number in line.split()]
                scaledDCS[energy] += numbers_in_line


        #the last number in each energy entry is the CS for that energy, it must be filtered out

        CS = {energy:value[-1] for energy, value in scaledDCS.items()}
          
        self.CS = CS

        energy = list(CS.keys())
        CS     = array(list(CS.values()))

        self.imfp = CubicSpline(energy, CS*1e-24*self.N)

        
        self.scaledDCS = {energy:dcs[:-1] for energy, dcs in scaledDCS.items()}

##    def imfp(self, E):
##        """Returns inverse of mean free path."""
##        
##        if E in self.CS: return self.CS[E]*self.N
##
##        closest_match = min(self.CS.keys(), key = lambda k: abs(k-E))
##        return 1e-24*self.CS[closest_match]*self.N


    def X(self, E):
        """Return callable - X(k) where k = W/E. and maximum height of graph for sampling."""
        if E in self.scaledDCS: DCS = self.scaledDCS[E]
        else:
            closest_match = min(self.scaledDCS.keys(), key = lambda k: abs(k-E))
            DCS = self.scaledDCS[closest_match]
            
        return (max(DCS), CubicSpline(self.k, DCS))

    
    def p(self, E):
        """Return callable. p(theta) where theta is angle of emited photon."""
        pass







   #     elastic = (self.getData('Electron and positron elastic cross sections'),
  #                 self.getData('Electron elastic differential cross sections'),
   #                self.N)






class Elastic(Electron):
    def __init__(self, elastic):
        #Electron and positron elastic cross sections
        #self.N = elastic[-1]

        textCS     = elastic[0]
        textDCS    = elastic[1]
        self.N     = elastic[2]
        
        self.E     = []
        self.sigma = []

        for line in textCS[1:]:
            list_it   = line.split()
            to_number = list(map(float, list_it))
            self.E += [to_number[0]]
            self.sigma += [to_number[1:4]]

        self.E      = array(self.E)
        self.sigma  = array(self.sigma)
        self.L      = 1/(self.sigma*self.N)


        C1 = 0.01
        self.L_hard   = [max(L[0], C1*L[1]) for L in self.L]
        self.mfp      = CubicSpline(self.E, self.L[:, 0])
        self.mfp_hard = CubicSpline(self.E, self.L_hard)
        #################################################


        
        

        
        #energy = array(list(self.CS.keys()))
        #cross  = array(list(self.CS.values()))

        #self.imfp = CubicSpline(energy, cross*self.N)

        

            

        DCS    = elastic[1]
        starts = []
        for i, line in enumerate(DCS):
            if line[0:3] == ' -1':
                starts += [i]

        starts1 = starts[0:-1]
        starts2 = starts[1:]

        
        self._DCS = {}

        for start,end in zip(starts1, starts2):
            joining = ''
            for i in range(start, end):
                if i == start: joining += DCS[i][3:]
                else:          joining += DCS[i]

            joined  = joining.split()
            numbers = list(map(float, joined))
            self._DCS[numbers[0]] = numbers[1:]

    def DCS(self, E):
        """Give energy in eV."""
        DCS = get_array(self._DCS, E)*self.N
        N   = len(DCS)

        thetas = arange(0, pi, pi/N)
        func = CubicSpline(thetas, DCS)
        return (max(DCS), func)

    def theta_c(self, E):
        from scipy.optimize import fsolve
        from scipy.integrate import quad
        
        c         = 1/self.mfp_hard(E)/2/pi
        integrand = self.DCS(E)[1]
        F = lambda x: c - quad(integrand, x, pi)[0]
        return fsolve(F,  0.02)
        
        


    def P(l, x):
        l += 1
        coefs     = [0]*l
        coefs[-1] = 1
        return legval(x, coefs)

    def F_soft(self, E, s,  theta):
        def summand(l, E):
            A = (2*l+1)/4/pi
            B = -P(1, cos(theta))*s/self.



    def plot(self, E):
        """Plot DCS closest to E and prints CS closest to E."""
        DCS = get_array(self._DCS, E)
        #CS = get_array(self.CS, E)
        N   = len(_DCS)
        thetas = arange(0, pi, pi/N)

        plot(thetas, DCS)



    def show(self, E):
        DCS = get_array(self.DCS, E)
        CS  = get_array(self.CS, E)

        print(DCS)
        print(f"Total Cross Section = {CS}")



        





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


















class Photon(Material):
    def __init__(self, path, density):
        df = read_csv(path + "photon.txt", sep=" ")
        data = array(df)[:,0:4]*density
        E = data[:, 0]

        coefs = data[:, 1]
        self.rayleigh      = Rayleigh(E, coefs)

        coefs = data[:, 2]
        self.compton       = Compton(E, coefs)

        coefs = data[:, 3]
        self.photoelectric = Photoelectric(E, coefs)


    def imfp(self, E):
        index  = searchsorted(self.E, E, side="left")
        result = .5*(self.coefs[index-1] + self.coefs[index])
        return result
        

print("imported materials")

class Rayleigh(Photon):
    def __init__(self, E, coefs):
        self.coefs = coefs
        self.E = E
        pass
    
    def DCS(self, theta):
        return (3/4)*(1 + cos(2*theta))

class Compton(Photon):
    def __init__(self, E, coefs):
        self.coefs = coefs
        self.E = E
        pass
    
    def DCS(self, E, theta):
        return 1/(1+E/0.511 * (1-cos(theta)))

class Photoelectric(Photon):
    def __init__(self, E, coefs):
        self.coefs = coefs
        self.E = E
        pass





















if __name__=='__main__':
    water = Material("Water")
    #brain = Material("Brain")

    def imfp(E):
        return water.electron.brem.imfp(E)
    
    from numpy import *
    x = 0.1
    k = geomspace(x, 1, 32)
    k = 1 + x- k
    k = flip(k)

    #Y = [imfp(E)/water.density for E in arange(1e3, 100e3, 1000)]
    #DCS = water.electron.brem.DCS(1e6)
    #plot(Y)
    #yscale('log')
    #show()
    
    img = imread("graph.png")
    figure()
    xlim(0,1), ylim(0,13)
    imshow(img, extent=[0,1,0,13], aspect='auto')
    
    
    #k = load("K.npy")
    
    
    
    #k = arange(32)/32
    X1 = water.electron.brem.scaledDCS[1e3]
    X2 = water.electron.brem.scaledDCS[10e3]
    X3 = water.electron.brem.scaledDCS[100e3]

    X4 = water.electron.brem.scaledDCS[1e6]
    X5 = water.electron.brem.scaledDCS[10e6]
    X6 = water.electron.brem.scaledDCS[100e6]


    scatter(k, X1, s=1)
    scatter(k, X2, s=1)
    scatter(k, X3, s=1)
    scatter(k, X4, s=1)
    scatter(k, X5, s=1)
    scatter(k, X6, s=1)

    

    show()
