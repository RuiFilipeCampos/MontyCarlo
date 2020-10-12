from pandas import read_csv
from numpy import *

class Materials:
    def __init__(self, path, density):
        self.path     = path
        self.rho      = density
        self.photon   = Photon(path, density)
        self.electron = Electron(path)

class Photon:
    def __init__(self, path, density):
        df = read_csv(path + "photon.txt", sep=" ")
        self.coef = array(df)[:,0:4]*density
    
    def get_coefs(self, E):
        index = searchsorted(self.coef[:,0], E, side="left") #finds index where element can be inserted without changing order
        result = (self.coef[index-1, :] +
                  self.coef[index  , :])/2 #average the entries
        return(result[1:])

class Electron:
    def __init__(self, path):
        self.path = path
        df = read_csv(path + "photon.txt", sep=" ") ## change this to electrons!!!!
        self.coef = array(df)[:,0:4]
        
    def get_coefs(self, E):
        index = searchsorted(self.coef[:,0], E, side="left") #finds index where element can be inserted without changing order
        result = (self.coef[index-1, :] +
                  self.coef[index  , :])/2 #average the entries
        return(result[1:])

    def get_phaseFunc(self):
        #return func that can be evaluated at theta and phi
        pass




