from numpy import * #array, geomspace, flip, load, searchsorted

from . import electron as el
from .photon.photon import Photon


import os
directory = os.path.dirname(__file__)


def count_spaces(string):
    """Count the number of spaces at the beggining of a string."""
    return len(string) - len(string.strip(' '))


class Material:

    def __init__(self, Z, A, photon, electron):
        self.Z = Z
        self.A = A
        self.photon = photon
        self.electron = electron
        
    def create(self, density):
        """Create a material by providing this function with a density value."""
        self.photon.final_init(density)
        #self.electron.final_init(density)
        return self

    def __mul__(self, other):
        self.A = other*self.A

        
        photon = self.photon*other
        electron = self.electron#*other
        
        return Material(self.Z, self.A, photon, electron)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        A = self.A + other.A
        
        photon   = self.photon + other.photon
        electron = self.electron# + other.electron
        return Material(self.Z, self.A, photon, electron)
    
    def __radd__(self, other):
        return self.__add__(other)
        
        


class Element(Material):
    def __init__(self, Z, full_init = False):
        self.Z = Z
        
        #Evaluated Atomic Database Library
        EADL_path      = directory + r"/EADL/" + str(Z) + ".txt"
        with open(EADL_path, "r") as file:
            text = file.readlines()
            text = [line.strip('\n') for line in text]
            line = text[0]
            self.A = float(line[13:24])
            
        
        #Evaluated Photon Database Library
        EPDL_path      = directory + r"/EPDL/" + str(Z) + ".txt"
        self.EPDL_dict = self.get_bookmarked_text(EPDL_path)
        self.photon    = Photon(self, self.EPDL_dict)
        

        #Evaluated Electron Database Library
        EEDL_path      = directory + r"/EEDL/" + str(Z) + ".txt"
        self.EEDL_dict = self.get_bookmarked_text(EEDL_path)
        self.electron = None
        #self.electron  = ph.Electron(self.EEDL_dict)
        
        
        

    def get_bookmarked_text(self, path):
        """Reads EPICS file format and returns a dict with flags as keys
        and data as a list of strings."""
        with open(path, "r") as file:
            text = file.readlines()
            text = [line.strip('\n') for line in text]

            bookmarks = [0]
            
            for n, line in enumerate(text):
                if line == "                                                                       1":
                    bookmarks += [n + 1]
                    
            #gather all bookmarked text into a dict
            bookmark_ = bookmarks[0:-1]
            _ookmarks = bookmarks[1:]
        
            bookmarked_text = {}
            for i, j in zip(bookmark_, _ookmarks):
                line1, line2 = text[i], text[i+1]

                #on line 1
                Yi = float(line1[7:9])    #particle identifier
                Yo = float(line1[10:12])  #secondary particle designator
                Iflag = float(line1[31])  #interpolation flag
                
                #on line 2
                C  = float(line2[0:2])    #reaction descriptor
                I  = float(line2[2:5])    #reaction property
                S  = float(line2[5:8])    #reaction modifier
                X1 = float(line2[22:32])  #subshell designator
                
                flags = (Yi, C, S, X1, Yo, I)

                flags = tuple(map(int, flags))
                bookmarked_text[flags] = (Iflag, text[i+2:j-1])
                
        return bookmarked_text
    
    #to be removed I guess
    def getData(self, search):
        for key in self.bookmarked_text:
            if search in key:
                return self.bookmarked_text[key]

    def printKeys(self):
        for key in self.bookmarked_text:
            print(key)


print("> Imported materials!")
