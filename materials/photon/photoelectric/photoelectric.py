#External Imports
from numpy import *
from scipy.interpolate import *
from scipy.integrate import *

#Internal Imports
from ....tools.CStools import TotalCrossSection
from ....tools.distributions import UnivariateDistribution
from ....tools.data import getAxis



shell_dict = {1:"K", 
              2:"L",  3:"L1",  4:"L23",   5:"L2", 6:"L3",
              7:"M",  8:"M1",  9:"M23",  10:"M2", 11:"M3", 12:"M45", 13:"M4", 14:"M5",
              15:"N", 16:"N1", 17:"N23", 18:"N2", 19:"N3", 20:"N45", 21:"N4", 22:"N5", 23:"N67", 24:"N6", 25:"N7",
              26:"O", 27:"O1", 28:"O23", 29:"O2", 30:"O3", 31:"O45", 32:"O4", 33:"O5", 34:"O67", 35:"O6", 36:"O7", 37:"O89", 38:"O8", 39:"O9"}



class Photoelectric:

    def __init__(self, upper_dir, CS, bookmarked_text):
        print("> *** photon/photoelectric: creating TotalCrossSection")
        self.CS = TotalCrossSection(upper_dir, CS)
        self.CS.final_init(upper_dir.density)























        

        print("> *** photon/photoelectric: TO DO -> fit subshell cross sect")

        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax  = plt.gca()

        ax.set_title(f"Integrated Cross Sections for subshells of Element {upper_dir.Z}.")
        ax.set_xlabel("Energy(MeV)")
        ax.set_ylabel("Cross Section (barn)")
        ax.set_yscale("log")
        ax.set_xscale("log")


        #getting CS data of shells
        xAxisCol, yAxisCol = [], []
        for key in bookmarked_text: #[(7, 73, 91, ?, 0, 0)]
            a, b, c, d, e, f = key
            if (a, b, c, e, f) == (7, 73, 91, 0, 0):
                xAxis, yAxis = getAxis(bookmarked_text[key][1])
                try:
                    label = shell_dict[d]
                except:
                    label = str(d)
                ax.plot(xAxis, yAxis, label=label)
                xAxisCol += [array(xAxis)]
                yAxisCol += [array(yAxis)]
                    # construct shells
        ax.legend()
        plt.show()
            #create new element
        self.elements = {Element(upper_dir.Z, xAxisCol, yAxisCol)}










    
    #COMPOSITION PHASE 2*H + O
    ######################################################################
    def __add__(self, other):
        new_self = Photoelectric.__new__(Photoelectric)


        #intersection = self.element & other.elements

        #for element in intersection:
        #    element.multiply(2)
        
        #new_self.elements = self.elements | other.elements
        new_self.CS       = self.CS + other.CS

        return new_self
    
    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        new_self = Photoelectric.__new__(Photoelectric)
        new_self.CS = self.CS*other

        #for element in self.elements:
            #element.multiply(other)

        #new_self.elements = self.elements

        return new_self
    
    def __rmul__(self, other):
        return self.__mul__(other)
    ######################################################################






    #FINALIZE
    #######################################################################
    def final_init(self, density):
        self.CS.final_init(density)

        #for element in self.elements:
            #element.final_init(density)
    ##########################################################################

    def __iter__(self):
        yield from self.elements




















class Element:
    #####################################################
    def __init__(self, Z, xAxisCol, yAxisCol):
        self.Z = Z
        self.yAxisCol = array(yAxisCol)
        self.xAxisCol = array(xAxisCol)

    ######################################################


    #######################################################
    #DEBUG
    def _scatterRAW(self):
        import matplotlib.pyplot as plt

        fig = plt.figure()

        for X, Y in zip(self.xAxisCol, self.yAxisCol):
            plt.plot(X, Y)

        plt.xlabel("Energy(MeV)")
        plt.ylabel("Cross Section(barns)")
        plt.title(f"CS of shells for element {self.Z}")
        plt.show()









    #enabling set operations
    #######################################################
    def __hash__(self):
        return self.Z

    def __eq__(self, other):
        return self.Z == other.Z
    ########################################################







    #construction phase
    ########################################################
    def multiply(self, other):
        self.yAxisCol = other*self.yAxisCol
    ########################################################


    ########################################################
    def final_init(self, density):

        #Na = 6.0221409 check this stuff - NIST might be using older definition of Na
        self.shells = []
        for xAxis, yAxis in zip(self.xAxisCol, self.yAxisCol):
            self.shells += [interp1d(xAxis, yAxis, bounds_error=False, fill_value = 0.)]

        #Na = 6.02214076E23
        #number_density = Na * density/self.upper_dir.A * 1E-24
        #self.imfp = lambda x: number_density*self.interpolation(x)
    ########################################################

    def __iter__(self):
        yield from self.shells