# distutils: language = c++

__doc__ = """

MAIN SOURCE OF INFO:
    PENELOPE 11  - https://drive.google.com/file/d/1aIfsVbZUbIwjDzlvOwWBA_lN2XTu2VGd/
    PENELOP 2018 - https://drive.google.com/file/d/1rb_wKkICOyL5UMuG4chuRxBQHuqR_8q1/
"""

__author__ = "Rui Campos"

print("Importing .materials.materials")


cimport cython 
from .logger import MaterialLogger


# THIS SHOULD BE AN IMPORT FROM .types 
class map(dict):
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError


from libc.stdlib cimport malloc, free

from .._init cimport LIMS
from .._init cimport EAX
from .._init import eax

cdef extern from "<math.h>" nogil:
    double frexp(double x, int* exponent)

from ..settings import __photonCUTOFF__, __electronCUTOFF__

from .photon.photon import Photon
from .electron.main cimport Electron

import numpy as np
import matplotlib.pyplot as plt

from . import database as db

def makeAlias(X, Y):
    X, Y = X.copy(), Y.copy()
    N = len(Y)
    Y = N*Y
    points = []
    while len(Y) > 0:
        
        i_min, i_max = np.argmin(Y), np.argmax(Y)
        
        
        ymax = Y[i_max]
        ymin = Y[i_min]
        
        xmin = X[i_min]
        xmax = X[i_max]
        
        dy = 1 - ymin
        
        Y[i_max] -= dy
        
        point = [xmin, ymin, xmax]
        points.append(point)
        
        X = np.delete(X, i_min)
        Y = np.delete(Y, i_min)

    
    points = sorted(points)
    return np.array(points)


def makeLinLin(x, y):
    m = np.diff(y)/np.diff(x)
    
    #y  = m*x - m*x[i] + y[i]
    #m*(x - x[i]) + y[i]
    return - m*x[:-1] + y[:-1], m




from ..settings import __montecarlo__
path = __montecarlo__/'materials'
path = str(path)

with open(path + "/compton_profiles/p-biggs.dat") as file:
    raw_grid = file.readlines()

gridAU = [float(x) for x in raw_grid]

import numpy as np
alpha = 1/137
gridAU = np.array(gridAU, dtype = float)*137
 ##################################### better value plzzz
cdef double[::1] CPGRID = gridAU

#print(np.array(CPGRID))


def rebuildShell(this):
    cdef Shell self
    self = <Shell>Shell.__new__(Shell)
    self.Z = this.Z
    self.index = this.index
    self.binding_energy = this.binding_energy;
    self.PHELa = this.PHELa
    self.PHELb = this.PHELb
    self.cCUMUL = this.cCUMUL
    self.cINVCUMUL = this.cINVCUMUL
    self.cumul = this.cumul

    return self


from scipy.interpolate import CubicSpline
cdef class Shell:
    def __reduce__(self):
        this = map()        
        this.Z = self.Z
        this.index = self.index
        this.binding_energy = self.binding_energy;
        this.PHELa = np.array(self.PHELa) 
        this.PHELb = np.array(self.PHELb)
        this.cCUMUL = np.array(self.cCUMUL)
        this.cINVCUMUL = np.array(self.cINVCUMUL)
        this.cumul = np.array(self.cumul)

        return rebuildShell, (this,)
    
    def __init__(self, formula,  Z, designator, int index, double N, double binding_energy): #most of the work is going to be here
        

        self.Z = Z
        self.index = index
        self.binding_energy = binding_energy;
        
        
        
        
        # PHOTOELECTRIC EFFECT
        ######################################################################
        TABLE = db.EPDL[Z-1][(7, 73, 91, designator, 0, 0)]
        X, Y = TABLE.X, TABLE.Y
        
        X = X*1e6 #changing to eV
        X = np.array(  list( dict.fromkeys(list(X)) ))
        Y = np.array(  list( dict.fromkeys(list(Y)) ))        
        spline = CubicSpline(X, N*Y*1e-24, extrapolate = False)
        Yeax = spline(eax)

        test = np.isinf(Yeax)
        if np.any(test):
            print("err")
            raise ValueError("materials.pyx :: .Shell: infinite values in one of the cross sections")
        
        np.nan_to_num(Yeax, nan=0, copy = False )

        self.PHELa, self.PHELb = makeLinLin(eax, Yeax)
        from numpy import array as arr
        formula.log.add_to_plot(formula.__fig, eax[:-1], arr(self.PHELa) + eax[:-1]*arr( self.PHELb), label = f"Shell #{index} | {designator} | Uk = {binding_energy}eV")

        ######################################################################
        

        
        # COMPTON PROFILE
        ######################################################################
        
        ######################################################################


        

    cdef double sample_compton_profile(self, mixmax_engine *genPTR, double pz_max):
        # print("")
        # print("starting compton profile sampling")
        # print("SHELL INDEX:", self.index)

        # print("pz_max", pz_max)
        
        cdef int START, END, MID
        
        START = 0
        END   = 29
        
        cdef double xMID
        
        #do binary search 
        while START <= END:
            #find middle
            MID = START + (END - START)//2 #prevents overflow somehow 
            
            xMID = CPGRID[MID]
            
            if pz_max == xMID: #found the value
                END = MID
                break
                #return MID
            
            if pz_max < xMID: # discard right side
                END = MID - 1 #do not include mid
                continue
            
            START = MID + 1
            


        # eval cumul 
        
        cdef double rc = self.cCUMUL[3, END]
        cdef double dx = pz_max - CPGRID[END]
        
  
        rc += self.cCUMUL[2, END]*dx
        
        dx *= dx
        rc += self.cCUMUL[1, END]*dx
        
        dx *= dx
        rc += self.cCUMUL[0, END]*dx
        
        rc = genPTR.get_next_float() * rc
        
        # eval invcumul
        START = 0
        END   = 29
        

        #do binary search 
        while START <= END:
            #find middle
            MID = START + (END - START)//2 #prevents overflow somehow 
            
            xMID = self.cumul[MID]
            
            if rc == xMID: #found the value
                END = MID
                break
                #return MID
            
            if rc < xMID: # discard right side
                END = MID - 1 #do not include mid
                continue
            
            START = MID + 1 
            
            

        
        dx = rc - self.cumul[END]
        rc = self.cINVCUMUL[3, END]
        
        rc += self.cINVCUMUL[2, END]*dx
        
        dx *= dx
        rc += self.cINVCUMUL[1, END]*dx
        
        dx *= dx
        rc += self.cINVCUMUL[0, END]*dx
        
        
    

        return rc
         
        
        
    def set_sampling_arrays(self, cCUMUL, cINVCUMUL, cumul):
        self.cCUMUL = cCUMUL
        self.cINVCUMUL = cINVCUMUL
        self.cumul = cumul

    def log_compton_profile(self, formula):
        formula.log.add_plot(self.cumul, self.cumul, title = f"Cumul", xlabel = "None", ylabel = "None")

        
    def construct_sampler(self, profile):
        from scipy.integrate import cumtrapz, trapz
        
        norm = trapz(profile, CPGRID)
        PDF = profile/norm
        
        # constructing cumul
        cumul = cumtrapz(PDF, CPGRID, initial = 0)
        
        spline = CubicSpline(CPGRID, cumul)
        self.cCUMUL = spline.c

        # constructing invCumul
        spline = CubicSpline(cumul, CPGRID)
        self.cINVCUMUL = spline.c
        self.cumul = cumul

        
    cdef void ionize(self, mixmax_engine *genPTR, PARTICLES* particles, double E):
        print(self.index)
        
        particles.ELECTRONS.push_back(E - self.binding_energy)
        
        print("RUNNING RELAX")
        self.rATOM.run(0, particles, genPTR)
        

    


def reconstruct_Atom(d):
    #cdef Atom self = Atom(d['Z'], d['CUT_OFF'], d['N'], pickle = True)    
    cdef Atom self
    self = <Atom>Atom.__new__(Atom)
    super(Atom, self ).__init__(d['Z'], d['CUT_OFF'])
    self.ALIAS = d['ALIAS']
    self.PHELa = d['PHELa'] 
    self.PHELb = d['PHELb'] 
    self.arrSHELLS = d['arrSHELLS'] 
    self.S = d['S']
    return self
    
#@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef class Atom(crAtom):
    
    def __reduce__(self):
         d = dict()
         d['Z'] = self.Z
         d['CUT_OFF'] = self.CUT_OFF
         d['ALIAS'] = np.array(self.ALIAS)
         d['PHELa'] = np.array(self.PHELa)
         d['PHELb'] = np.array(self.PHELb)
         d['arrSHELLS'] = np.array(self.arrSHELLS)
         d['S'] = self.S
         return reconstruct_Atom, (d,)
         
        
    
    def __init__(self, formula,  Z,  double CUT_OFF, double N, pickle = False):
        formula.log.add_header(f"Atom Z = {Z}", level = "h3")
        formula.log.add_paragraph("Input Information:")
        formula.log.add_attribute("Relaxation Model Cut Off: ", CUT_OFF)


        self.Z = Z
        self.CUT_OFF = CUT_OFF
        # RELAXATION MODEL
        ######################################################################
        super().__init__(Z, CUT_OFF)
        ######################################################################
        if pickle:
            return

        formula.log.add_attribute("Number of Shells: ", self.Nsh)

        
        
        # INITIALIZING SHELLS
        ######################################################################
        cdef int i
        temp = []
        cdef Shell shell


        formula.__fig = formula.log.new_plot()

        for i in range(self.Nsh):
            shell = Shell(formula, Z, self.DESIGNATORS[i], i, N, self.BE[i])
            #shell.rATOM = self.rATOM
            temp.append(shell)

        formula.log.finish_plot(formula.__fig, xlabel = "Energy (eV)", ylabel = "imfp (cm^-1)", logscale = True)
        del formula.__fig
         #, )

        self.arrSHELLS = np.array(temp)
        ######################################################################
        
        
        
        
        # PHOTOELECTRIC EFFECT 
        ######################################################################
        formula.log.add_header(f"PHOTOELECTRIC", level = "h4")

        PHELa = np.array(self.arrSHELLS[0].PHELa)
        PHELb = np.array(self.arrSHELLS[0].PHELb)
        
        for i in range(1, self.Nsh):
            PHELa += np.array(self.arrSHELLS[i].PHELa)
            PHELb += np.array(self.arrSHELLS[i].PHELb)
        
        self.PHELa = PHELa
        self.PHELb = PHELb

        from numpy import array as arr
        formula.log.add_plot(eax[:-1], arr(self.PHELa) + eax[:-1]*arr( self.PHELb), title = f"Atomic IMFP Photoelectric", xlabel = "Energy (eV)", ylabel = "imfp (cm^-1)")

        ######################################################################
        
        
        
        
        # GOS MODEL  
        ######################################################################
        
        ######################################################################
        
        
        
        
        # COMPTON  
        ######################################################################
        formula.log.add_header(f"COMPTON", level = "h4")

        # ATOM STUFF
        Stable = db.EPDL[Z-1][(7, 93, 0, 0, 0, 942)]
        h = 4.135667696e-15 #eV s
        m_e = 9.10938e-31 #kg
        c = 299792458 #m/s

        C = h/(m_e*c)
        CONST = C**2 * 1e-4
        
        
        self.S = CSa(CONST*Stable.X**2, Stable.Y/Z)
        
        
        probs = np.array(self.number_el)
        probs = probs/sum(probs)
        index = np.arange(0, len(probs))
        
        self.ALIAS = makeAlias(index, probs)

        formula.log.add_paragraph("REMINDER: add an histogram for prob of choosing atoms")
        
        
        # SHELL STUFF
        compton_profiles = self.get_profile(formula)
        from scipy.integrate import cumtrapz, trapz

        shell_indexes = [[0], # K
                         [1, 2, 3], # L
                         [4, 5, 6, 7, 8],# M
                         [9, 10, 11, 12, 13, 14, 15],# N
                         [16, 17, 18, 19, 20, 21, 22, 23, 24],# O
                         [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35],# P
                         [36, 37, 38] # Q
                         ]       
        
        
        formula.log.add_paragraph("> Constructing compton profile samplers")

        print(Z, self.Nsh)
        if Z < 36:
            subshell_index = 0
            formula.log.add_paragraph(f"     subshell_index = {subshell_index}")

            for shell_index, profile in enumerate(compton_profiles):
                formula.log.add_paragraph(f"            subshell_index = {subshell_index}")
                formula.log.add_paragraph(f"            shell_index = {shell_index}")

                subshells = shell_indexes[shell_index]
                
            
                norm = trapz(profile, CPGRID)
                PDF = profile/norm
                
                # constructing cumul
                cumul = cumtrapz(PDF, CPGRID, initial = 0)
                
                spline = CubicSpline(CPGRID, cumul)
                cCUMUL = spline.c
        
                # constructing invCumul
                spline = CubicSpline(cumul, CPGRID)
                cINVCUMUL = spline.c
                
                
                # shell = self.arrSHELLS[shell_index]
                # shell.set_sampling_arrays(cCUMUL, cINVCUMUL, cumul)
                
                
                
                for _ in range(len(subshells)):
                    
                    if subshell_index >= self.Nsh:
                        break
                    
                    shell = self.arrSHELLS[subshell_index]
                    shell.set_sampling_arrays(cCUMUL, cINVCUMUL, cumul)
                    shell.log_compton_profile(formula)
                    subshell_index += 1
        else:
                        
            for shell_index, profile in enumerate(compton_profiles):
                print(shell_index)
                
            
                norm = trapz(profile, CPGRID)
                PDF = profile/norm
                
                # constructing cumul
                cumul = cumtrapz(PDF, CPGRID, initial = 0)
                
                spline = CubicSpline(CPGRID, cumul)
                cCUMUL = spline.c
        
                # constructing invCumul
                spline = CubicSpline(cumul, CPGRID)
                cINVCUMUL = spline.c
                
                
                shell = self.arrSHELLS[shell_index]
                shell.set_sampling_arrays(cCUMUL, cINVCUMUL, cumul)
                shell.log_compton_profile(formula)
            if shell_index < self.Nsh-1:
                for si in range(shell_index, self.Nsh):
                    shell = self.arrSHELLS[si]
                    shell.set_sampling_arrays(cCUMUL, cINVCUMUL, cumul)
                    shell.log_compton_profile(formula)
            
            
            
            
            
            
            
            # print(self.Z, shell.index)

            # shell.construct_sampler(profile)
            
        # if self.Z == 8:
        #     shell = self.arrSHELLS[shell_index + 1]
        #     shell.construct_sampler(profile)
            
        ######################################################################
        
        
        
    def get_profile(self, formula):
        with open(path + f"/compton_profiles/profile-{self.Z}.dat") as file:
            raw_profile = file.readlines()
            raw_profile = [line.strip('\n') for line in raw_profile]
    
        profiles = []
    
        data = []
        for line in raw_profile:
            new_line = [float(x) for x in line.split()]
            data += new_line
    
        profiles = []
        profile = []
        
        profile = data[0:31]
        import numpy as np
        x = np.arange(0, len(profile))

        fig = formula.log.new_plot()

        #print("PROFILES OF", self.Z)
        for i in range(len(data)//31 ):
            print(i)
            profile = data[31*i:31*i+31]
            formula.log.add_to_plot(fig, x, profile)
            profiles.append(profile)
        formula.log.finish_plot(fig, title = "Compton Profiles", xlabel = "Bound Electron Momentum (a.u.)", ylabel = "Compton Profile")

        print("NUMBER OF PROFILES", len(profiles))
        import numpy as np
        return np.array(profiles)

        
        
        
    # cdef Shell choose_shell(self, mixmax_engine *genPTR, double E):
    #     """
    #     Randomly choose a shell based on occupation numbers.
    #     """
    #     if self.Nsh == 1: return self.arrSHELLS[0]
    #     cdef Shell shell
        
        
    #     cdef double R = self.Nsh*genPTR.get_next_float()
    #     cdef int N = <int> R
        
    #     if R - N < self.ALIAS[N, 1]:
    #         return self.ALIAS[<int> self.ALIAS[N, 0]]
    #     return self.arrATOMS[<int> self.ALIAS[N, 2]]
    
    
    cdef void PHELchoose(Atom self,int index, double E, mixmax_engine* genPTR, PARTICLES* particles):
        """
        Choose a shell based on photoelectric partial CS.
        """
        
        cdef double r = genPTR.get_next_float()*(self.PHELa[index] + E*self.PHELb[index])
        cdef int i
        cdef double cumul = 0
        cdef Shell shell
        for i in range(self.Nsh):
            cumul += self.arrSHELLS[i].PHELa[index] + self.arrSHELLS[i].PHELb[index] * E
            if r < cumul:
                shell = self.arrSHELLS[i]
                particles.ELECTRONS.push_back(E - shell.binding_energy)
                #print(genPTR.get_next_float())

               # print("RUNNING RELAX")
                self.run(shell.index, particles, genPTR)
                #shell.ionize(genPTR, particles, E)
                break
        else:
            raise RuntimeError("PHOTOELECTRIC DID NOT CHOOSE A SHELL")
            
    cdef void ionize(self, int shell_index, mixmax_engine* genPTR, PARTICLES* particles):
        self.run(shell_index, particles, genPTR)
        
        
        
def reconstruct_Molecule(d):
    #cdef Material self
    #self = <Material>Material.__new__(Material)

    cdef Molecule self = Molecule({}, 0, pickle = True)
    self.atomALIAS = d['atomALIAS']
    self.Nat = d['Nat'] 
    self.arrATOMS = d['arrATOMS'] 
    self.arrNi = d['arrNi'] 
    self.PHELa = d['PHELa'] 
    self.PHELb = d['PHELb'] 
    self.ALIAS =  d['ALIAS'] 
    self.Nsh = d['Nsh'] 
    return self




cdef class Molecule:   ################################## dont forget to choose where and how to multiply by N
    # cdef double Nat
    # cdef Atom* arrATOMS
    # cdef double* arrNi
    # cdef double[::1] PHELa, PHELb
    def __reduce__(self):
        d = dict()
        
        d['atomALIAS'] = np.array(self.atomALIAS)
        d['Nat'] = self.Nat
        d['arrATOMS'] = np.array(self.arrATOMS)
        d['arrNi'] = np.array(self.arrNi)
        d['PHELa'] = np.array(self.PHELa)
        d['PHELb'] = np.array(self.PHELb)
        d['ALIAS'] = np.array(self.ALIAS)
        d['Nsh'] = self.Nsh
        return reconstruct_Molecule, (d, )
        
    def __init__(self, formula, CUT_OFF, pickle = False):
        if pickle: return

        formula.log.add_header("Molecule", level = "h2")

        ### RANDOM IONIZATION
        self.ALIAS = formula.molecule_data.ALIAS
        self.Nsh = formula.molecule_data.Nsh
        
        
        ### PHOTOELECTRIC
        self.Nat = len(formula)
        formula.log.add_attribute("Number of Atoms", self.Nat)
        formula.log.add_attribute("Number of Shells", self.Nsh)
        #self.arrATOMS = <Atom*>malloc(self.Nat*sizeof(Atom))
        #self.arrNi = <double*>malloc(self.Nat*sizeof(double))
        
        cdef int i = 0
        
        temp = []
        probs = []
        ZZ = []
        arrNi = []
        for Z, Ni in formula.items():
            temp.append( Atom(formula, Z, CUT_OFF, formula.N))
            probs.append(Z*Ni)
            ZZ.append(Z)
            arrNi.append(formula[Z])
            #self.arrNi[i] = formula[Z] 
            i += 1
        self.arrATOMS = np.array(temp)
        self.arrNi = np.array(arrNi, dtype = float)
        probs = np.array(probs)
        probs = probs/sum(probs)
        
        PHELa = self.arrNi[0]*np.array(self.arrATOMS[0].PHELa)
        PHELb = self.arrNi[0]*np.array(self.arrATOMS[0].PHELb)
        
        for i in range(1, self.Nat):
            PHELa += self.arrNi[i]*np.array(self.arrATOMS[i].PHELa)
            PHELb += self.arrNi[i]*np.array(self.arrATOMS[i].PHELb)
        
        self.PHELa = PHELa
        self.PHELb = PHELb


        formula.log.add_paragraph("Finishing Molecule")
        from numpy import array as arr
        formula.log.add_plot(eax[:-1], arr(self.PHELa) + eax[:-1]*arr( self.PHELb), title = f"Molecular IMFP Photoelectric", xlabel = "Energy (eV)", ylabel = "imfp (cm^-1)")

        
        ## selecting atoms
        index = np.arange(0, len(probs))
        self.atomALIAS = np.array(makeAlias(index, probs), order = "F")


        
        
    cdef Atom get(self, double Z):
        cdef Atom atom
        for atom in self.arrATOMS:
            if atom.Z == Z:
                return atom
        # else:
        #     raise KeyError(f"Atom w/ Z = {Z} was not found in the molecule.")

    cdef void PHELionize(Molecule self, int index, double E, mixmax_engine *genPTR, PARTICLES* particles):
        
        
        
        cdef double r = genPTR.get_next_float()*( self.PHELa[index] + self.PHELb[index]*E)
        cdef double cumul = 0
        cdef int i
        for i in range(self.Nat):
            cumul += (self.arrATOMS[i].PHELa[index] + E*self.arrATOMS[i].PHELb[index])*self.arrNi[i]
            if r < cumul:
                (<Atom> self.arrATOMS[i]).PHELchoose(index, E, genPTR, particles)
                break
        else:
            raise RuntimeError("PHOTOELECTRIC  DID NOT CHOOSE AN ATOM")
   
    cdef void ionize(self, mixmax_engine *genPTR, PARTICLES* particles):
        
        cdef double r = genPTR.get_next_float()*self.Nsh
        cdef int N = <int> r
        
        if r - N < self.ALIAS[N, 0]: #accpet
            N = <int> self.ALIAS[N, 1]
            (<Atom> self.arrATOMS[N]).ionize(<int> self.ALIAS[N, 2], 
                                                      genPTR, 
                                                      particles)
        else:
            N = <int> self.ALIAS[N, 3]
            (<Atom> self.arrATOMS[N]).ionize(<int> self.ALIAS[N, 4], 
                                                      genPTR, 
                                                      particles)
            
            
    cdef void ionize_particular(self, mixmax_engine *genPTR, PARTICLES* particles, int atom_index, int shell_index):
            (<Atom> self.arrATOMS[atom_index]).ionize(shell_index, 
                                                      genPTR, 
                                                      particles)
            
    cdef Atom choose_atom(self, mixmax_engine *genPTR):
        cdef double R = self.Nat*genPTR.get_next_float()
        cdef int N = <int> R
        if R - N < self.atomALIAS[N, 1]:
            return self.arrATOMS[<int> self.atomALIAS[N, 0]]
        return self.arrATOMS[<int> self.atomALIAS[N, 2]]
        
        

        
        

#msg = "Can't modify initialized material!"



class dynamic_dict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

from .database import EADL, Na
def calculate(formula):
        
    
        Am = 0
        
        for Z in formula:
            Am += formula[Z]*EADL[Z-1]['Aw']
        
        
        N  = formula.density * Na / Am     
        
        return N, Am


# cdef class Molecule:

#     """
#     Container class for relaxation models and other information.
#     """
#     def __init__(self, formula):
#         self.formula = {}
#         self.rformula = {}
        
#         for Z in formula:
#             atom = Atom(Z)
#             self.formula[atom] = formula[Z]
#             self.rformula[Z] = atom
            
#         #self.formula = {Atom(Z):formula[Z]  for Z in formula}
#         #self.rformula = {formula[Z]:Atom(Z) for Z in formula}
        
#     def __getitem__(self, Z):
#         return self.rformula[Z]

cdef dict __CACHE__








def rebuildMaterial(this):

    cdef Material self
    self = <Material>Material.__new__(Material)
    self.name = this.name
    self.C1 = this.C1
    self.C2 = this.C2
    
    self.Wcc = this.Wcc
    self.Wcr = this.Wcr
    self.N = this.N
    self.Am = this.Am
    self.formula = this.formula

    self.molecule = this.molecule

    self.density = this.density
    
    self.photon   = this.photon
    self.electron = this.electron
    self.positron = this.positron
    return self



def Mat(formula, density, name = "Untitled", 
                 C1 = 0.1, C2 = 0.1,
                 Wcr = 10e3, Wcc = 100e3):
    """Create a new `Material` instance or read it from cache if it already has been compiled.
    """

    # Safety Check
    for Z in formula:
        if isinstance(formula[Z], float) or isinstance(formula[Z], int):
            raise ValueError(f"Coefficient of the element `{Z}` must be numeric type.")

        if isintance(Z, int):
            continue

        if isinstance(Z, float):
            formula[int(Z)] = formula.pop(Z)
            continue
		
        raise ValueError("Elements must be a numeric type.")
	
	# Check if it exists, return from cache if so		
    saved = f"{name}_C1{C1}_C2{C2}_Wcc{Wcc}_Wcr{Wcr}"

    the_path = "mat/" + saved
    from os import path
    if path.exists(the_path):
        import pickle
        with open(the_path, "rb") as f:
            mat = pickle.load(f)
        return mat

    # else, compile and return a new material
    return Material(formula, density, name = name, C1 = C1, C2 = C2, Wcr = Wcr, Wcc = Wcc)



cdef class Material:

    """
    Objects of this class are given to Vol.
    """

    def __reduce__(self):
        this = map()
        this.name = self.name
        this.C1 = self.C1
        this.C2 = self.C2
        
        this.Wcc = self.Wcc
        this.Wcr = self.Wcr
        this.N = self.N
        this.Am = self.Am
        this.formula = self.formula
        this.density = self.density

        this.molecule = self.molecule
#
        
        #
        this.photon   = self.photon
        this.electron = self.electron
        this.positron = self.positron
        return rebuildMaterial, (this,)



    def __init__(self, formula, density, name = "Untitled", 
                 C1 = 0.1, C2 = 0.1,
                 Wcr = 10e3, Wcc = 100e3):
        


        self.name = name

        print(f"Compiling data for material '{name}'")

        self.formula = formula
        
        formula = dynamic_dict(formula)
        formula.log = MaterialLogger(formula, density, name = name, C1 = C1, C2 = C2, Wcr = Wcr, Wcc = Wcc)

        formula.log.add_paragraph("> Setting general information into 'formula' dynamic_dict")
        formula.molecule_data = db.MoleculeDATA(formula)
        
        formula.density = density
        formula.C1 = C1
        formula.C2 = C2
        self.C1 = C1
        self.C2 = C2
        
        self.Wcc = Wcc
        self.Wcr = Wcr
        
        formula.Wcr = Wcr 
        formula.Wcc = Wcc
        
        N, Am = calculate(formula)
        formula.log.add_paragraph(f"Number Density = {N}")
        formula.log.add_paragraph(f"Atomic Weight = {Am}")

        formula.N = N
        formula.Am = Am


        self.N = N
        self.Am = Am   
        
        formula.log.add_paragraph("> Constructing Molecule -> *Atom -> *Shell structure and preparing relaxation models.")
        self.molecule = Molecule(formula, min(__photonCUTOFF__, __electronCUTOFF__))
       # self.photon.molecule = self.molecule
        formula.relax = self.molecule
        formula.molecule = self.molecule

        self.density = density
        
 

        formula.log.add_header("COMPILING PHOTON DATA")
        self.photon   = Photon(formula, density)
        

        
   
        
        
       
        formula.log.add_header("COMPILING ELECTRON DATA")

        self.electron = Electron(formula)
        
        formula.log.add_header("COMPILING POSITRON DATA")

        self.positron = Positron(formula)
        

        
        formula.log.add_paragraph("> Writing to html file.")
        formula.log.write(self.name + ".html")

        import pickle
        to_save = f"{name}_C1{C1}_C2{C2}_Wcc{Wcc}_Wcr{Wcr}"
        with open(f"mat/{to_save}", "wb") as f:
            pickle.dump(self, f)

        

        
    
    def __repr__(self):
        rep =  f"<{self.name}: "
        for Z, x in self.formula.items():
            rep += f"{x}x(Z = {Z}) | "
        rep += f"{self.density} g/cm^3  | "
        rep += f"C1 = {self.C1} C2 = {self.C2} Wcc = {self.Wcc} Wcr = {self.Wcr} >"
        
        return rep
    
    def getEl(self):
        return self.electron


    def plot_photonsIMFP(self):
        import matplotlib.pyplot as plt
        import numpy as np
        arr = np.array
        # coherent
        
        
        
     #   plt.plot(eax[:-1], totalCS[0] + totalCS[1]*eax[:-1])

        plt.plot(eax[:-1], arr(self.photon.coherent.imfpA) + eax[:-1]*arr( self.photon.coherent.imfpB), label = "coh")
        plt.plot(eax[:-1], arr(self.photon.incoherent.imfpA) + eax[:-1]*arr( self.photon.incoherent.imfpB), label = "incoh")
        plt.plot(eax[:-1], arr(self.photon.pairproduction.imfpA) + eax[:-1]*arr( self.photon.pairproduction.imfpB), label = "pair")
        plt.plot(eax[:-1], arr(self.photon.tripletproduction.imfpA) + eax[:-1]*arr( self.photon.tripletproduction.imfpB), label = "trip")
        plt.plot(eax[:-1], arr(self.molecule.PHELa) + eax[:-1]*arr( self.molecule.PHELb), label = "photoelectric")
        plt.legend()
        plt.xscale("log")
        plt.xlabel("Energy (eV)")
        plt.ylabel("Inverse Mean Free Path (cm^-1)")
        plt.yscale("log")
        plt.show()

 

    def plot_photonsMFP(self):
        import matplotlib.pyplot as plt
        import numpy as np
        arr = np.array
        # coherent
        
        
        
     #   plt.plot(eax[:-1], totalCS[0] + totalCS[1]*eax[:-1])
        imfp = arr(self.photon.coherent.imfpA) + eax[:-1]*arr( self.photon.coherent.imfpB)
        plt.plot(eax[:-1], 1/imfp, label = "coh")

        imfp = arr(self.photon.incoherent.imfpA) + eax[:-1]*arr( self.photon.incoherent.imfpB)
        plt.plot(eax[:-1], 1/imfp, label = "incoh")

        imfp = arr(self.photon.pairproduction.imfpA) + eax[:-1]*arr( self.photon.pairproduction.imfpB)
        select = imfp != 0
        imfp = imfp[select]
        x = (eax[:-1])[select]
        plt.plot(x, 1/imfp, label = "pair")

        imfp = arr(self.photon.tripletproduction.imfpA) + eax[:-1]*arr( self.photon.tripletproduction.imfpB)
        select = imfp != 0
        imfp = imfp[select]
        x = (eax[:-1])[select]
        plt.plot(x, 1/imfp, label = "trip")



        imfp = arr(self.molecule.PHELa) + eax[:-1]*arr( self.molecule.PHELb)
        select = imfp != 0
        imfp = imfp[select]
        x = (eax[:-1])[select]
        plt.plot(eax[:-1], 1/imfp, label = "photoelectric")


        plt.legend()
        plt.xscale("log")
        plt.xlabel("Energy (eV)")
        plt.ylabel("Mean Free Path (cm)")
        plt.yscale("log")
        plt.show()

        
    def EXTRACT_ALL(self):
        return {"el": self.electron.EXTRACT_ALL()}
    
    
    
    #Eventually I should find a good naming convention
    @property
    def el(self): return self.electron
    
    @property
    def po(self): return self.positron
	
    @property
    def ph(self): return self.photon
    
    @property
    def mol(self): return self.molecule
        

    
        
