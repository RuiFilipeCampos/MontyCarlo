#distutils: language = c++
# cython: profile=True
# cython: annotate=True

from libc.math cimport log, sqrt
cimport cython

cdef double ELECTRON_REST_ENERGY = 0.51099895000*1e6 #eV
cdef double  _2ELECTRON_REST_ENERGY    = 2 *ELECTRON_REST_ENERGY

from libcpp.vector cimport vector
from ..._init import eax

from .GOS cimport CMolecule
from .GOS cimport CAtom
from .GOS cimport CShell

class MAP(dict):
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

import numpy as np


cdef struct LinLin:
    double *a
    double *b
    
# cdef LinLin _makeLinLin(double[:] x, double[:] y):
#     cdef LinLin interp
#     res = makeLinLin(np.array(x), np.array(y))
#     cdef double a = res[0]
#     cdef double b = res[1]
    
#     interp.a = &a[0]
#     interp.b = &b[0]
#     return interp
   
def makeLinLin(x, y):
    m = np.diff(y)/np.diff(x)
     
    #y  = m*x - m*x[i] + y[i]
    #m*(x - x[i]) + y[i]
    
    res = np.array([- m*x[:-1] + y[:-1], m])
    
    return res

getLinLin = makeLinLin



def rebuildgosMolecule(this):
    cdef gosMolecule self
    self = <gosMolecule> gosMolecule.__new__(gosMolecule)
    self.Nat = this.Nat

    self.totalCS = this.totalCS
    self.number_density = this.number_density

    self.gosATOMS = this.gosATOMS
    return self






cdef class gosMolecule:
    
    
    def __reduce__(self):
        this = MAP()
        this.Nat = self.Nat
        import numpy as np

        this.totalCS = np.array(self.totalCS)
        this.number_density = self.number_density


        this.gosATOMS = np.array(self.gosATOMS)
        return rebuildgosMolecule, (this, )


    def __init__(gosMolecule self, CMolecule cmolecule, formula):
        self.Nat = cmolecule.N + 1 #to include CB shell
        
        print("initializing memview")
        self.totalCS = makeLinLin(eax, cmolecule.SIGMA0)
        print("done")
        gosATOMS = []
        cdef CAtom catom
        
        
        
        
        self.number_density = formula.N
        
        #cb shell first, most likely one
        #catom = cmolecule.cb   #cb inherits from CAtom
        
        
        
        
        gosATOMS.append(gosCBShell(cmolecule.cb, cmolecule.SIGMA0cb))
        
        
        for catom in cmolecule.ATOMS:
            gosATOMS.append(gosAtom(catom, formula))
            
        self.gosATOMS = np.array(gosATOMS)
        
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef void sample(gosMolecule self, mixmax_engine* genPTR, int index, double E, PARTICLES *particles):
        
        
        cdef double r = genPTR.get_next_float()*(self.totalCS[0, index] + E*self.totalCS[1, index])
        
        
        cdef gosAtom atom
        cdef int i = 0
        cdef double cumul = 0
        for i in range(self.Nat):
            atom = self.gosATOMS[i]
            cumul += atom.totalCS[0, index] + E*atom.totalCS[1, index]
            if r < cumul:
                atom.sample(genPTR, index, E, particles)
                
                break
        
    
    
    def cut(self, double Wcc, shells_for_softIMFP):
        """
        
        """
        newATOMS = []
        currentATOMS = np.array(self.gosATOMS)
        cdef gosAtom atom
        
        atom = self.gosATOMS[0]
        if (<gosCBShell> atom).Wk > Wcc: 
            return False

        shells_for_softIMFP.append([(<gosCBShell> atom).fk, (<gosCBShell> atom).Wk])
        
        
            
        cdef int i = 0
        for atom in currentATOMS[1:]:
            empty = atom.cut(Wcc, shells_for_softIMFP)
            if not empty: 
                i += 1
                newATOMS.append(atom)


        
        if i == 0:
            self.totalCS = np.array(self.totalCS)*0
            return True
        self.Nat = len(newATOMS)
        self.gosATOMS = np.array(newATOMS)
        
        new_totalCS = np.array(self.totalCS)*0
        for atom in self.gosATOMS:
            new_totalCS += np.array(atom.totalCS)
            
        self.totalCS = np.array(new_totalCS) #- np.array(atom.totalCS)
        print(np.array(self.totalCS))
        
        
        return False
    
    
    
    def makearray(self):
        #cdef double [:, ::1]
 
        
        cdef double a, b, Uk, Wk
        arr = [list() for _ in range((len(eax)-1))]
        arr_atoms = list()
        
        cdef gosAtom atom
        cdef gosShell shell
     #   ionizationCS[identify which shell, A OR B, ENERGY_index]
        cdef int lenposs = 0
        cdef int atom_index, shell_index
        for atom_index in range(self.Nat):
            
            print("molec.Nat", self.Nat)
            atom = self.gosATOMS[atom_index]
            arr_atoms.append(atom.ATOMptr)
            for shell_index in range(atom.Nsh):
                print("atom.Nsh", atom.Nsh)
                shell = atom.gosSHELLS[shell_index]
                Wk = shell.Wk
                
                if shell.Nsh == 0:
                    for i in range(len(eax) - 1):
                        CS = np.array(shell.totalCS)
                        
                        CS_E = CS[:, i]
    
                        a = CS_E[0]*self.number_density
                        b = CS_E[1]*self.number_density
                        p1, p2, p3 = shell.p1, shell.p2, shell.p3

                        arr[i] += [a*p1[i], b*p1[i], Wk, 0, 0., 0.]
                        arr[i] += [a*p2[i], b*p2[i], Wk, 0, 0., 0.]
                        arr[i] += [a*p3[i], b*p3[i], Wk, 0, 0., 0.]
                        continue
                
                for j in range(shell.Nsh):
                    print("shell.Nsh", shell.Nsh)
                    CS = np.array(shell.ionizationCS[j])
                    Uk = shell.BE[j]
                    lenposs += 6*3
                    p1, p2, p3 = shell.p1, shell.p2, shell.p3
                    print(Uk)
                    for i in range(len(eax) - 1):
                        CS_E = CS[:, i]
                        a = CS_E[0]*self.number_density
                        b = CS_E[1]*self.number_density
                        
                        arr[i] += [a*p1[i], b*p1[i], Wk, Uk, <double> atom_index, <double> shell.INDEX[j]]
                        arr[i] += [a*p2[i], b*p2[i], Wk, Uk, <double> atom_index, <double> shell.INDEX[j]]
                        arr[i] += [a*p3[i], b*p3[i], Wk, Uk, <double> atom_index, <double> shell.INDEX[j]]
                        
                        
                    # for j in range(shell.Nsh):
                    #     a = shell.ionizationCS[j, 0, i]
                    #     b = shell.ionizationCS[j, 1, i]
                    #     Uk = shell.BE[j]
                    #     arr[i] += [a, b, Wk, Uk]
        
        arr = np.array(arr, order = "C")
        arr_atoms = np.array(arr_atoms, order = "C")
        return arr, len(arr), arr_atoms
    
    
    def __getitem__(self, Z):
        cdef gosAtom atom
        for atom in self.gosATOMS:
            if atom.Z == Z:
                return atom
        
    def plot( self, units = "cm2"):
        import matplotlib.pyplot as plt
        import numpy as np
        
        
        fig = plt.figure()
        totalCS = np.array(self.totalCS)
        if units == "cm2":
            
            plt.plot(eax[:-1], totalCS[0] + totalCS[1]*eax[:-1])
        
            plt.ylabel("MOLECULE (cm^2) ")
        elif units == "barn":
            plt.plot(eax[:-1], (totalCS[0] + totalCS[1]*eax[:-1])*1e24)
            plt.ylabel("MOLECULE (barn) ")
            
            
            
        plt.xlabel("Energy of Incident Electron (eV)")
        plt.xscale("log"); plt.yscale("log");
        plt.grid(which = "both")
        plt.show()
            
            
        
def rebuildgosAtom(this):
    cdef gosAtom self
    self = <gosAtom> gosAtom.__new__(gosAtom)
    self.Z = this.Z
    self.Nsh = this.Nsh

    self.totalCS = this.totalCS
    self.gosSHELLS = this.gosSHELLS
    return self
            
            
cdef class gosAtom:
    #cdef gosShell[::1] gosSHELLS
    #cdef double[:, ::1] totalCS
    


    def __reduce__(self):
        this = MAP()
        this.Z = self.Z
        this.Nsh = self.Nsh

        from numpy import array
        this.totalCS = array(self.totalCS)
        this.gosSHELLS = array(self.gosSHELLS)

        return rebuildgosAtom, (this, )


    def __init__(gosAtom self, CAtom catom, formula):
        self.Z = catom.Z
        self.Nsh = catom.N
        
        cdef Molecule molecule = <Molecule> formula.molecule
        
        #self.Nat = molecule.Nat
        #self.rATOM = molecule.get(self.Z)
        cdef int ii 
        for ii in range(molecule.Nat):
            if (<Atom> molecule.arrATOMS[ii]).Z == self.Z:
                self.rATOM =  (<Atom> molecule.arrATOMS[ii]).rATOM
                self.ATOMptr =molecule.arrATOMS[ii]
                #print((<Atom> molecule.arrATOMS[ii]).Z)
                break
        else: 
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>NO ATOMS OFUND")
            
            raise RuntimeError("NO ATOM FOUND")
        
        
        
        totalCS = np.array(catom.SIGMA0) #needed for determining normalzation procedure
        self.totalCS = makeLinLin(eax, totalCS)
        
        
        
        # initing shells
        
        gosSHELLS = []
        cdef CShell cshell
        
        #outer shell first -> most likely one
        cshell = catom.SHELLS[-1]
        gosSHELLS.append(gosOuterShell(cshell, catom.Z))
        
        
        
        
        inner_totalCS = np.zeros(len(totalCS))
        
        cdef gosInnerShell gos_shell
        
        for cshell in catom.SHELLS[0:-1]:
            gos_shell = gosInnerShell(cshell, catom.Z, inner_totalCS)
            gos_shell.rATOM = self.rATOM
            gos_shell.ATOMptr = self.ATOMptr
            gosSHELLS.append(gos_shell)
        
        
        
        # RENORMALIZATION 
       # cdef LinLin outer_totalCS = _makeLinLin(eax, totalCS - inner_totalCS)

        (<gosOuterShell> gosSHELLS[0]).set_totalCS(makeLinLin(eax, totalCS - inner_totalCS))

        self.gosSHELLS = np.array(gosSHELLS)
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True)   
    cdef void sample(gosAtom self, mixmax_engine* genPTR, int index, double E, PARTICLES *particles):

        cdef double r = genPTR.get_next_float()*(self.totalCS[0, index] + E*self.totalCS[1, index])
        
        
        cdef gosShell shell
        cdef int i = 0
        cdef double cumul = 0
        for i in range(self.Nsh):
            
            shell = self.gosSHELLS[i]
            cumul += shell.totalCS[0, index] + E*shell.totalCS[1, index]
            if r < cumul:
                #if shell.Wk < E: (22424.126634367873, 22325.69824366879)
                shell.sample(genPTR, index, E, particles)
                break
                #else: print(shell.Wk, E) (22:424.126634367873, 22:325.69824366879)
                
        #print(r, cumul)    (1.117703167661862e-22, 5.158809176477914e-22)
            
            
    def __getitem__(self, shell_index):
        return self.gosSHELLS[shell_index]
        
    def plot(gosAtom self, units = "cm2"):
        import matplotlib.pyplot as plt
        import numpy as np
        
        
        fig = plt.figure()
        totalCS = np.array(self.totalCS)
        if units == "cm2":
            
            plt.plot(eax[:-1], totalCS[0] + totalCS[1]*eax[:-1])
        
            plt.ylabel("ATOM (cm^2) ")
        elif units == "barn":
            plt.plot(eax[:-1], (totalCS[0] + totalCS[1]*eax[:-1])*1e24)
            plt.ylabel("ATOM (barn) ")
            
            
            
        plt.xlabel("Energy of Incident Electron (eV)")
        plt.xscale("log"); plt.yscale("log");
        plt.grid(which = "both")
        plt.show()
    
    def cut(self, double Wcc, shells_for_softIMFP):
        cdef gosShell shell
        newSHELLS = []
        
        
   
        
        for shell in self.gosSHELLS:
            if shell.Wk < Wcc:
                shells_for_softIMFP.append([shell.fk, shell.Wk])
                continue
            newSHELLS.append(shell)
            
        if len(newSHELLS) == 0: return True
        self.Nsh = len(newSHELLS)
        self.gosSHELLS = np.array(newSHELLS)
        
         #- np.array(shell.totalCS)
        new_totalCS = 0*np.array(self.totalCS)
        
        for shell in self.gosSHELLS:
            new_totalCS += np.array(shell.totalCS)
        self.totalCS = np.array(new_totalCS)
        
        return False
        
        
from scipy.interpolate import CubicSpline

cdef class gosShell:
    # cdef double Wk;
    # cdef:
    #     double[:, ::1] totalCS;
    #     double[:, ::1] pCLOSE 
    #     double[:, ::1] pLFAR  
    #     double[:, ::1] pTFAR  
    
    #cdef double[::1] sigma0A, sigma0B
    
    def plot(gosShell self, units = "cm2"):
        import matplotlib.pyplot as plt
        import numpy as np
        
        
        fig = plt.figure()
        totalCS = np.array(self.totalCS)
        if units == "cm2":
            
            plt.plot(eax[:-1], totalCS[0] + totalCS[1]*eax[:-1])
        
            plt.ylabel("Shell Cross Section (cm^2) ")
        elif units == "barn":
            plt.plot(eax[:-1], (totalCS[0] + totalCS[1]*eax[:-1])*1e24)
            plt.ylabel("Shell Cross Section (barn) ")
            
            
            
        plt.xlabel("Energy of Incident Electron (eV)")
        plt.xscale("log"); plt.yscale("log");
        plt.grid(which = "both")
        plt.show()
        
        
    cdef void sample(self, mixmax_engine* genPTR, int index, double E, PARTICLES *particles):
        print("GOS :: sample method was not overriden in gosShell")
        raise RuntimeError("GOS :: sample method was not overriden in gosShell")



def rebuildgosOuterShell(this):
    cdef gosOuterShell self
    self = <gosOuterShell> gosOuterShell.__new__(gosOuterShell)


    self.Wk = this.Wk    
    self.fk = this.fk   
    self.Nsh = this.Nsh
    from numpy import array

    self.p1     = array(this.p1    )  
    self.p2     = array(this.p2    )  
    self.p3     = array(this.p3    )   
    
    self.pCLOSE = array(this.pCLOSE, order = "F")  #  = np.array(makeLinLin(eax, CLOSE/norm), order = "F")
    self.pLFAR  = array(this.pLFAR , order = "F")  #   = np.array(makeLinLin(eax, LFAR/norm), order = "F")
    self.pTFAR  = array(this.pTFAR , order = "F")  #  = np.array(makeLinLin(eax, TFAR/norm), order = "F")

    self.totalCS = array(this.totalCS, order = "F")  #np.array(interp, order = "F")
    return self
       
cdef class gosOuterShell(gosShell):


    def __reduce__(self):
        this = MAP()
        this.Wk = self.Wk    
        this.fk = self.fk   
        this.Nsh    = self.Nsh
        
        from numpy import array
        this.p1     = array(self.p1    )  
        this.p2     = array(self.p2    )  
        this.p3     = array(self.p3    )   
        
        this.pCLOSE = array(self.pCLOSE)  #  = np.array(makeLinLin(eax, CLOSE/norm), order = "F")
        this.pLFAR  = array(self.pLFAR )  #   = np.array(makeLinLin(eax, LFAR/norm), order = "F")
        this.pTFAR  = array(self.pTFAR )  #  = np.array(makeLinLin(eax, TFAR/norm), order = "F")

        this.totalCS = array(self.totalCS)  #np.array(interp, order = "F")
        return rebuildgosOuterShell, (this,)


    def __init__(self, CShell cshell, double Z):
        self.Wk = cshell.Wk
        self.fk = cshell.fk
        self.Nsh = 0
        #self.totalCS = makeLinLin(eax, totalCS)

        #these guys already come evaluated at eax
        CLOSE = np.array(cshell.close.SIGMA0)
        LFAR  = np.array(cshell.distant.SIGMA0L)
        TFAR  = np.array(cshell.distant.SIGMA0T)
         
        norm = CLOSE + LFAR + TFAR
        norm[norm == 0] = 1
        
        
        # self.pCLOSE = makeLinLin(eax, CLOSE/norm)
        # self.pLFAR  = makeLinLin(eax, LFAR/norm)
        # self.pTFAR  = makeLinLin(eax, TFAR/norm)
        
        self.p1 = CLOSE/norm
        self.p2 = LFAR/norm
        self.p3 = TFAR/norm
        
        self.pCLOSE = np.array(makeLinLin(eax, CLOSE/norm), order = "F")
        self.pLFAR  = np.array(makeLinLin(eax, LFAR/norm), order = "F")
        self.pTFAR  = np.array(makeLinLin(eax, TFAR/norm), order = "F")
        
        
        
        
    def set_totalCS(self, double[:, ::1] interp):
        self.totalCS = np.array(interp, order = "F")
        
    def plot(gosShell self, units = "cm2"):
        import matplotlib.pyplot as plt
        import numpy as np
        
        
        fig = plt.figure()
        totalCS = np.array(self.totalCS)
        if units == "cm2":
            
            plt.plot(eax[:-1], totalCS[0] + totalCS[1]*eax[:-1])
        
            plt.ylabel("Shell Cross Section (cm^2) ")
        elif units == "barn":
            plt.plot(eax[:-1], (totalCS[0] + totalCS[1]*eax[:-1])*1e24)
            plt.ylabel("Shell Cross Section (barn) ")
            
            
            
        plt.xlabel("Energy of Incident Electron (eV)")
        plt.xscale("log"); plt.yscale("log");
        plt.grid(which = "both")
        plt.show()
        
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True)   
    cdef void sample( self, mixmax_engine* genPTR, int index, double E, PARTICLES *particles):

        r = genPTR.get_next_float()
        
        
        cumul = self.pCLOSE[0, index] + E*self.pCLOSE[ 1, index]
        

        
  
            
        if r < cumul: 
            self.sampleCLOSE(E, genPTR, particles)
            return
        
        #cumul += self.pLFAR[index, 0] + E*self.pLFAR[index, 1]
            
        if r < cumul + self.pLFAR[0, index] + E*self.pLFAR[1, index]:
            self.sampleLFAR(E, genPTR, particles)
            return
        
        # new energy and cos for projectile
        particles.ELECTRONS.push_back(E - self.Wk)
        particles.ELECTRONS.push_back(1.)
        
        # energy and cos for secondary electron
        particles.ELECTRONS.push_back(self.Wk)
        particles.ELECTRONS.push_back(1.)

        #cdef double Esec = self.Wk - self.BE[i]
        #cdef double newE = E - self.Wk
        return #(newE, 1., Esec, -1.)
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef void sampleCLOSE(self, double E, mixmax_engine *genPTR, PARTICLES *particles):
        
        cdef double r
        cdef double kc = self.Wk/E
        cdef double a = (E/(E + ELECTRON_REST_ENERGY))**2
        cdef double CONST = 1+ 5*a*kc
        cdef double CONST2 = (1 - 2*kc)
        cdef double CONST3 = CONST2/(5*a*kc)
    
        while 1:
            r = CONST*genPTR.get_next_float()
            
            if r < 1: r = kc / (1 - r*CONST2)
            else: r = kc + (r - 1)*CONST3
            
            if r < kc: continue
            if r > .5: continue
        
         #   else: P = (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(k*(1-k))))
                
            if genPTR.get_next_float()*(1+5*a*r*3) < r*r*  (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(r*(1-r))  )):
                break
        
        
        # EFFECT ON TARGET
        cdef double W = r*E
        cdef double cos2 = (E - W)*(E + _2ELECTRON_REST_ENERGY)/E/(E - W + _2ELECTRON_REST_ENERGY)
        particles.ELECTRONS.push_back(E - W)
        particles.ELECTRONS.push_back(sqrt(cos2))
        
        # SECONDARY ELECTRON
        cos2 = (W)*(E + _2ELECTRON_REST_ENERGY)/E/(W + _2ELECTRON_REST_ENERGY)
        particles.ELECTRONS.push_back(W)
        particles.ELECTRONS.push_back(sqrt(cos2))

        
        
        #cos2 = W * (E + _2ELECTRON_REST_ENERGY)/E/(W + _2ELECTRON_REST_ENERGY)
        #cdef double cos_sec = sqrt(cos2)
        
        
        
        
        
        
                # new energy and cos for projectile
                
                

        
        # # energy and cos for secondary electron
        # particles.ELECTRONS.push_back(self.Wk - self.BE[i])
        # particles.ELECTRONS.push_back(cos_sec)
        
        
        # if self.relax:
        #     return (STATE.E - W, cos, W-self.Uk, cos_sec)         
        # return (STATE.E - W, cos, W, cos_sec)
    
        
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True)   
    cdef double Q_minus(self, double E):
        """eq 3.83"""
        cdef double dE = E - self.Wk
        return sqrt( (self.momentum(E) - self.momentum(dE))**2  + ELECTRON_REST_ENERGY**2  ) \
               - ELECTRON_REST_ENERGY
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef double momentum(self, double E):
        return sqrt(  E * (E + _2ELECTRON_REST_ENERGY)  ) 
    
    
    
            #cdef double newE = STATE.E - self.Wk 
       # cdef double Esec = self.Wk - self.Uk if self.relax else self.Wk
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef void sampleLFAR(self, double E, mixmax_engine *genPTR, PARTICLES *particles):
        
        
        cdef double Qm = self.Q_minus(E) #
       # if Qm / _2ELECTRON_REST_ENERGY < 0:
        #    print(E, Qm, _2ELECTRON_REST_ENERGY, Qm / _2ELECTRON_REST_ENERGY, self.Wk)
        
        cdef double Qs = Qm/(1 + Qm / _2ELECTRON_REST_ENERGY)

         
        #print(_2ELECTRON_REST_ENERGY)
        

        cdef double A = Qs/self.Wk *(1 + self.Wk/_2ELECTRON_REST_ENERGY )
        cdef double B = Qs/_2ELECTRON_REST_ENERGY
        cdef double Q = Qs / (A**(genPTR.get_next_float()) - B)

        cdef double p2E = E * (E + _2ELECTRON_REST_ENERGY)
        cdef double p2d = (E - self.Wk) * ((E - self.Wk) + _2ELECTRON_REST_ENERGY)   #self.p2(E - self.Wk)
        cdef double p2Q = Q * (Q + _2ELECTRON_REST_ENERGY)

        cdef double cos = p2E + p2d - p2Q
        #cos = .5*cos/(p2E*p2d)**.5
        particles.ELECTRONS.push_back(E - self.Wk)
        particles.ELECTRONS.push_back(.5*cos/(p2E*p2d)**.5)
        
        
        
        
        
        cdef double Wk2 = self.Wk*self.Wk
        cdef double X = (E/ELECTRON_REST_ENERGY + 1)**2
        A = Wk2 / ((X-1)/X) / p2Q

        B = (p2Q - Wk2)/(2*self.Wk*(E + _2ELECTRON_REST_ENERGY))
        
        particles.ELECTRONS.push_back(self.Wk)
        particles.ELECTRONS.push_back(A * (1 + B)**2)




def rebuildgosCBShell(this):
    cdef gosCBShell self
    self = <gosCBShell> gosCBShell.__new__(gosCBShell)


    self.Wk = this.Wk    
    self.fk = this.fk   
    from numpy import array

  
    
    self.pCLOSE = this.pCLOSE  #  = np.array(makeLinLin(eax, CLOSE/norm), order = "F")
    self.pLFAR  = this.pLFAR   #   = np.array(makeLinLin(eax, LFAR/norm), order = "F")
    self.pTFAR  = this.pTFAR   #  = np.array(makeLinLin(eax, TFAR/norm), order = "F")

    self.totalCS = this.totalCS  #np.array(interp, order = "F")

    return self



cdef class gosCBShell(gosAtom):
    cdef double Wk;
    cdef double fk
    cdef:
        double[:, ::1] pCLOSE 
        double[:, ::1] pLFAR  
        double[:, ::1] pTFAR  
    
    
    def __reduce__(self):
        this = MAP()
        this.Wk = self.Wk    
        this.fk = self.fk   
         
        from numpy import array
        this.pCLOSE = array(self.pCLOSE)  #  = np.array(makeLinLin(eax, CLOSE/norm), order = "F")
        this.pLFAR  = array(self.pLFAR )  #   = np.array(makeLinLin(eax, LFAR/norm), order = "F")
        this.pTFAR  = array(self.pTFAR )  #  = np.array(makeLinLin(eax, TFAR/norm), order = "F")

        this.totalCS = array(self.totalCS)  #np.array(interp, order = "F")
        return rebuildgosCBShell, (this,)


    def __init__(self, CShell cbshell, totalCS):
        self.Wk = cbshell.Wk
        self.fk = cbshell.fk
        
        self.totalCS = np.array(makeLinLin(eax, cbshell.SIGMA0))

        #these guys already come evaluated at eax
        CLOSE = np.array(cbshell.close.SIGMA0)
        LFAR  = np.array(cbshell.distant.SIGMA0L)
        TFAR  = np.array(cbshell.distant.SIGMA0T)
         
        norm = CLOSE + LFAR + TFAR
        norm[norm == 0] = 1
        
        
        self.pCLOSE = makeLinLin(eax, CLOSE/norm)
        self.pLFAR  = makeLinLin(eax, LFAR/norm)
        self.pTFAR  = makeLinLin(eax, TFAR/norm)
    
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef void sample( self, mixmax_engine* genPTR, int index, double E, PARTICLES *particles):

        r = genPTR.get_next_float()
        
        
        cumul = self.pCLOSE[0, index] + E*self.pCLOSE[ 1, index]
        

        
  
            
        if r < cumul: 
            self.sampleCLOSE(E, genPTR, particles)
            return
        
        #cumul += self.pLFAR[index, 0] + E*self.pLFAR[index, 1]
            
        if r < cumul + self.pLFAR[0, index] + E*self.pLFAR[1, index]:
            self.sampleLFAR(E, genPTR, particles)
            return
        
        # new energy and cos for projectile
        particles.ELECTRONS.push_back(E - self.Wk)
        particles.ELECTRONS.push_back(1.)
        
        # energy and cos for secondary electron
        particles.ELECTRONS.push_back(self.Wk)
        particles.ELECTRONS.push_back(1.)

        #cdef double Esec = self.Wk - self.BE[i]
        #cdef double newE = E - self.Wk
        return #(newE, 1., Esec, -1.)
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef void sampleCLOSE(self, double E, mixmax_engine *genPTR, PARTICLES *particles):
        
        cdef double r
        cdef double kc = self.Wk/E
        cdef double a = (E/(E + ELECTRON_REST_ENERGY))**2
        cdef double CONST = 1+ 5*a*kc
        cdef double CONST2 = (1 - 2*kc)
        cdef double CONST3 = CONST2/(5*a*kc)
    
        while 1:
            r = CONST*genPTR.get_next_float()
            
            if r < 1: r = kc / (1 - r*CONST2)
            else: r = kc + (r - 1)*CONST3
            
            if r < kc: continue
            if r > .5: continue
        
         #   else: P = (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(k*(1-k))))
                
            if genPTR.get_next_float()*(1+5*a*r*3) < r*r*  (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(r*(1-r))  )):
                break
        
        
        # EFFECT ON TARGET
        cdef double W = r*E
        cdef double cos2 = (E - W)*(E + _2ELECTRON_REST_ENERGY)/E/(E - W + _2ELECTRON_REST_ENERGY)
        particles.ELECTRONS.push_back(E - W)
        particles.ELECTRONS.push_back(sqrt(cos2))
        
        # SECONDARY ELECTRON
        cos2 = (W)*(E + _2ELECTRON_REST_ENERGY)/E/(W + _2ELECTRON_REST_ENERGY)
        particles.ELECTRONS.push_back(W)
        particles.ELECTRONS.push_back(sqrt(cos2))

        
        
        #cos2 = W * (E + _2ELECTRON_REST_ENERGY)/E/(W + _2ELECTRON_REST_ENERGY)
        #cdef double cos_sec = sqrt(cos2)
        
        
        
        
        
        
                # new energy and cos for projectile
                
                

        
        # # energy and cos for secondary electron
        # particles.ELECTRONS.push_back(self.Wk - self.BE[i])
        # particles.ELECTRONS.push_back(cos_sec)
        
        
        # if self.relax:
        #     return (STATE.E - W, cos, W-self.Uk, cos_sec)         
        # return (STATE.E - W, cos, W, cos_sec)
    
        
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef double Q_minus(self, double E):
        """eq 3.83"""
        cdef double dE = E - self.Wk
        return sqrt( (self.momentum(E) - self.momentum(dE))**2  + ELECTRON_REST_ENERGY**2  ) \
               - ELECTRON_REST_ENERGY
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef double momentum(self, double E):
        return sqrt(  E * (E + _2ELECTRON_REST_ENERGY)  ) 
    
    
    
            #cdef double newE = STATE.E - self.Wk 
       # cdef double Esec = self.Wk - self.Uk if self.relax else self.Wk
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef void sampleLFAR(self, double E, mixmax_engine *genPTR, PARTICLES *particles):
        
        
        cdef double Qm = self.Q_minus(E) #
       # if Qm / _2ELECTRON_REST_ENERGY < 0:
        #    print(E, Qm, _2ELECTRON_REST_ENERGY, Qm / _2ELECTRON_REST_ENERGY, self.Wk)
        
        cdef double Qs = Qm/(1 + Qm / _2ELECTRON_REST_ENERGY)

         
        #print(_2ELECTRON_REST_ENERGY)
        

        cdef double A = Qs/self.Wk *(1 + self.Wk/_2ELECTRON_REST_ENERGY )
        cdef double B = Qs/_2ELECTRON_REST_ENERGY
        cdef double Q = Qs / (A**(genPTR.get_next_float()) - B)

        cdef double p2E = E * (E + _2ELECTRON_REST_ENERGY)
        cdef double p2d = (E - self.Wk) * ((E - self.Wk) + _2ELECTRON_REST_ENERGY)   #self.p2(E - self.Wk)
        cdef double p2Q = Q * (Q + _2ELECTRON_REST_ENERGY)

        cdef double cos = p2E + p2d - p2Q
        #cos = .5*cos/(p2E*p2d)**.5
        particles.ELECTRONS.push_back(E - self.Wk)
        particles.ELECTRONS.push_back(.5*cos/(p2E*p2d)**.5)
        
        
        
        
        
        cdef double Wk2 = self.Wk*self.Wk
        cdef double X = (E/ELECTRON_REST_ENERGY + 1)**2
        A = Wk2 / ((X-1)/X) / p2Q

        B = (p2Q - Wk2)/(2*self.Wk*(E + _2ELECTRON_REST_ENERGY))
        
        particles.ELECTRONS.push_back(self.Wk)
        particles.ELECTRONS.push_back(A * (1 + B)**2)









def rebuildgosInnerShell(this):
    cdef gosInnerShell self
    self = <gosInnerShell> gosInnerShell.__new__(gosInnerShell)
    self.Wk  = this.Wk  #cshell.Wk
    self.fk  = this.fk  #cshell.fk
    self.Nsh = this.Nsh #0

    self.INDEX = np.array(this.INDEX)
    self.BE =    np.array(this.BE)
    self.Ncol = this.Ncol

    self.ionizationCS = np.array(this.ionizationCS )
    self.totalCS = np.array(this.totalCS, order = "F")

    
    self.p1 = self.p1
    self.p2 = self.p2
    self.p3 = self.p3
    
    self.pCLOSE = np.array(this.pCLOSE, order = "F")
    self.pLFAR  = np.array(this.pLFAR , order = "F")
    self.pTFAR  = np.array(this.pTFAR, order = "F")
    return self




cimport numpy as cnp
cdef class gosInnerShell(gosShell):
    cdef public double[:, : , ::1] ionizationCS
    cdef public double[::1] BE
    cdef public int[::1] INDEX
    cdef public int Ncol
  #  cdef public int Nsh
    cdef public double kc
    cdef public double Uk
    #cdef public cnp.ndarray p1, p2, p3
    
    def __reduce__(self):
        this = MAP()
        this.Wk  = self.Wk  #cshell.Wk
        this.fk  = self.fk  #cshell.fk
        this.Nsh = self.Nsh #0

        this.INDEX = np.array(self.INDEX)
        this.BE =    np.array(self.BE)
        this.Ncol = self.Ncol

        this.ionizationCS = np.array(self.ionizationCS )
        this.totalCS = np.array(self.totalCS, order = "F")

        
        this.p1 = self.p1
        this.p2 = self.p2
        this.p3 = self.p3
        
        this.pCLOSE = np.array(self.pCLOSE, order = "F")
        this.pLFAR  = np.array(self.pLFAR , order = "F")
        this.pTFAR  = np.array(self.pTFAR, order = "F")
        return rebuildgosInnerShell, (this, )

        
    def __init__(self, CShell cshell, double Z, to_sum):
        self.Wk = cshell.Wk
        self.fk = cshell.fk
        self.Nsh = 0
        # getting ionization cross sections
        from ..database import EEDL
        
        
        totalCS = np.zeros(len(eax))
        ionizationCS = []
        self.INDEX = np.array(cshell.INDEX)
        self.BE = np.array(cshell.BE)
        self.Ncol = cshell.Ncollection
        

        for designator in cshell.designators:
            self.Nsh += 1
            table = EEDL[Z-1][(9, 81, 91, designator, 0, 0)]
            
            #ionizationIND.append(designator_to_index[designator])
            
            
            
            E = table.X    # MeV
            CS = table.Y   # barn
            
            E = E*1e6 # to eV
            CS = CS*1e-24 # to cm^2
            
            
            
            
            spline = CubicSpline(np.log(E[1:]), np.log(CS[1:]), extrapolate = False)
            
            
            
            CS = spline(np.log(eax))
            
            
            CS = np.exp(CS)
            np.nan_to_num(CS, nan=0, copy = False )
            
            totalCS += CS
            interp = makeLinLin(eax, CS)
                        
            ionizationCS.append(interp)
        #cdef float[::1] arr = <float [:cpp_vector.size()]>cpp_vector.data()
        self.ionizationCS = np.array(ionizationCS)
        to_sum += totalCS
        self.totalCS = np.array(makeLinLin(eax, totalCS), order = "F")

        #these guys already come evaluated at eax
        CLOSE = np.array(cshell.close.SIGMA0)
        LFAR  = np.array(cshell.distant.SIGMA0L)
        TFAR  = np.array(cshell.distant.SIGMA0T)
         
        norm = CLOSE + LFAR + TFAR
        norm[norm == 0] = 1
        
        self.p1 = CLOSE/norm
        self.p2 = LFAR/norm
        self.p3 = TFAR/norm
        
        self.pCLOSE = np.array(makeLinLin(eax, CLOSE/norm), order = "F")
        self.pLFAR  = np.array(makeLinLin(eax, LFAR/norm), order = "F")
        self.pTFAR  = np.array(makeLinLin(eax, TFAR/norm), order = "F")
        
        
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True)
    cdef inline void sample( self, mixmax_engine* genPTR, int index, double E, PARTICLES *particles):
        cdef double r = genPTR.get_next_float()*(self.totalCS[0, index] + E*self.totalCS[1, index])
        
        
        #which shell is going to be ionized ?
        
        
        cdef double cumul = 0

        #cdef int j
        cdef int i   
        for i in range(self.Nsh):
            cumul += self.ionizationCS[i, 0, index] + E*self.ionizationCS[i, 1, index]
            if r < cumul:
                self.Uk = self.BE[i]
                i = self.INDEX[i]
                
                break
        
        
        #print(self.rATOM)
       # print(self.ATOMptr.Z)
        self.ATOMptr.run(i, particles, genPTR)
        
        # which regime is going to happen?
        
        if E < self.Wk:
            particles.ELECTRONS.push_back(0)
            particles.ELECTRONS.push_back(0)
        
            # energy and cos for secondary electron
            particles.ELECTRONS.push_back(E - self.Uk)
            particles.ELECTRONS.push_back(0)
            return
        
        
                #double[:, ::1] pCLOSE 
       # double[:, ::1] pLFAR  
       # double[:, ::1] pTFAR 
        
        
        
        r = genPTR.get_next_float()
        
        
        cumul = self.pLFAR[0, index] + E*self.pLFAR[1, index]
        

  
            
        if r < cumul: 
            self.sampleLFAR(E, genPTR, particles)
            return
        
        #cumul += self.pLFAR[index, 0] + E*self.pLFAR[index, 1]
            
        if r < cumul + self.pCLOSE[0, index] + E*self.pCLOSE[ 1, index]:
            self.sampleCLOSE(E, genPTR, particles)
            return
        
        # new energy and cos for projectile
        particles.ELECTRONS.push_back(E - self.Wk)
        particles.ELECTRONS.push_back(1.)
        
        # energy and cos for secondary electron
        particles.ELECTRONS.push_back(self.Wk - self.Uk)
        particles.ELECTRONS.push_back(1.)

        #cdef double Esec = self.Wk - self.BE[i]
        #cdef double newE = E - self.Wk
        return #(newE, 1., Esec, -1.)
        #return
        
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True)
    cdef void sampleCLOSE(self, double E, mixmax_engine *genPTR, PARTICLES *particles):
        
        cdef double r
        cdef double kc = self.Wk/E
        cdef double a = (E/(E + ELECTRON_REST_ENERGY))**2
        cdef double CONST = 1+ 5*a*kc
        cdef double CONST2 = (1 - 2*kc)
        cdef double CONST3 = CONST2/(5*a*kc)
    
        while 1:
            r = CONST*genPTR.get_next_float()
            
            if r < 1: r = kc / (1 - r*CONST2)
            else: r = kc + (r - 1)*CONST3
            
            if r < kc: continue
            if r > .5: continue
        
         #   else: P = (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(k*(1-k))))
                
            if genPTR.get_next_float()*(1+5*a*r*3) < r*r*  (1/(r*r) + 1/(1-r)**2 - 1/(r*(1-r)) + a * (1 + 1/(r*(1-r))  )):
                break
        
        
        # EFFECT ON TARGET
        cdef double W = r*E
        cdef double cos2 = (E - W)*(E + _2ELECTRON_REST_ENERGY)/E/(E - W + _2ELECTRON_REST_ENERGY)
        particles.ELECTRONS.push_back(E - W)
        particles.ELECTRONS.push_back(sqrt(cos2))
        
        # SECONDARY ELECTRON
        cos2 = (W)*(E + _2ELECTRON_REST_ENERGY)/E/(W + _2ELECTRON_REST_ENERGY)
        particles.ELECTRONS.push_back(W - self.Uk)
        particles.ELECTRONS.push_back(sqrt(cos2))

        
        
        #cos2 = W * (E + _2ELECTRON_REST_ENERGY)/E/(W + _2ELECTRON_REST_ENERGY)
        #cdef double cos_sec = sqrt(cos2)
        
        
        
        
        
        
                # new energy and cos for projectile
                
                

        
        # # energy and cos for secondary electron
        # particles.ELECTRONS.push_back(self.Wk - self.BE[i])
        # particles.ELECTRONS.push_back(cos_sec)
        
        
        # if self.relax:
        #     return (STATE.E - W, cos, W-self.Uk, cos_sec)         
        # return (STATE.E - W, cos, W, cos_sec)
    
        
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True)  
    cdef double Q_minus(self, double E):
        """eq 3.83"""
        cdef double dE = E - self.Wk
        return sqrt( (self.momentum(E) - self.momentum(dE))**2  + ELECTRON_REST_ENERGY**2  ) \
               - ELECTRON_REST_ENERGY
               
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef double momentum(self, double E):
        return sqrt(  E * (E + _2ELECTRON_REST_ENERGY)  ) 
    
    
    
            #cdef double newE = STATE.E - self.Wk 
       # cdef double Esec = self.Wk - self.Uk if self.relax else self.Wk
    @cython.boundscheck(False)
    @cython.wraparound(False) 
    @cython.initializedcheck(False)
    @cython.cdivision(True) 
    cdef inline void sampleLFAR(self, double E, mixmax_engine *genPTR, PARTICLES *particles):
        
               # cdef double dE = E - self.Wk
        #return sqrt( (self.momentum(E) - self.momentum(dE))**2  + ELECTRON_REST_ENERGY**2  ) \
              # - ELECTRON_REST_ENERGY
              
        cdef double p2E = E * (E + _2ELECTRON_REST_ENERGY)
        cdef double p2d = (E - self.Wk) * ((E - self.Wk) + _2ELECTRON_REST_ENERGY)   #self.p2(E - self.Wk)
        
        
        cdef double Qs = sqrt( sqrt( p2E )  - p2d  + ELECTRON_REST_ENERGY*ELECTRON_REST_ENERGY  ) - ELECTRON_REST_ENERGY
        
        
        
        
       # if Qm / _2ELECTRON_REST_ENERGY < 0:
        #    print(E, Qm, _2ELECTRON_REST_ENERGY, Qm / _2ELECTRON_REST_ENERGY, self.Wk)
        
        Qs = Qs/(1 + Qs / _2ELECTRON_REST_ENERGY)

         
        #print(_2ELECTRON_REST_ENERGY)
        

       # cdef double A = Qs/self.Wk *(1 + self.Wk/_2ELECTRON_REST_ENERGY )
        #cdef double B = Qs/_2ELECTRON_REST_ENERGY
        Qs = Qs / ((   Qs/self.Wk *(1 + self.Wk/_2ELECTRON_REST_ENERGY )      )**(genPTR.get_next_float()) - Qs/_2ELECTRON_REST_ENERGY)


        cdef double p2Q = Qs * (Qs + _2ELECTRON_REST_ENERGY)

      #  cdef double cos = p2E + p2d - p2Q
        #cos = .5*cos/(p2E*p2d)**.5
        particles.ELECTRONS.push_back(E - self.Wk)
        particles.ELECTRONS.push_back(.5*(p2E + p2d - p2Q)/sqrt(p2E*p2d))
        
        
        
        
        
        #cdef double Wk2 = self.Wk*self.Wk
        #cdef double X = (E/ELECTRON_REST_ENERGY + 1)**2
        Qs = (E/ELECTRON_REST_ENERGY + 1)
        #A = Wk2 / ((X-1)/X) / p2Q

        #B = (p2Q - Wk2)/(2*self.Wk*(E + _2ELECTRON_REST_ENERGY))
        
        particles.ELECTRONS.push_back(self.Wk - self.Uk)
        particles.ELECTRONS.push_back( self.Wk*self.Wk / ((Qs*Qs-1)/(Qs*Qs) / p2Q * (1 + (p2Q - self.Wk*self.Wk)/(2*self.Wk*(E + _2ELECTRON_REST_ENERGY)))**2))
        
        #cdef double cos_sec = A * (1 + B) **2
        #return (newE, cos, Esec, cos_sec)

        
        
        
        # SECONDARY ELECTRON
        # cos2 = (W)*(E + _2ELECTRON_REST_ENERGY)/E/(W + _2ELECTRON_REST_ENERGY)
        # particles.ELECTRONS.push_back(self.Wk - self.Uk)
        # particles.ELECTRONS.push_back(sqrt(cos2))
        
        
        
        
        
        
   
    def plot(gosShell self, units = "cm2"):
        import matplotlib.pyplot as plt
        import numpy as np
        
        cdef double[:, ::1] CS
        fig = plt.figure()
        totalCS = np.array(self.totalCS)
        if units == "cm2":
            
            plt.plot(eax[:-1], totalCS[0] + totalCS[1]*eax[:-1])
            
            for CS in self.ionizationCS:
                y = np.array(CS)
                plt.plot(eax[:-1], y[0] + y[1]*eax[:-1])
                
        
            plt.ylabel("Shell Cross Section (cm^2) ")
        elif units == "barn":
            plt.plot(eax[:-1], (totalCS[0] + totalCS[1]*eax[:-1])*1e24)
            for CS in self.ionizationCS:
                y = np.array(CS)
                plt.plot(eax[:-1], (y[0] + y[1]*eax[:-1])*1e24)
                
        
            plt.ylabel("Shell Cross Section (cm^2) ")
            
            
            plt.ylabel("Shell Cross Section (barn) ")
            
            
        plt.title(f"Wk = {self.Wk}")
        plt.xlabel("Energy of Incident Electron (eV)")
        plt.xscale("log"); plt.yscale("log");
        plt.grid(which = "both")
        plt.show()


        

        
        
        
        
        
        