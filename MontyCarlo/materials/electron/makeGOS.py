# -*- coding: utf-8 -*-
"""
Created on Sun Jan 24 18:51:27 2021

@author: Rui Campos
"""

from numpy import *
from scipy.optimize import fsolve


ELECTRON_REST_MASS      = 0.51099895000e6
#taken from PENELOPE
EXCITATION = [19.2, 41.8, 40.0, 63.7,
 76.0,81.0,82.0, 95.0, 115.0, 137.0, 149.0, 156.0, 166.0, 173.0, 173.0, 180.0, 174.0, 188.0, 190.0, 191.0, 216.0,
 233.0, 245.0, 257.0, 272.0, 286.0, 297.0, 311.0, 322.0,
 330.0, 334.0, 350.0, 347.0, 348.0, 343.0, 352.0, 363.0,
 366.0, 379.0, 393.0, 417.0, 424.0, 428.0, 441.0, 449.0,
 470.0, 470.0, 469.0, 488.0, 488.0, 487.0, 485.0, 491.0,
 482.0, 488.0, 491.0, 501.0, 523.0, 535.0, 546.0, 560.0,
 574.0, 580.0, 591.0, 614.0, 628.0, 650.0,
 658.0, 674.0, 684.0, 694.0, 705.0, 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0, 800.0, 810.0, 823.0,
 823.0, 830.0, 825.0, 794.0, 827.0, 826.0, 841.0, 847.0, 878.0, 890.0, 902.0, 921.0, 934.0, 939.0, 952.0, 966.0, 980.0]



def getData(Z):
    from MontyCarlo.settings import __montecarlo__
    from MontyCarlo.tools.data import getAxis, getTable
    
    file_name = str(Z) + ".txt"
    file_path = __montecarlo__/'materials'/'EADL'/file_name
    del file_name
    path = str(file_path)
    del file_path
   # directory = str(__materials__)
    #path = directory + "\\EADL\\" + str(Z) + ".txt"
    
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

            line0 = text[0]
            Z  = int(line0[0:3])
            Aw = float(line0[13:24])


            bookmarked_text = {}

            for i, j in zip(bookmark_, _ookmarks):
                    line1, line2 = text[i], text[i+1]

                   #on line 1
                    Yi = float(line1[7:9])    #particle identifier
                    Yo = float(line1[10:12])  #secondary particle designator

                    #on line 2
                    C  = float(line2[0:2])    #reaction descriptor
                    I  = float(line2[2:5])    #reaction property
                    S  = float(line2[5:8])    #reaction modifier
                    X1 = float(line2[22:32])  #subshell designator

                    flags = (Yi, C, S, X1, Yo, I)

                    flags = tuple(map(int, flags))
                    bookmarked_text[flags] = text[i+2:j-1]
                    
    data = {Id:getTable(bookmarked_text[Id]) for Id in bookmarked_text}
    
    return Aw, data, EXCITATION[Z-1]





class Molecule:
    def __init__(self, formula, density):
        self.formula = formula
        self.density = density
        
        #PREPPING FOR A LONG INIT
        self.ATOMS = []
                
        #CONSTRUCTING ATOMS AND RELEVENT QUANTITIES
        for Zx, x in formula.items():
           # atom = Atom(Zx, x, formula.relax[Zx])
            atom = Atom(Zx, x)

            self.ATOMS.append(atom)  
        
    def __repr__(self):
        repres = ""
        for atom in self:
            repres += f"    {atom.x}x{atom} \n"
        return repres
    def __iter__(self):
        yield from self.ATOMS
        
    def __getitem__(self, i):
        return self.ATOMS[i]

class Atom:
    SHELLS = []
    
    def __init__(self, Z, x):
        self.x = x
        self.Z = Z
        self.Aw, data, self.I = getData(Z)
        
        # GETTING SHELL DATA 
        DATA             = data[(0, 91, 0, 0, 0, 912)]
        DESIGNATORS      = [int(design) for design in DATA[:, 0]]
        EL_NUMBERS       = DATA[:, 1]
        BINDING_ENERGIES = data[(0, 91, 0, 0, 0, 913)][:, 1]
        
        # CREATING SHELLS
        SHELLS = []
        for i in range( len(DESIGNATORS) ):

            
            SHELLS.append(
                          Shell(index = i,
                                i  = DESIGNATORS[i]      , 
                                fk = EL_NUMBERS[i]       ,  
                                Uk = BINDING_ENERGIES[i]*1e6 ,
                                )
                          )    
            
        self.SHELLS = SHELLS   
    
        
    
    def __iter__(self):
        yield from self.SHELLS

    def __repr__(self):
        
        repres = f"<Atom Z={self.Z}, Aw = {self.Aw} amu, I = {self.I} eV> \n"
        for shell in self:
            repres += f"    {shell} \n"

        
        
        return repres
    
    def __len__(self):
        return len(self.SHELLS)
    
    def __getitem__(self, aslice):
        return self.SHELLS[aslice]
    
    def pop(self, i):
        return self.SHELLS.pop(i)




class Shell:
    Wk = None
    def __init__(self,  index = 0, i = 0, fk = 0, Uk = 0):
        self.INDEX = [index]
        self.BE = [Uk]
        self.designators = [i]
        self.i = i
        self.fk = fk
        self.Uk = Uk
    def __repr__(self):
        return f"<Shell #{self.i}, fk = {self.fk}, Uk = {self.Uk} eV, Wk = {self.Wk} eV>"
    

    
    def __bool__(self):
        return not self.empty





def select_shells(molecule):
    """
    Remove outer shells of all atoms in molecule. Skip hydrogen and helium.
    """
    
    shells = []
    for atom in molecule:
        if len(atom) < 2:
            continue
            
        #select shell that will contribute to plasmon oscillations
        shells.append(atom.pop(-1))
    return shells 


Na = 6.0221409e+23

H_BAR                   = 6.5821e-16 #eV s
MASS_ELECTRON           = 9.1094e-28 #g
ELECTRON_CHARGE         = 4.8032e-10 #statC
from numpy import pi


def makeCB(molecule, selected_shells):
    
    
    
    Am = sum(atom.Aw*atom.x for atom in molecule)    #total atomic mass
    N = molecule.density * Na / Am                   #number density
    
    Ztot = sum(atom.x*atom.Z for atom in molecule)          #total number of electrons
        
    omega2 = 4*pi*N*Ztot*H_BAR**2 * ELECTRON_CHARGE**2 / MASS_ELECTRON
    omega = sqrt(omega2)
    
    

    #CREATOMG CONDUCTION BAND
    #cb = select_shells(molecule)
    fcb = sum(shell.fk for shell in selected_shells)
    Wcb = (fcb/Ztot)**.5 * omega
    
    cb = Shell(i = 0, Uk = 0, fk = fcb)
    cb.Wk = Wcb
    
    
    #INSERTING USEFUL QUANTITIES
    molecule.omega = omega
    molecule.N = N
    molecule.Ztot = Ztot
    return cb



def calculate_ressonance(molecule, cb):
    from scipy.optimize import fsolve
    from numpy import exp, log
    
    
    ZlnI = sum(atom.x*atom.Z*log(atom.I) for atom in molecule)
    
    I = exp(ZlnI/molecule.Ztot)
    
    molecule.ZlnI = ZlnI
    molecule.I = I
    
    
    A = cb.fk*log(cb.Wk) - ZlnI if cb else -ZlnI
    B  = 2*molecule.omega**2 / 3 / molecule.Ztot
    
    def eqn(a):
        _sum = 0
        for atom in molecule:
            
            from_atom = 0
            for shell in atom:
                val = (a*shell.Uk)**2 + B*shell.fk
                from_atom += shell.fk * log(val)
            
            _sum += atom.x*from_atom
            
            
            
                #val = (a*shell.Uk)**2 + B*shell.fk*atom.x
                #_sum += atom.x*shell.fk*log( (a*shell.Uk)**2 + B*shell.fk*atom.x )
                #_sum += atom.x*( shell.fk*log( (a*shell.Uk)**2 + B*shell.fk*atom.x )
    
        return A +.5*_sum
    
    res = fsolve(eqn, exp(1/2))
    
    a = res[0]
    molecule.a = a
    print(f"        > a = {a}")
    
    
    
    
    omega = molecule.omega
    Ztot = molecule.Ztot
    
    
    for atom in molecule:
        for shell in atom:
            shell.Wk = ((shell.Uk*a)**2 + 2/3 * shell.fk * omega**2  / Ztot)**.5
            
    

    _ZlnI = cb.fk*log(cb.Wk) if cb else 0
    
    for atom in molecule:
        for shell in atom:
            _ZlnI += atom.x*shell.fk*log(shell.Wk)
            

    err = (ZlnI - _ZlnI)/ZlnI *100
    print("                err = ", err, "%")
    
    return res[0]
    
    


from numpy import array

def makeF(molecule, cb):

    
    SHELLS = []
    for atom in molecule:
        for shell in atom:
            SHELLS.append(shell)
    SHELLS.append(cb)
    

    fk = array([shell.fk for shell in SHELLS])
    Wk = array([shell.Wk for shell in SHELLS])
    
    omega = molecule.omega
    Ztot  = molecule.Ztot
    
    def F(L):
        return omega**2 / Ztot * sum(fk/(Wk**2 + L**2))
    
    return F, fk, Wk


def makeDelta(molecule, fk, Wk):
    omega = molecule.omega
    Ztot = molecule.Ztot
    
    def delta(L, beta2):
        return sum(fk*log(1 + L**2/Wk**2))/Ztot - L**2/omega**2 * (1 - beta2)
    
    return delta

def getLfrom(F, beta2):
    def eqn(L):
        return F(L) - 1 + beta2
    
    return fsolve(eqn, 5)[0]

def calculate_ZlnI(molecule, cb):
    ZlnI = cb.fk*log(cb.Wk) if cb else 0
    for atom in molecule:
        for shell in atom:
            ZlnI += shell.fk * log(shell.Wk)
    return ZlnI


def newShellFrom(group):
    if not group:
        return []
    
    if len(group) == 1:
        return group
    
    
 
    
    ZlnI = 0
    Ztot = 0
    Uk = 0
    for shell in group:
        #print(shell)
        Uk += shell.Uk
        Ztot += shell.fk
        ZlnI += shell.fk*log(shell.Wk)
    
    Uk = Uk/len(group) #Uk is not needed, just keeping it for book keeping
    Wk = exp(ZlnI/Ztot)
    
    shell = Shell(i = group[0].i, fk = Ztot, Uk = Uk)
    shell.designators = []
    shell.INDEX = []
    shell.BE = []
    for _shell in group:
        shell.designators += _shell.designators
        shell.INDEX += _shell.INDEX
        shell.BE += _shell.BE
        
    
    
    
    shell.Wk = Wk
    return [shell]

def newMolecule(molecule):
    
    for atom in molecule:
        if atom.Z == 1: continue
        if atom.Z <= 27: #only K shell has binidng energy ABOVE 1keV, group everything else!
            newSHELLS = [atom[0]] #K shell
            newSHELLS += newShellFrom(atom[1:])
            atom.SHELLS = newSHELLS
            continue
        if atom.Z <= 51: #condition to classify L shells as inner
            newSHELLS = [atom[0]]
            
            innerL, outerL = [], []
            for L in atom[1:4]:
                if L.Uk < 1e3: outerL.append(L)
                else: innerL.append(L)
            
            if len(innerL) == 1:
                newSHELLS += innerL[0]
                newSHELLS += newShellFrom(outerL + atom[4:])
            else:
                newSHELLS += newShellFrom(innerL)
                newSHELLS += newShellFrom(outerL + atom[4:])
            
            atom.SHELLS = newSHELLS
            continue

        if atom.Z <= 84: #condition to classify M shells as inner
            newSHELLS = [atom[0]] # keeping the K shell
            newSHELLS += newShellFrom(atom[1:4]) # if we are classifying M shells as inner, L shells are automatically inner

            innerM, outerM = [], []
            for M in atom[4:9]:
                if M.Uk < 1e3: outerM.append(M)
                else: innerM.append(M)

            if len(innerM) == 1:
                newSHELLS += innerM[0]
                newSHELLS += newShellFrom(outerM + atom[9:])
            else:
                newSHELLS += newShellFrom(innerM)
                newSHELLS += newShellFrom(outerM + atom[9:])
            
            atom.SHELLS = newSHELLS
            continue

        #everyone else has N shells with binding energy above 1keV




        newSHELLS = [atom[0]] # K shell 
        newSHELLS += newShellFrom(atom[1:4]) # L shell
        newSHELLS += newShellFrom(atom[4:9]) # M shell
        inner, outer = [], []
        for shell in atom[9:]:
            if shell.Uk < 1e3: outer.append(shell)
            else: inner.append(shell)
         
        newSHELLS += newShellFrom(inner)
        newSHELLS += newShellFrom(outer)
        atom.SHELLS = newSHELLS
        continue


        #newSHELLS = [atom[0]]
        #newSHELLS += newShellFrom(atom[1:4])
        #newSHELLS += newShellFrom(atom[4:9])
        #newSHELLS += newShellFrom(atom[9:])
        
        #atom.SHELLS = newSHELLS
        
def pyGOS(formula, density):
    from builtins import sum
    molecule = Molecule(formula, density)
    shells = select_shells(molecule)

    


        

    cb = makeCB(molecule, shells)
    
    if shells == []:
        cb.empty = True
    else:
        cb.empty = False
        
    calculate_ressonance(molecule, cb)
    
    F, fk, Wk = makeF(molecule, cb)
    delta = makeDelta(molecule, fk, Wk)
    delta1 = lambda beta2: delta(getLfrom(F, beta2), beta2)
    #we have F and delta! getLfrom

    old = calculate_ZlnI(molecule, cb)
    newMolecule(molecule)
    new = calculate_ZlnI(molecule, cb)


    formula.log.add_paragraph(f"Ressonance energy of cb shell: {cb.Wk}")
    for atom in molecule:
        formula.log.add_paragraph(f"-----> Z = {atom.Z}")
        for shell in atom:
            formula.log.add_paragraph(f"------------> Wk = {shell.Wk}")

    #print(old, new)
    #print( (old - new)/old * 100, "%" )
    
    return molecule, cb, delta1
        
        
        
