# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 12:35:38 2020

@author: Rui Campos
"""



from numpy import * #array, geomspace, flip, load, searchsorted
from numpy.random import *

from ..particles.photons import choose
from ..settings import __montecarlo__
from ..tools.data import getAxis, getTable
from ..tools.interpol1 import LinLinInterpolation

from . import database as db

__materials__ = __montecarlo__/'materials'
directory = str(__materials__)


class Vacancy:
    def __init__(self, position):
        self.current_shell = position
        #self.current_shell.Nel -= 1
        self.current_shell.Nvac += 1
        
        
    @property
    def __hash__(self):
        return self.current_shell.i
    
    
    def move(self, shell):
        #self.current_shell.Nel += 1
        self.current_shell.Nvac -= 1
        #shell.Nel -= 1
        shell.Nvac += 1
        self.current_shell = shell
    
    def sampleTransition(self):
        return self.current_shell.chooseTransition()
    
    def __lt__(self, other):
        return self.__hash__ < other.__hash__
    
    def __gt__(self, other):
        return self.__hash__ > other.__hash__
    
    def __eq__(self, other):
        return self.__hash__ == other.__hash__
    
    def __int__(self):
        return self.__hash__





class Shell:
    def __init__(self, *args, **kwargs):
        
        self.Nvac = 0
        self.args = args
        self.i, self.Nel, self.binding_energy, self.KE, self.avgR, self.Z = args
        
        
        self.CS = db.EPDL[self.Z-1][(7, 73, 91, self.i, 0, 0)].getLinLinInterpol()
        
        
        if not self.Nel.is_integer():
            self.Nmax = 1
        else:
            self.Nmax = self.Nel
        
        self.kwargs = kwargs
        
        self.rTrans = kwargs['rTrans']
        self.nrTrans = kwargs['nrTrans']
        
        if self.rTrans is not None and self.nrTrans is not None:
            self.rTrans = insert(self.rTrans, 1, 0, axis = 1)

            #self.trans = append(self.rTrans, self.nrTrans)
            self.trans = concatenate((self.rTrans, self.nrTrans), axis = 0)
            
            self.Emax = max(self.trans[:, 3])
        
            self.setCumul()
            
            def chooseTransition():
                i = searchsorted(self.cumul, rand(), side = 'right') - 1
                j, k, _, E = self.trans[i, :]
                return j, k, E
            
            
            self.chooseTransition = chooseTransition
            
        else: 
            self.Emax = 0
            self.chooseTransition = lambda : (0, 0, 0)
            self.trans = None
        
        # if self.rTrans is not None and self.nrTrans is not None:
        #     self.rTj = self.rTrans[:, 0]
        #     self.rTProb   = self.rTrans[:, 1]
        #     self.rTEnergy = self.rTrans[:, 2]
        # else:
        #     self.rTProb =[0]
        #     self.rTj = []
        
        # if self.nrTrans is not None:
        #     self.nrTj      = self.nrTrans[:, 0]
        #     self.nrTk      = self.nrTrans[:, 1]
        #     self.nrTProb   = self.nrTrans[:, 2]
        #     self.nrTEnergy = self.nrTrans[:, 3]
        # else:
        #     self.nrTProb = [0]
        #     self.nrTj = []

 
        
    def setTransitions(self, atom):
        
        if self.trans is None:
            return None
        
        TRANS = []
        for trans in self.trans:
            j, k, p, E = trans
            TRANS += [[atom[j], atom[k], p, E]]
        self.trans = array(TRANS)
        
    def setCumul(self):
        probs = []

            
        probs = self.trans[:, 2]
        
        self.cumul = [sum(probs[0:i]) for i in range(len(probs))]
        
        
    def getNel(self):
        Nel = self.Nel - self.Nvac
        if Nel < 0: return 0
        return Nel
        
        
            
            
    def chooseTransition(self):
        
        i = searchsorted(self.cumul, rand())
        j, k, _, E = self.trans[i, :]
        
        return j, k, E
        

class Atom:
    def __init__(self, Z):
        self.Z = Z
        self.vacancies = []
        self.path = directory + "\\EADL\\" + str(Z) + ".txt"
        self.Aw, self.EADL_dict = self.getBookmarkedText()
        
        self.data = {Id:getTable(self.EADL_dict[Id]) for Id in self.EADL_dict}
        
        number_el = self.data[(0, 91, 0, 0, 0, 912)]
        shells = [int(design) for design in number_el[:, 0]]
        
        number_el = number_el[:, 1]
        
        
        binding_energy = self.data[(0, 91, 0, 0, 0, 913)][:, 1]
        KE    = self.data[(0, 91, 0, 0, 0, 914)][:, 1]
        avgR  = self.data[(0, 91, 0, 0, 0, 915)][:, 1]
        #rLW  = self.data[(0, 91, 0, 0, 0, 916)]
        #nrLW = self.data[(0, 91, 0, 0, 0, 917)]
    
        self.SHELLS = []
        for i in range(len(shells)):
            
            self.SHELLS += [Shell(shells[i], number_el[i],  binding_energy[i], KE[i], avgR[i], Z,
                                  rTrans  = self.data.get((0, 92, 91, shells[i], 7, 931), None),
                                  nrTrans = self.data.get((0, 92, 91, shells[i], 9, 932), None)
                                  )]
        
        for shell in self.SHELLS:
            shell.setTransitions(self)
            
        self.CS = db.EPDL[Z-1][(7, 73, 0, 0, 0, 0)].getLinLinInterpol()
    def __iter__(self):
        yield from self.SHELLS
        
    def printShells(self):
        for shell in self.SHELLS:
            print(f"------{shell.i}--------|| ", shell.Nvac*"* | " + int(shell.Nel - shell.Nvac)*"e | ")
            #print("#VACANCIES = ", shell.Nvac)
            #print("#ELECTRONS = ", shell.Nel)




    def _introduceVacancy(self, shell):
        self.vacancies += [Vacancy(shell)]
                
    def introduceVacancy(self, i):
        for shell in self.SHELLS:
            if shell.i == i:
                self.vacancies += [Vacancy(shell)]
    
    def __getitem__(self, i):
        for shell in self.SHELLS:
            if shell.i == i:
                return shell
        return None
        
    
        
    def run(self):
        spectrum = {'el':[], 'ph':[]}
        
        while len(self.vacancies) > 0:
            self.vacancies.sort()
            #print("___________________________________________________")
            #self.printShells()

            activeVac = self.vacancies[0]
            if activeVac.current_shell.i > 15:
                break
            
            if activeVac.current_shell.trans is None:
                self.vacancies.remove(activeVac)
                continue
                
            
            possible_transitions = []
            
            for transition in activeVac.current_shell.trans:
                shell_i, shell_j = transition[0], transition[1]
                
                
                if shell_j is None:
                    if shell_i.Nel - shell_i.Nvac > 0:
                        possible_transitions += [transition]
                        continue
                elif shell_i == shell_j:
                    if shell_i.Nel - shell_i.Nvac > 1:
                        possible_transitions += [transition]
                else:
                    if  shell_i.Nel - shell_i.Nvac > 0 \
                    and shell_j.Nel - shell_j.Nvac > 0:
                        possible_transitions += [transition]

                    
                
                    
            if len(possible_transitions) == 0:
                self.vacancies.remove(activeVac)
                continue
            
            if len(possible_transitions) == 1:
                j, k, _, E = possible_transitions[0]
                
            else:
                possible_transitions = array(possible_transitions)
                probs = possible_transitions[:, 2]
                probs = probs/sum(probs)
                
                cumul = [0.]
                C = 0
                
                for p in probs:
                    C += p
                    cumul += [C]
        
                #cumul = [sum(probs[0:i]) for i in range(len(probs))]
                for _ in range(100_000):
                    
                    i = searchsorted(cumul, rand(), side ='left')
                    j, k, _, E = possible_transitions[i-1, :]
                    #j, k, E = activeVac.sampleTransition()
                    
                    p2 = (k.Nel - k.Nvac) / k.Nel if k is not None else 1
                    p1 = (j.Nel - j.Nvac) / j.Nel
    
                    if rand() < p1*p2:
                        break
                else:
                    print(possible_transitions)
                    print(activeVac.current_shell.i)
                    print("p1*p2 = ", p1p2)
                    raise RuntimeError("exceeded number of iter")
                
            if j is None: #or j > 45:
                self.vacancies.remove(activeVac)
            elif k is None:
                spectrum['ph'] += [E]
                activeVac.move(j)
  
            else:
                activeVac.move(j)
                self._introduceVacancy(k)
                spectrum['el'] += [E]

        #self.printShells()                
        self.reset()
        return spectrum
        
        
        
    def reset(self):
        for shell in self.SHELLS:
            shell.Nvac = 0
    
    
    def getBookmarkedText(self):
        path = self.path

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
                        
        return Aw, bookmarked_text