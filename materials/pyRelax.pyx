#distutils: language = c++
# cython: c_string_type=unicode, c_string_encoding=utf8
# cython: profile=False

print("pyRelax")

from .cppRelaxAPI cimport Shell as rShell
from .cppRelaxAPI cimport Atom as rAtom
from .cppRelaxAPI cimport Transition as rTransition
from .cppRelaxAPI cimport setADRESS_RAD
from .cppRelaxAPI cimport setADRESS_NONRAD
from .._random.random cimport genPTR
from .cppRelaxAPI cimport PARTICLES
from libc.string cimport memcpy 


from libcpp.deque cimport deque
from libcpp.vector cimport vector


from cython.operator cimport dereference as deref

import numpy as np #array, geomspace, flip, load, searchsorted
from numpy.random import *

from ..settings import __montecarlo__, __photonCUTOFF__, __electronCUTOFF__

cdef double photonCUTOFF = __photonCUTOFF__
cdef double electronCUTOFF = __electronCUTOFF__*1e-6


from ..tools.data import getAxis, getTable
from ..tools.CubicInverseTransform import makeAlias

#from ..tools.interpol1 import LinLinInterpolation
#from ..tools cimport search

#from ..particles.electrons cimport Electron
#from ..particles.photons cimport Photon


from . import database as db
from libc.stdlib cimport malloc, free
__materials__ = __montecarlo__/'materials'


cdef str directory = str(__materials__)

from numpy import searchsorted, concatenate, append



cimport cython
from libcpp.string cimport string



#from ..tools.interpol1 cimport LinLinInterpolation

#cimport numpy as cnp

import numpy as np   

cdef class Atom:
    def __init__(self, Z, double cut_off):
        cdef int i

        #data processing in python
        self.Z = Z
        self.path = directory + "\\EADL\\" + str(Z) + ".txt"
        self.Aw, self.EADL_dict = self.getBookmarkedText()
        
        self.data = {Id:getTable(self.EADL_dict[Id]) for Id in self.EADL_dict}
        
        number_el = self.data[(0, 91, 0, 0, 0, 912)]
        DESIGNATORS = [int(design) for design in number_el[:, 0]]
        DESIGNATORS = np.array(DESIGNATORS)
        
        number_el = number_el[:, 1]
        number_el = np.ceil(number_el)
        binding_energy = self.data[(0, 91, 0, 0, 0, 913)][:, 1]*1e6
        
        binding_energy = np.array(binding_energy)
        
        # FOR PYTHON
        self.BE = binding_energy
        self.number_el = number_el
        ################################

        #permutation = binding_energy.argsort()
        
        #binding_energy = np.flip(binding_energy[permutation])
        #number_el = np.flip(number_el[permutation])
        #DESIGNATORS = np.flip(DESIGNATORS[permutation])

        
       # for x, y in zip(number_el, binding_energy):
      #      print(x, y)
        
        ## MERGING SHELLS WITH FRACTIONAL ELECTRONS
        #  find ocurrences
        # fractional_part = number_el - np.floor(number_el)
        # condition = fractional_part > 0.
        
        # translate = dict()
        # translate_indexes = dict()
        # print(number_el)
        # print(len(number_el))
        
        # #if there are fractional electrons, find them
        # if sum(condition) != 0:
        #     #construct a set, it will automatically remove repeated values
        #     binding_energies = {Uk for Uk in binding_energy[condition]}
        #     binding_energies = list(binding_energies)
        #     binding_energies = np.array(binding_energies)
            
        #     indexes = np.arange(0, len(binding_energy), 1)
        #     indexes = indexes[condition]
            
        #     collect = dict()
        #     for m in indexes:
        #         collect[binding_energy[m]] = m
        #     print(collect)
            
        #     for Uk in binding_energies:
        #         print(Uk)
        #         #shell_locations = binding_energy[binding_energy == Uk]
                
        #         #indexes = np.arange(0, len(binding_energy), 1)
        #         #indexes = indexes[binding_energy == Uk]
                
        #         i0 = indexes[0] #first shell is kept
        #         designator0 = DESIGNATORS[i0]
                
        #         translate[designator0] = designator0
        #         translate_indexes[i0] = i0
                
        #         for i in indexes[1:]:
        #             translate[DESIGNATORS[i]] = designator0
        #             translate_indexes[i] = i
                
        #         #just checking if everything is alright
        #         print("These numbers should be consecutive:", indexes)
        #         for m in indexes:
        #             print(Uk, binding_energy[m], m)
        #             if not binding_energy[m] == Uk: print("assert??")
        #             assert binding_energy[m] == Uk

        #         new_Nel = sum(number_el[i0: i0+len(indexes)])
        #         number_el[i0] = new_Nel
                
        #         #perform surgery
        #         to_remove = indexes[1:]
        #         print(len(number_el))
        #         number_el = np.delete(number_el, to_remove)
        #         binding_energy = np.delete(binding_energy, to_remove)
        #         DESIGNATORS = np.delete(DESIGNATORS, to_remove)
        #         print(len(number_el))
        # #hopefully this has preserverd the order of the data
        # self.translate = translate
        #guarantee c contiguity
        DESIGNATORS    = np.ascontiguousarray(DESIGNATORS)
        number_el      = np.ascontiguousarray(number_el)
        binding_energy = np.ascontiguousarray(binding_energy)
                
        
        self.Nsh = len(DESIGNATORS)
        
        cdef int Ntransitions = 0
        #I have no choice  :((((((( 
        for i in range(self.Nsh):
            designator = DESIGNATORS[i]
            rTrans  = self.data.get((0, 92, 91, designator, 7, 931), None)
            nrTrans = self.data.get((0, 92, 91, designator, 9, 932), None)
            if rTrans is None and nrTrans is None: continue
            Ntransitions += len(rTrans) + len(nrTrans)
            
        self.rTRANSITIONarr = <rTransition*>malloc(Ntransitions * sizeof(rTransition))
       # print(Ntransitions)
        
        
        #print(number_el)
        
        
        #initialize all shells - store their references
        self.rSHELLarr = <rShell*>malloc(self.Nsh * sizeof(rShell))
        
        
        if not self.rSHELLarr: raise MemoryError()
        
        
        self.DESIGNATORS = np.array(DESIGNATORS, dtype = float)
        
        
        
        
        #cdef double *frac = <double*>malloc(Ntr * sizeof(double))
        
        
        cdef double[::1] frac
        cdef double test_frac
        cdef int kl
        cdef int ll
        cdef int size
        
        
        size = <int>(sum(number_el + 1))
        
        
        self.temp_frac = <double*>malloc(size * sizeof(double))
      #  print(<int>self.temp_frac)
        
        if self.temp_frac == NULL: raise MemoryError()
       # else: print(">>>  rTRANSITIONarr  :: sucessful allocation ")
        cdef rShell _SHELL
        size = 0
        cdef bint dontSIMULATE
        for i in range(self.Nsh):
            
           # print("size", size)
            # what about those 1.154 electrons in a shell? :s
            #frac = np.arange(0, size + 1 , dtype = np.float ) / size
            #print(frac)
            
            for ll in range(size, size + int(number_el[i]) + 1):
            #for ll in range(i,  i + size + 1):
                
                #print(frac[ll - i])
                #self.temp_frac[ll] = frac[ll - i]
                self.temp_frac[ll] = <double>(ll-size)/int(number_el[i])
                
               # print("in loop frac", (ll-size)/int(number_el[i]), ll)
               #print("stored:", self.temp_frac[ll])
            
            # #print(np.array(frac))
            kl = size + int(number_el[i])
          #  print("final index? kl = ", kl)
            
            test_frac = deref(&self.temp_frac[0] + kl )
           # print("test_frac", test_frac)
           # print("BE", binding_energy[i], cut_off, binding_energy[i] <= cut_off)
            
       
            self.rSHELLarr[i] = rShell(&self.temp_frac[0] + kl ,  binding_energy[i] <= cut_off)
            
            
            #self.rSHELLarr[i] = rShell(&self.temp_frac[0] + kl , binding_energy[i], self.Nsh - i - 1)
            
            # test_frac = deref(self.rSHELLarr[i].frac)
            # print("test_frac", test_frac)
            # print("")
            
            # if i != 0: 
            #     test_frac = deref(self.rSHELLarr[i-1].frac)
            #     print("prev test_frac", test_frac)
            #     print("")
            # print("ADRESS", <int>self.rSHELLarr[i].frac)
            size += int(number_el[i]) + 1
            
            
            

        
        
        
        self.rSHELLarr[self.Nsh-1].LAST = 1
        
        for i in range(self.Nsh):
            test_frac = deref(self.rSHELLarr[i].frac)
           # print(test_frac)
            #print("test_frac", test_frac)
           # print("ADRESS", <int>self.rSHELLarr[i].frac)
           # print("")
        #self.temp_frac = <double*>malloc(size * sizeof(double))
        
        &self.rSHELLarr[0]
        
        self.rATOM = rAtom(&self.rSHELLarr[0], self.Nsh, Ntransitions)
        
        

 

        # use list of references to set transitions
        cdef rShell* rSHELLptr
        cdef rShell* SHELLptr1

    
        cdef double p, E
       

        
        self.NTR = []
        self.PROBS = []
        cdef int Ntr
        cdef int ii
        
        designator = None
        
        
        #self.to_free =  <rTransition**>malloc(self.Nsh * sizeof(rTransition*))
        
        

         
        
        size = 0
        cdef int ind = 0
        
      #  print("ashdjhsakjhdsakjh -------------")
        for i in range(self.Nsh):
            previous_designator = designator
            
            designator = DESIGNATORS[i]
            
            if designator == previous_designator: 
                print("this happened, line 334")
                continue
            
            
            rSHELLptr = self.fetchFD(DESIGNATORS[i])

            
             
            
            
           # print("FRAC:", deref(rSHELLptr.frac))
            rTrans  = self.data.get((0, 92, 91, designator, 7, 931), None)
            nrTrans = self.data.get((0, 92, 91, designator, 9, 932), None)
            
            if rTrans is None and nrTrans is None:
                self.NTR.append(0)
                rSHELLptr.dontSIMULATE = 1
                #rSHELLptr.Nt = 0
                
                previous_designator = designator
                
                continue
            
            
            
            
            
            
            
            
            rTrans = np.insert(rTrans, 1, 0, axis = 1)

            _trans = append(rTrans, nrTrans)
            _trans = concatenate((rTrans, nrTrans), axis = 0)

            
            probs = _trans[:, 2]
            probs = probs/sum(probs)
            self.PROBS.append(probs)
 
            indexes = np.arange(0, len(probs))
            ALIASES = makeAlias(indexes, probs)
            
            ## construct array of transitions for this shell
            
            Ntr = len(_trans)
            self.NTR.append(Ntr)
           
           
            #rTRANSITIONarr = <rTransition*>malloc(Ntr * sizeof(rTransition))
            #self.to_free.push_back(rTRANSITIONarr)
            

            RAD = []
            for ii in range(Ntr):
                
                #ptr = rTRANSITIONarr[i]
                
                trans = _trans[ii]
                j, k, p, E = trans[0], int(trans[1]), trans[2], trans[3]
                E = 1e6*E
                if k == 0:
                    #rSHELLptr.constructRadiative( ii, E, self.fetchFD(j))
                    
                    
                    #self.rTRANSITIONarr[size + ii] = constructRadiative(E, self.fetchFD(j))
                    RAD.append(0)
                    
                    
                    
                    #rSHELLptr.newTransition(SHELLptr1, &NULL_SHELL, p, E, True)
                    

                    self.rTRANSITIONarr[size + ii].j = self.fetchFD(j)
                    #rTRANSITIONarr[ii].k = &NULL_SHELL
                    self.rTRANSITIONarr[size + ii].E = E
                    
                    setADRESS_RAD(  &self.rTRANSITIONarr[size + ii]  )
                    
                    
                    #self.rTRANSITIONarr[size + ii].perform =  getRADadress()
                    
                    continue
                #rSHELLptr.constructnonRadiative(ii, E, self.fetchFD(j), self.fetchFD(k))
                #self.rTRANSITIONarr[size + ii] = constructnonRadiative(E, self.fetchFD(j), self.fetchFD(k))
                RAD.append(1)
                
                
                #rSHELLptr.newTransition(self.fetchFD(j), self.fetchFD(k), p, E, False)
                self.rTRANSITIONarr[size + ii].j = self.fetchFD(j)
                self.rTRANSITIONarr[size + ii].k = self.fetchFD(k)
                self.rTRANSITIONarr[size + ii].E = E
                setADRESS_NONRAD(  &self.rTRANSITIONarr[size + ii]  )
            
            #print(ALIASES)
            for ii in range(Ntr):
                #rSHELLptr.set_alias(ii, ALIASES[ii][1],  ALIASES[ii][2],  RAD[ii])
                self.rTRANSITIONarr[size + ii].p = <double>ALIASES[ii][1]
                self.rTRANSITIONarr[size + ii].t = &(self.rTRANSITIONarr[size + <int>ALIASES[ii][2]])   
                
            
            rSHELLptr.setTRANSITIONS(&self.rTRANSITIONarr[size], Ntr - 1)
            size += Ntr
            
            
        for i in range(self.Nsh):
            test_frac = deref(self.rSHELLarr[i].frac)
            #print(test_frac)
        
            
            
            
    def __dealloc__(self):
      
        free(self.rTRANSITIONarr)
        free(self.temp_frac)
        free(self.rSHELLarr)

        
                

                    
    def sampleT(self, int shell_index, int N):
        
        from numpy.random import rand
        cdef rShell* shell = self.fetchFI(shell_index)
        
        import time
        cdef double t0, tf
        cdef int i
        cdef double r
        cdef int[::1] sample = np.arange(0, N)
        
        t0 = time.perf_counter_ns()
        for i in range(N):
            shell.sample_transition()
        tf = time.perf_counter_ns()
        print(f"{(tf-t0)/N}ns per run") 
        
    def sample(self, int shell_index, int N):
        
        from numpy.random import rand
        cdef rShell* shell = self.fetchFI(shell_index)
        
        import time
        cdef double t0, tf
        cdef int i
        cdef double r
        cdef int[::1] sample = np.arange(0, N)
        
        t0 = time.perf_counter_ns()
        for i in range(N):
            sample[i] = <int> shell.sample_transition()
        tf = time.perf_counter_ns()
        print(f"{(tf-t0)/N}ns per run") 
        return sample

            
        


        
    cdef rShell* fetchFD(self, int designator):
         """Takes a EADL shell designator and returns a reference to an rShell object"""
         cdef int i
         for i in range(self.Nsh):
             if self.DESIGNATORS[i] == designator:
                 return self.rATOM.fetchFI(i)
         else: 
             raise ValueError(">>>  pyRelax.pyx > .fetchFD > SHELL NOT FOUND designator = " + str(designator))
         
         
     

    cdef rShell* fetchFI(self, int index):
         return self.rATOM.fetchFI(index)
     
        
     
        
    def __len__(self): return int(self.Nsh)
     
     

        
        
        
    def __str__(self):
        rep = f"<Atom Z = {self.Z}, Aw = {self.Aw}, Nsh = {self.Nsh}> \n"
        
  
        
        for i in range(self.Nsh):
            rep += f"    <Shell #{i} #el = {self.number_el[i]}, E = {self.BE[i]}> \n"
        
        return rep
        
        
        to_print = ""
        to_print += f"<Atom Z={self.Z}> \n"
        cdef rShell shell
        cdef rShell *ptr
        cdef int i
        cdef double f
        for i in range(self.Nsh):
            print(i, self.Nsh)
            ptr = self.fetchFI(i)
            
            shell = deref(ptr)
            print(shell.LAST)
            
            #if shell.LAST: return to_print
            if shell.dontSIMULATE: continue
            f = deref(shell.frac)
            print(f)
            

            print(shell.dontSIMULATE)
            
            print(shell.Nt)
            
            
            
            print(shell.trSTART.E)

            
            #to_print += f"    <Shell #{i} @{<int>ptr} #frac = {deref(shell.frac)}  >\n"

            # to_print += <str>(shell.list_transitions())
             
        return to_print
    
    
    
    def EMPTYrunT(self, int shell_index, double CUT_OFF, int N):
        import time
        cdef double t0, tf
        cdef int i
                
        t0 = time.perf_counter_ns()
        for i in range(N):
            pass
        tf = time.perf_counter_ns()
        print(f"{(tf-t0)/N}ns per run")    
    
    cdef void run(self, int index, PARTICLES* particles, mixmax_engine *genPTR):

        #print("RUNNING RELAX")
        self.rATOM.run(index, particles, genPTR)
        
        
        
    def run2(self, int shell_index, int N):
        import time
        cdef double t0, tf
        cdef int i
        cdef PARTICLES particles

        
        t0 = time.perf_counter_ns()
        for i in range(N):
            self.rATOM.run(shell_index, &particles, genPTR)

            
        tf = time.perf_counter_ns()
        print(f"{(tf-t0)/N}ns per run")
        return (tf-t0)/N
            
    def get_spectrum(self, int shell_index, int N = 100_000):
        
        import time
        cdef double t0, tf
        cdef int i
        cdef PARTICLES particles

        cdef double E
        cdef int size 
        
        photon = {}
        electron = {}
        cdef int _
        t0 = time.perf_counter_ns()
        for i in range(N):

            
            self.rATOM.run(shell_index, &particles, genPTR)
            size = particles.PHOTONS.size()
        
            
            for _ in range(size):
                
                E = particles.PHOTONS.back()
                particles.PHOTONS.pop_back()


                try:
                    photon[E] += 1
                except KeyError:
                    photon[E] = 1
                        
            size = particles.ELECTRONS.size()
            for _ in range(size):
                E = particles.ELECTRONS.back()
                particles.ELECTRONS.pop_back()
                try:
                    electron[E] += 1
                except KeyError:
                    electron[E] = 1
                    
        en = []
        ph = []
        for key, value in photon.items():
            en.append(key)
            ph.append(value)
        
        en = np.array(en)
        ph = np.array(ph)
        
        permute = en.argsort()
        en = en[permute]
        ph = ph[permute]
        
        en2 = []
        el = []
        
        for key, value in electron.items():
            en2.append(key)
            el.append(value)        
        en2 = np.array(en2)
        el = np.array(el)
        
        permute = en2.argsort()
        en2 = en2[permute]
        el = el[permute]
        
                
        tf = time.perf_counter_ns()
        return en, ph, en2, el
            
        
        print("SAMPLE SIZE:", size)
        
        print(f"{(tf-t0)/N}ns per run")
        
        
        
        return
        # cdef vector[SEC] res = deref(self.rATOM.run(shell_index, CUT_OFF))
        
        # cdef int N= res.size()
        # cdef int i
        # to_return = []
        # for i in range(N):
        #     to_return += [(res[i].TYPE, res[i].E)]
            
        # return to_return


        
        #for shell in self.SHELLS:
            #shell.setTransitions(self)
            
        #self.CS = db.EPDL[Z-1][(7, 73, 0, 0, 0, 0)].getLinLinInterpol()
        

        
        

    
    def __repr__(self):
        rep = f"<Atom Z = {self.Z}, Aw = {self.Aw}, Nsh = {self.Nsh}> \n"
        
        if self.Nsh > 5:
            rep += f"    <Shell #{0} #el = {self.number_el[0]}, E = {self.BE[0]}> \n"
            rep += f"    <Shell #{1} #el = {self.number_el[1]}, E = {self.BE[1]}> \n"
            rep += f"        ... + {self.Nsh - 4} shells +  ...       \n"
            rep += f"    <Shell #{self.Nsh-2} #el = {self.number_el[-2]}, E = {self.BE[-2]}> \n"
            rep += f"    <Shell #{self.Nsh-1} #el = {self.number_el[-1]}, E = {self.BE[-1]}> \n"
            return rep
        
        for i in range(self.Nsh):
            rep += f"    <Shell #{i} #el = {self.number_el[i]}, E = {self.BE[i]}> \n"
        
        
        
        
        return rep
    

        
   #  cdef void _printShells(self):
   #      cdef Shell shell
   #      for shell in self.SHELLS:
   #          print("nVAC:", shell.Nvac)
   #          print(f"------{shell.i}--------|| ", shell.Nvac*"* | " + int(shell.Nel - shell.Nvac)*"e | ")
   #          #print("#VACANCIES = ", shell.Nvac)
   #          #print("#ELECTRONS = ", shell.Nel)
    
   #  def printShells(self):
   #      self._printShells()




   #  cdef void _introduceVacancy(Atom self, Shell shell):
   #      cdef int i
   #      cdef Vacancy newVac = Vacancy._new(shell)
   #      cdef Vacancy oldVac
   #      for i in range(len(self.vacancies)):
   #          oldVac = self.vacancies[i]
   #          if oldVac.current_shell.i < newVac.current_shell.i:
   #              self.vacancies.insert(i, newVac)
   #              break
   #      else:
   #          self.vacancies.append(newVac)
        
    
   #  def ionize(self, Shell shell):
   #      self._introduceVacancy(shell)
    
   #  def introduceVacancy(self, int i):
   #      cdef Shell shell
   #      for shell in self.SHELLS:
   #          if shell.i == i:
   #              self.vacancies.append(Vacancy._new(shell))
                
   #  # cdef void _introduceVacancy(self, int i):
   #  #     cdef Shell shell
   #  #     for shell in self.SHELLS:
   #  #         if shell.i == i:
   #  #             self.vacancies.append(Vacancy._new(shell))
    
   #  cpdef printVacancies(self):
   #      print(self.vacancies)
    
    
   #  def __getitem__(self, int i):
   #      cdef Shell shell
   #      for shell in self.SHELLS:
   #          if shell.i == i:
   #              return shell
   #      return None
        
    
   # # def run(self):
   # #     return self._run()
    
        
        
        
   #  cdef void reset(self):
   #      cdef Shell shell
   #      self.vacancies = []
   #      for shell in self.SHELLS:
   #          shell.Nvac = 0
    
    
    cpdef getBookmarkedText(self):
        cdef str path = self.path

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
    
    



























    # cdef list _run(self, Volume current_region, double x, double y, double z):
    #     self.nSECONDARY = 0
    #     #cdef dict spectrum = {'el':[], 'ph':[]}
    #     cdef list secondary = []
    #     cdef Vacancy activeVac
    #     cdef list possible_transitions
    #     cdef Shell shell_i, shell_k, shell_j
        
    
    #     cdef object probs
        
    #     cdef list cumul 
    #     cdef double C, P
        
    #     cdef int _, i
    #     cdef Shell j, k
    #     cdef double E, p2, p1
    #     cdef double[:, :] possible_transitionsArr
        
    #     cdef Transition transition, CHOSEN_TRANSITION
        
    #     while len(self.vacancies) > 0:
    #         #self.vacancies.sort()
            
    #         #print("___________________________________________________")
    #         #self.printShells()

    #         activeVac = self.vacancies[0]
            
    #         if activeVac.current_shell.i > 15:
    #             break
            
    #         if activeVac.current_shell.trans == []:
    #             self.vacancies.remove(activeVac)
    #             continue
                
            
    #         possible_transitions = []
    #         #print(activeVac.current_shell.trans)
    #         for transition in activeVac.current_shell.trans:
    #             #print(transition)
    #             if transition.k is False:
    #                 if transition.j.Nel - transition.j.Nvac > 0:
    #                     possible_transitions.append(transition)
                        
    #             elif transition.j == transition.k:
    #                 if transition.j.Nel - transition.j.Nvac > 1:
    #                     possible_transitions.append(transition)
    #             else:
    #                 if  transition.j.Nel - transition.j.Nvac > 0 \
    #                 and transition.k.Nel - transition.k.Nvac > 0:
    #                     possible_transitions.append(transition)

                    
                
    #         #print(possible_transitions)
    #         if possible_transitions == []:
    #             self.vacancies.remove(activeVac)
    #             continue
            
    #         if len(possible_transitions) == 1:
    #             CHOSEN_TRANSITION = possible_transitions[0]
                
    #         else:
                
    #             C = 0.
    #             for transition in possible_transitions:
    #                 C += transition.p
                
    #             probs = [transition.p/C for transition in possible_transitions]
                
    #             cumul = [0.]
    #             C = 0.
    
    #             for p in probs:
    #                 C += p
    #                 cumul += [C]
                

    #             #cumul = [sum(probs[0:i]) for i in range(len(probs))]
    #             for _ in range(100_000):
    #                 i = search._sortedListDOUBLE(cumul, np.random.rand(), 0, len(cumul)-1)
                    
    #                 transition = possible_transitions[i-1]
                    
    #                 #j, k, _, E = possible_transitions[i-1, :]
    #                 #j, k, E = activeVac.sampleTransition()
                    
    #                 p2 = (transition.k.Nel - transition.k.Nvac) / transition.k.Nel \
    #                      if transition.k is not False else 1
                    
    #                 p1 = (transition.j.Nel - transition.j.Nvac) / transition.j.Nel
    
    #                 if np.random.rand() < p1*p2:
    #                     CHOSEN_TRANSITION = transition
    #                     break
    #             else:
    #                 print(possible_transitions)
    #                 print(activeVac.current_shell.i)
    #                 print("p1*p2 = ", p1p2)
    #                 raise RuntimeError("exceeded number of iter")
                
    #         #if j is None: #or j > 45:
    #         #    self.vacancies.remove(activeVac)
    #         #print(CHOSEN_TRANSITION)
            
            
    #         if CHOSEN_TRANSITION.RADIATIVE is True:
    #             if CHOSEN_TRANSITION.E > photonCUTOFF:
                    
    #                 secondary.append(Photon._newISOTROPIC(CHOSEN_TRANSITION.E, 
    #                                                       x, y, z,
    #                                                       current_region))
    #                 self.nSECONDARY += 1
  
    #             activeVac.move(CHOSEN_TRANSITION.j)
  
    #         else:
    #             activeVac.move(CHOSEN_TRANSITION.j)
    #             self._introduceVacancy(CHOSEN_TRANSITION.k)
                
    #             if CHOSEN_TRANSITION.E > electronCUTOFF:
    #                 self.nSECONDARY += 1
    #                 secondary.append(Electron._newISOTROPIC(CHOSEN_TRANSITION.E*1e6, 
    #                                                       x, y, z,
    #                                                       current_region))
  
                    
    #             #spectrum['el'] += [CHOSEN_TRANSITION.E]
                

    #     #self.printShells()                
    #     self.reset()
    #     return secondary



















# cdef class Vacancy:

    
#     @staticmethod
#     def new(Shell current_shell):
#         return Vacancy._new(current_shell)
    
#     @staticmethod
#     cdef Vacancy _new(Shell current_shell):
#         self = <Vacancy>Vacancy.__new__(Vacancy)
        
#         self.current_shell = current_shell
#         #self.current_shell.Nel -= 1
#         self.current_shell.Nvac += 1
#         return self
        
    
    
    

         
#     @property
#     def __hash__(Vacancy self):
#         return self.current_shell.i
    
    
#     cdef void move(Vacancy self, Shell shell):
#         #self.current_shell.Nel += 1
        
#         cdef Shell current_shell = self.current_shell
        
#         current_shell.Nvac -= 1
#         #shell.Nel -= 1
#         shell.Nvac += 1
        
#         self.current_shell = shell


# cdef class Shell:
#     shells[i], number_el[i],  binding_energy[i],  
#     def __init__(self, i, fi, Ui, rTrans, nrTrans):
#         self.rSHELL = Shell(i, fi, Ui)
        
        
        
#         #self.Nvac = 0
        
#         #self.args = args
#         #self.i, self.Nel, self.binding_energy, self.KE, self.avgR, self.Z = args
#         #self.CS = db.EPDL[self.Z-1][(7, 73, 91, self.i, 0, 0)].getLinLinInterpol()

#         #self.kwargs = kwargs
#         self.rTrans = rTrans
#         self.nrTrans = nrTrans
        
#         if (self.rTrans  is not None) and \
#            (self.nrTrans is not None):
            
#             self.rTrans = np.insert(self.rTrans, 1, 0, axis = 1)

#             self._trans = append(self.rTrans, self.nrTrans)
#             self._trans = concatenate((self.rTrans, self.nrTrans), axis = 0)
#             self.transFLAG = True
#             self.rSHELL.transFLAG = True
#         else:
#             self.transFLAG = False
#             self.rSHELL.transFLAG = False

 
        
#     cpdef void setTransitions(self, Atom atom):
        
#         if self.transFLAG is False:
#             return
        
        
        
#         cdef list TRANS = []
#         cdef double[:] trans
#         cdef int i, k
#         cdef double p, E
#         cdef rShell* shell
#         cdef rShell* shell1, shell2
        
#         for trans in self._trans:
#             j, k, p, E = trans[0], int(trans[1]), trans[2], trans[3]
#             if k == 0:
#                 shell = atom[j]
#                 self.rSHELL.newTransition(shell, shell, p, E, True)
                
#                 #TRANS.append(Transition._new(shell, shell, p, E, True))
#             else:
                
#                 self.rSHELL.newTransition(atom[j], atom[k], p, E, True)
                
#                 #TRANS.append(Transition._new(atom[j], atom[k], p, E, False))
            
#         #self.trans = TRANS
    
    

  