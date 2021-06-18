# cython: profile=False


from libc.math cimport log10,pi, log
import numpy as np
from ...tools.interpol1 cimport LogLinInterpolation

cimport cython

from numpy import array, logspace



cdef double  ELECTRON_REST_MASS      = 0.51099895000e6            
cdef double  _2ELECTRON_REST_MASS    = 2 *ELECTRON_REST_MASS
cdef double  _4ELECTRON_REST_MASS_2  = 4*ELECTRON_REST_MASS**2
cdef double  SPEED_OF_LIGHT          = 2.99792458e10
cdef double  H_BAR                   = 6.5821e-16 
cdef double  MASS_ELECTRON           = 9.1094e-28 
cdef double  ELECTRON_CHARGE         = 4.8032e-10

ctypedef double adimensional

cdef double CONST = 2*pi *ELECTRON_CHARGE**4 / MASS_ELECTRON  *6.242e+11

cdef double Wcc = 0.



from libc.stdlib cimport rand, RAND_MAX


cimport cython
@cython.cdivision(True)
cdef double urand():
    cdef double r = rand()
    return r / (RAND_MAX + 1)

#from numpy.random import rand as urand




cdef struct dynS:
    double E, X, v2, beta2, delta

cdef dynS STATE




cdef class EmptyInterp(LogLinInterpolation):
    def __init__(self):
        pass
    
    cdef double _eval(self, double x):
        return 0




@cython.initializedcheck(False)
@cython.cdivision(True)
cdef class CMolecule:





    def __init__(self, object pyMolecule, object cb, object delta):
        self.SIGMA0 = []
        cdef list ATOMS = []
        
        STATE.E = 0.
        STATE.X = 0.
        STATE.beta2 = 0.
        STATE.delta = 0.
        STATE.v2 = 0.
        
        
        print("        > creating shells")
        cdef object atom
        for atom in pyMolecule:
            ATOMS.append(CAtom(atom))
            
        self.N = len(ATOMS)
        self.number_density = pyMolecule.N
        
        self.ATOMS = array(ATOMS)
        
        if cb.empty == True:
            self.cb = EmptyShell()
            self.delta = EmptyInterp()
            return
        else: self.cb = CShell(cb)
        
        print("        > making density corrections")
        def delta2(double E):
            X = (E/ELECTRON_REST_MASS + 1)**2
            beta2 = (X-1)/X
            return delta(beta2)
        
        print("        > note: KE = logspace(1, 9, 5_000)")
        KE = logspace(1, 9, 5_000)
        dd = array([delta2(E) for E in KE])

        import numpy
        self.delta = LogLinInterpolation(numpy.log10(KE), dd)
        

    def __call__(self, i): return self.ATOMS[i]
        
    # def py_sample(self):
    #     return self.sample()
    
    cdef void reset(self):
        self.sigma0, self.sigma1, self.sigma2 = 0, 0, 0
        
    cdef void add(self, double s0, double s1, double s2):
        self.sigma0 += s0
        self.sigma1 += s1 
        self.sigma2 += s2
    
    
    
    @cython.initializedcheck(False)
    @cython.boundscheck(True)
    @cython.wraparound(False)
    cdef void update(CMolecule self, double E, bint record):
        
        
        self.reset()
        
        
        STATE.E = E
        STATE.X = (E/ELECTRON_REST_MASS + 1)**2
        STATE.beta2 = (STATE.X-1)/STATE.X
        STATE.v2 = SPEED_OF_LIGHT**2 * STATE.beta2
        STATE.delta = self.delta._eval(E) #self.delta(beta2)
        #STATE.delta = 0
        
        
        cdef CAtom atom
        cdef int i
        for i in range(self.N):
            atom = self.ATOMS[i]
            atom.update(record)
            self.add(atom.sigma0, atom.sigma1, atom.sigma2)
        
        self.cb.update(record)
        self.add(self.cb.sigma0, self.cb.sigma1, self.cb.sigma2)
        self.imfp = self.number_density*self.sigma0
        self.SP = self.number_density*self.sigma1
        self.stragg = self.number_density*self.sigma2
        if record:
            self.SIGMA0.append(self.sigma0)
        
        
        
        
        

        
        
        
        
        
        
    # cdef (double, double, double, double) sample(CMolecule self, Volume current_region, double x, double y, double z):
        
    #     self.secondary = []
    #     self.nSECONDARY = 0
        
    #     cdef double cumul = self.cb.sigma0
    #     cdef double r = urand()*self.sigma0
        
        
    #     if cumul > r:
    #         return self.cb.sample()
        
        
    #     cdef CAtom atom
    #     cdef int i
    #     for i in range(self.N):
    #         atom = self.ATOMS[i]
    #         cumul += atom.sigma0
    #         if cumul > r:
    #             break
        

        
    #     atom.sample(current_region, x, y, z)
        
    #     self.secondary += atom.secondary
    #     self.nSECONDARY += atom.nSECONDARY
        
        
    #     return atom.last_sample
    
    def cut(self, Wcc):
        cdef CAtom atom
        
        new_atoms = []

        
        cdef int i
        for i in range(self.N):
            atom = self.ATOMS[i]
            if not atom.cut(Wcc):
                new_atoms.append(atom)

        
        if new_atoms == []: new_atoms = [EmptyAtom()]
            
        
        
        self.ATOMS = array(new_atoms)
        self.N = len(self.ATOMS)
        
        print("")
        print("")
        print("SHELLS LEFT IN GOS MODEL:")
        for i in range(self.N):
            atom = self.ATOMS[i]
            print(atom)

        
        if self.cb.Wk < Wcc:

            self.cb = EmptyShell()
            
            print("NO CB SHELL")
        else:
            print(self.cb)
            

        
    def py_update(self, double E):
        self.update(E, False)
    
    def generate_tables(self, *args, select = "total"):
        Eaxis = logspace(*args)
        sigma0, sigma1, sigma2 = [], [], []
        cdef double E
        if select == "total":
            for E in Eaxis:
                self.update(E, False)
                sigma0.append(self.sigma0)
                sigma1.append(self.sigma1)
                sigma2.append(self.sigma2)
            sigma0, sigma1, sigma2 = map(array, (sigma0, sigma1, sigma2) )
            return Eaxis, sigma0, sigma1, sigma2
        cdef CAtom atom
        if select == "atom":
            atom = self.ATOMS[0]
            for E in Eaxis:
                self.update(E, False)
                sigma0.append(atom.sigma0)
                sigma1.append(atom.sigma1)
                sigma2.append(atom.sigma2)
            sigma0, sigma1, sigma2 = map(array, (sigma0, sigma1, sigma2) )
            return Eaxis, sigma0, sigma1, sigma2
        
        cdef CShell shell
        if select == "K":
            shell = self.ATOMS[0].SHELLS[0]
            for E in Eaxis:
                self.update(E, False)
                sigma0.append(shell.sigma0)
                sigma1.append(shell.sigma1)
                sigma2.append(shell.sigma2)
            sigma0, sigma1, sigma2 = map(array, (sigma0, sigma1, sigma2) )
            return Eaxis, sigma0, sigma1, sigma2, shell.Wk
        
        if select == "L":
            shell = self.ATOMS[0].SHELLS[1]
            for E in Eaxis:
                self.update(E, False)
                sigma0.append(shell.sigma0)
                sigma1.append(shell.sigma1)
                sigma2.append(shell.sigma2)
            sigma0, sigma1, sigma2 = map(array, (sigma0, sigma1, sigma2) )
            return Eaxis, sigma0, sigma1, sigma2, shell.Wk
        
        if select == "M":
            shell = self.ATOMS[0].SHELLS[2]
            for E in Eaxis:
                self.update(E, False)
                sigma0.append(shell.sigma0)
                sigma1.append(shell.sigma1)
                sigma2.append(shell.sigma2)
            sigma0, sigma1, sigma2 = map(array, (sigma0, sigma1, sigma2) )
            return Eaxis, sigma0, sigma1, sigma2, shell.Wk
    
@cython.initializedcheck(False)
@cython.cdivision(True)            
cdef class CAtom:

    
    def __init__(self, object pyAtom):
        self.SIGMA0 = []

        #self.rAtom = pyAtom.rAtom
        
        self.Z = pyAtom.Z
        
        SHELLS = []
        for shell in pyAtom:
            SHELLS.append(CShell(shell))

        
        
        self.SHELLS = array(SHELLS)
        self.N = len(self.SHELLS)
        self.x = pyAtom.x
        cdef CShell shell1
        
        #outer shell does not relax
        if self.N >= 4:
        
            shell1 = SHELLS[3]
            shell1.close.relax = False
            shell1.distant.relax = False
        
    cdef void reset(self):
        self.sigma0, self.sigma1, self.sigma2 = 0, 0, 0
        
    cdef void add(self, double s0, double s1, double s2):
        self.sigma0 += s0
        self.sigma1 += s1 
        self.sigma2 += s2
        
    @cython.initializedcheck(False)
    @cython.boundscheck(True)
    @cython.wraparound(False)
    cdef void update(self, bint record):
        cdef CShell shell
        self.reset()
        
        cdef int i
        for i in range(self.N):
            shell = self.SHELLS[i]
            shell.update(record)
            self.add(shell.sigma0, shell.sigma1, shell.sigma2)
            
        self.sigma0 *= self.x
        self.sigma1 *= self.x
        self.sigma2 *= self.x
        if record:
            self.SIGMA0.append(self.sigma0)
        
        

        
        
        
        
    def __repr__(self):
        
        rep = ""
        cdef CShell shell
        for shell in self.SHELLS:
            rep += "\n" + shell.__repr__()
        return rep
        
    # cdef void sample(CAtom self, Volume current_region, double x, double y, double z):
    #     """
    #     Choose shell, ionize it using the pyRelax computational model and call the shells sampling method.
    #     """
        
    #     # CHOOSING SHELL 
    #     cdef double cumul = 0.
        
    #     cdef double r = urand()*self.sigma0
        
    #     cdef int i
        
    #     cdef CShell shell
        
    #     for i in range(self.N):
    #         shell = self.SHELLS[i]
    #         cumul += shell.sigma0
    #         if r < cumul:
    #             break
            
            
        
    #     # IONIZING SHELL
    #     #shell.sample()
        
    #     self.rAtom.introduceVacancy(shell.choose_shell())
    #     self.secondary = self.rAtom._run(current_region, x, y, z)
    #     self.nSECONDARY = self.rAtom.nSECONDARY
        
        
        
    #     # SAMPLING FROM SHELL
    #     self.last_sample = shell.sample()

    def cut(self, Wcc):
        cdef CShell shell
        
        new_shells = []

        cdef int i
        for i in range(self.N):
            shell = self.SHELLS[i]
            if shell.Wk >= Wcc:
                new_shells.append(shell)

        
        

        if new_shells == []:
            return True
                
        
        
        self.SHELLS = array(new_shells)
        self.N = len(self.SHELLS)
        return False


    def __call__(self, i): return self.SHELLS[i]
    

cdef class EmptyAtom(CAtom):
    def __init__(self):
        self.sigma0 = 0
        self.sigma1 = 0
        self.sigma2 = 0
        

        
    def __repr__(self):
        return "<EMPTY ATOM>"
    cdef void update(self, bint record):
        pass

cdef class EmptyShell(CShell):
    def __init__(self):
        self.SIGMA0 = []
        self.sigma0 = 0
        self.sigma1 = 0
        self.sigma2 = 0
        
        self.Wk = 0
        self.Uk = 0
        self.fk = 0
        
        
        self.distant = Distant._new  (0, 0, 0)
        
        
        
        
        self.close   = Close._new(0, 0, 0)
        
        
    
    cdef void update(self, bint record):
        
        if record:
            self.close.SIGMA0.append(0)
            self.distant.SIGMA0T.append(0)
            self.distant.SIGMA0L.append(0)

            self.SIGMA0.append(0)
        

@cython.initializedcheck(False)
@cython.cdivision(True)
cdef class CShell:


    

    def __init__(self, object pyShell):
        
        
        self.designators = pyShell.designators
        self.Ncollection = len(self.designators)
        
        self.SIGMA0 = []
        self.INDEX = pyShell.INDEX
        self.BE = pyShell.BE
        
        
        self.Wk = pyShell.Wk
        self.Uk = pyShell.Uk
        self.fk = pyShell.fk
        self.relax = True
        
        self.distant = Distant._new  (self.fk, self.Wk, self.Uk)
        self.close   = Close._new(self.fk, self.Wk, self.Uk)
        
    def __repr__(self):
        return f"<Shell{self.designators} Wk = {self.Wk}, Uk = {self.Uk}, fk = {self.fk}>"
        
        

        
        
    cdef void setRelax(self, bint relax):
        self.relax = relax
        self.close.relax = relax
        self.distant.relax = relax
        
    cdef void update(self, bint record):
        
        #cdef double X = (E/ELECTRON_REST_MASS + 1)**2
        #cdef double beta2 = (X-1)/X                              
        #cdef double v2 = SPEED_OF_LIGHT**2 * beta2
        
        
        cdef double const = CONST / STATE.v2 * self.fk
        #const = const #*6.242e+11*1e+24 #to barn*eV
        
        self.close.update(const, record)
        self.distant.update(const, record)
        
        self.sigma0 = self.close.sigma0 + self.distant.sigma0
        self.sigma1 = self.close.sigma1 + self.distant.sigma1
        self.sigma2 = self.close.sigma2 + self.distant.sigma2
        if record:
            self.SIGMA0.append(self.sigma0)
    
    # #(newE, cos, Esec, cos_sec)
    # cdef (double, double, double, double) sample(CShell self):
        
    #     cdef double cumul = self.close.sigma0
    #     cdef double r = urand()*self.sigma0
        
    #     if cumul > r:
    #         return self.close.sample()
        
    #     cumul += self.distant.L.sigma0
    #     if cumul > r:
    #         return self.distant.sampleL()
        
    #     return self.distant.sampleT()
    
    cdef int choose_shell(self):

 
        return self.designators[<int> (urand()*self.Ncollection)]
    

        
        
        
@cython.initializedcheck(False)
@cython.cdivision(True)       
cdef class Close:

    
    @staticmethod
    cdef Close _new(double fk, double Wk, double Uk):
        self = <Close>Close.__new__(Close)
        self.fk = fk
        self.Wk = Wk
        self.Uk = Uk
        self.SIGMA0 = []
        return self
        
    @cython.cdivision(True)      
    cdef void update(Close self, double const, bint record):
        
        
        if STATE.E/2 < self.Wk:
            self.sigma0 = 0
            self.sigma1 = 0
            self.sigma2 = 0
            if record:
                self.SIGMA0.append(0.)
            return
        
        cdef double s0o, s1o, s2o, s0f, s1f, s2f
        s0o, s1o, s2o = Close.J(self.Wk)
        s0f, s1f, s2f = Close.J(STATE.E/2)
        
        
        self.sigma0 = const*(s0f - s0o)
        self.sigma1 = const*(s1f - s1o)
        self.sigma2 = const*(s2f - s2o)
        
        if record:
            self.SIGMA0.append(self.sigma0)
        

    @staticmethod    
    cdef (double, double, double) J(double W):
        cdef double a, J0, J1, J2
        
        a = (STATE.E/(STATE.E + ELECTRON_REST_MASS))**2
        cdef double diff = STATE.E - W
        cdef double logW = log(W)
        cdef double logdiff = log(diff)
    
        J0 = -1/W + 1/diff + (1 - a)/STATE.E * log(diff/W) + a * W/ STATE.E**2
        
        J1 = log(W) + STATE.E /(diff) + (2-a)*log(diff) + a*W**2 / 2 / STATE.E**2
        
        J2 = (2-a)*W + (2*STATE.E**2 - W**2)/(diff) + (3-a)*STATE.E*log(diff) + a*W**3/3/STATE.E**2
    
        return J0, J1, J2
    
    # #(newE, cos, Esec, cos_sec)    
    # cdef (double, double, double, double) sample(Close self):
    #     cdef double gamma = STATE.X
    #     cdef double a = (STATE.E/(STATE.E + _2ELECTRON_REST_MASS))**2
        
    #     cdef double gp1 = (gamma + 1)**2
    #     self.b1 = a*(2*gp1 - 1 )/(gamma**2 - 1 )
    #     self.b2 = a*(3*gp1 + 1)/gp1
    #     self.b3 = a*2*gamma*(gamma - 1)/gp1
    #     self.b4 = a*(gamma - 1)**2 / gp1
        
        
        
        
        
    #     self.kc = max(self.Wk, Wcc)/STATE.E
    #     cdef double z
    #     cdef double k, k2, factor
        
    #     while 1:
    #         z = (1 + 5*a*self.kc*.5)*urand()
    #         if z < 1: k = self.kc / (z* (1 - 2*self.kc) )
    #         else :    k = self.kc + (z - 1)*(1 - 2*self.kc)/(5*a*self.kc)

    #         k2 = k**2
    #         factor = (1+5*a*k2)

    #         if urand()*factor < self.P(k)*k2:
    #             break
                
    #     cdef double W
    #     cdef double cos2
        
    #     W = k*STATE.E
    #     cos2 = (STATE.E - W)*(STATE.E + _2ELECTRON_REST_MASS)/STATE.E/(STATE.E - W + _2ELECTRON_REST_MASS)
    #     cdef double cos = cos2**.5
        
    #     #SECONDARY ELECTRON
    #     cos2 = W * (STATE.E + 2*0.511e6)/STATE.E/(W + 2*0.511e6)
    #     cdef double cos_sec = cos2**.5
    #     if self.relax:
    #         return (STATE.E - W, cos, W-self.Uk, cos_sec)         
    #     return (STATE.E - W, cos, W, cos_sec)
    
    cdef double P(self, double k):
        if k < self.kc:
            return 0.
        if k > 1:
            return 0.
        return 1/k**2 - self.b1 / k + self.b2 - self.b3*k + self.b4 *k**2
        
      
cdef struct SIGMA:
    double sigma0, sigma1, sigma2
    
    
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef class Distant:

    
    @staticmethod
    cdef Distant _new(double fk, double Wk, double Uk):
        self = <Distant>Distant.__new__(Distant)
        self.Wk = Wk
        self.fk = fk
        self.Uk = Uk

        self.Ldef.sigma0 = 0.; self.Ldef.sigma1 = 0.; self.Ldef.sigma2 = 0.
        self.Tdef.sigma0 = 0.; self.Tdef.sigma1 = 0.; self.Tdef.sigma2 = 0.
        
        self.SIGMA0T = []
        self.SIGMA0L = []
        return self
        
    cdef void update(self, double const, bint record):
        
        if STATE.E < self.Wk:
            self.L = self.Ldef
            self.T = self.Tdef
            if record:
                self.SIGMA0L.append(0.)
                self.SIGMA0T.append(0.)
            return
        
        cdef double Wk = self.Wk
  
        self.Qm = self.Q_minus(STATE.E)
        cdef double Q_minus = self.Qm


        cdef double logval = log( Wk  * (Q_minus + _2ELECTRON_REST_MASS)/(  Q_minus*(Wk + _2ELECTRON_REST_MASS)  ))
        cdef double common_mul = const*logval
        
        self.L.sigma0 = common_mul/Wk
        self.L.sigma1 = common_mul
        self.L.sigma2 = common_mul*Wk



        logval = log(STATE.X) - STATE.beta2 - STATE.delta
        common_mul = const*logval
        self.T.sigma0 = common_mul/Wk            #calculation of imfp
        self.T.sigma1 = common_mul               #calculation of stopping power
        self.T.sigma2 = common_mul*Wk            #calculation of energy straggling        
        
        self.sigma0 = self.T.sigma0 + self.L.sigma0
        self.sigma1 = self.T.sigma1 + self.L.sigma1
        self.sigma2 = self.T.sigma2 + self.L.sigma2
        
        
        if record:
            self.SIGMA0L.append(self.L.sigma0)
            self.SIGMA0T.append(self.T.sigma0)
                
    cdef double Q_minus(self, double E):
        """eq 3.83"""
        cdef double dE = E - self.Wk
        return ( (self.momentum(E) - self.momentum(dE))**2  + ELECTRON_REST_MASS**2  )**.5 \
               - ELECTRON_REST_MASS 

    cdef double momentum(self, double E):
        return (  E * (E + _2ELECTRON_REST_MASS)  ) ** .5
    
    cdef double p2(self, double E):
        return E * (E + _2ELECTRON_REST_MASS) 
    
    #(newE, cos, Esec, cos_sec)
    # cdef (double, double, double, double) sampleT(Distant self):
    #     cdef double Esec = self.Wk - self.Uk if self.relax else self.Wk
    #     cdef double newE = STATE.E - self.Wk
    #     return (newE, 1., Esec, -1.)
        
    
    #(newE, cos, Esec, cos_sec)
    # cdef (double, double, double, double) sampleL(Distant self):
    #     cdef double newE = STATE.E - self.Wk 
    #     cdef double Esec = self.Wk - self.Uk if self.relax else self.Wk
        
    #     cdef double Qs = self.Qm/(1 + self.Qm / _2ELECTRON_REST_MASS)
        
        
    #     cdef double r = urand()
    #     cdef double A = Qs/self.Wk *(1 + self.Wk/_2ELECTRON_REST_MASS )
    #     cdef double B = Qs/_2ELECTRON_REST_MASS
    #     cdef double Q = Qs / (A**r - B)
        
    #     cdef double p2E = self.p2(STATE.E)
    #     cdef double p2d = self.p2(STATE.E - self.Wk)
    #     cdef double p2Q = self.p2(Q)
        
    #     cdef double cos = p2E + p2d - p2Q
    #     cos = .5*cos/(p2E*p2d)**.5
        
    #     cdef double Wk2 = self.Wk**2
    #     A = Wk2 / STATE.beta2 / p2Q
    #     B = (p2Q - Wk2)/(2*self.Wk*(STATE.E + _2ELECTRON_REST_MASS))
    #     cdef double cos_sec = A * (1 + B) **2
        
        
    #     return (newE, cos, Esec, cos_sec)

        
        
