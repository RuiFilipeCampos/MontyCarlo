# cython: profile=True

__doc__ = """
    sources.pyx
        This module is VERY incomplete. Sources decide the flow of the simulation, they get to choose which primary particles get created, in which order
        they are simulated and even will be responsible for seeding them for multi processing.
"""

__author__ = "Rui Campos"


print("Importing `.sources`")


from .particles.particle cimport STATE

from .particles.positrons cimport Positron
from .particles.photons cimport Photon
from .particles.electrons cimport Electron
from .particles.particle cimport Particle

from .geometry.main cimport Volume

from collections import deque


from .external.mixmax_interface cimport mixmax_engine


from numpy.random import randint

cdef long int SEED = randint(100_000)
cdef mixmax_engine gen = mixmax_engine(0,0,0,SEED);

cdef mixmax_engine *genPTR = &gen;



import pyvista as pv
from libc.stdlib cimport srand
from numpy.random import randint

seed = randint(1e9)
import numpy as np
print("SEED:", seed)
srand(seed)
cimport numpy as cnp






cdef class Source:
    cdef mixmax_engine gen
    cdef mixmax_engine *genPTR
    cdef object container
    cdef int N

    cdef dict __dict__

    cpdef void runNOsec(self):

        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        


        while nCURRENT < self.N:
            print(f"{nCURRENT}/{nSIMULATE}")
 
            p = self.container.popleft()

            p._run(self.genPTR)

            self.container.extend(p.secondary)
            
            nCURRENT += 1


    cdef void _run(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        


        while nCURRENT < nSIMULATE:
            #print(f"{nCURRENT}/{nSIMULATE}")
 
            p = self.container.popleft()

            p._run(self.genPTR)

            self.container.extend(p.secondary)
            
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1

    cdef void _runPOS(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0

        self.pos_record = deque()

        self.pos_record.append(deque())
        self.pos_record.append(deque())
        self.pos_record.append(deque())


        while nCURRENT < nSIMULATE:
            print(f"{nCURRENT}/{nSIMULATE}")
 
            p = self.container.popleft()

            p._run(self.genPTR)


            if isinstance(p, Photon):
                self.pos_record[0].append(p.get_record_pos())
            elif isinstance(p, Electron):
                self.pos_record[1].append(p.get_record_pos())
            else:
                self.pos_record[2].append(p.get_record_pos())

            self.container.extend(p.secondary)
            
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1


    cpdef run(self, record_pos = True):
        if record_pos:
            self._runPOS()
        else: self._run()


    



cdef class IsotropicPoint(Source):
    cdef object particle, pos

    def __init__(self, particle, Volume region, double E = 1e6,  int N = 1000, pos = (0, 0, 0)):
        
        valid_options = ["photon", "positron", "electron"]

        assert particle in valid_options
        assert len(pos) == 3

        self.N = N
        self.pos = np.array(pos)

        self.particle = particle
        self.container = deque()



        cdef STATE state
        state.pos.x = self.pos[0]
        state.pos.y = self.pos[1]
        state.pos.z = self.pos[2]

        if not region.is_inside(state.pos):
            raise RuntimeError("Particles will not be inside the starting region.")
            import time
            time.sleep(1000)

        state.current_region = <void*> region

        self.gen = mixmax_engine(0, 0, 0, SEED);
        self.genPTR = &self.gen
        state.genPTR = self.genPTR


        state.E = E

        cdef int i

        if particle == "photon":
            for i in range(self.N):
                self.container.append(Photon._newISOTROPIC(state))
        elif particle == "electron":
            for i in range(self.N):
                self.container.append(Electron._newISOTROPIC(state))
        else:
            for i in range(self.N):
                self.container.append(Positron._newISOTROPIC(state))




cdef class UnitTestParticle:
    cdef dict __dict__
    cdef mixmax_engine gen
    cdef mixmax_engine* genPTR

    def __init__(self, particle, Volume region, double E = 1e6,  pos = (0, 0, 0), dire = (0, 0, 1)):

        valid_options = ["photon", "positron", "electron"]

        assert particle in valid_options
        assert len(pos) == 3 and len(dire) == 3


        self.pos = np.array(pos)
        self.dire = np.array(dire)
        self.dire = self.dire/np.sqrt(sum(self.dire**2))

        self.particle = particle
        self.container = deque()
        cdef long int SEED = randint(100_000)
        self.gen = mixmax_engine(0, 0, 0, SEED);
        self.genPTR = &gen;

        cdef STATE state
        state.pos.x = self.pos[0]
        state.pos.y = self.pos[1]
        state.pos.z = self.pos[2]

        if not region.is_inside(state.pos):
            raise RuntimeError("Particles will not be inside the starting region.")
            import time
            time.sleep(1000)

        state.dire.x = self.dire[0]
        state.dire.y = self.dire[1]
        state.dire.z = self.dire[2]

        state.genPTR = self.genPTR
        state.current_region = <void*> region

        state.E = E

        cdef int i
        cdef Particle p

        if particle == "photon":
            p = Photon._new(state)
            p.throwAZIMUTH()
            self.container.append(p)
        elif particle == "electron":
            p = Electron._new(state)
            p.throwAZIMUTH()
            self.container.append(p)
        else:
            p = Positron._new(state)
            p.throwAZIMUTH()
            self.container.append(p)



cdef class Beam(Source):
    cdef object particle, pos, dire

    def __init__(self, particle, Volume region, double E = 1e6,  int N = 1000, pos = (0, 0, 0), dire = (0, 0, 1)):
        
        valid_options = ["photon", "positron", "electron"]

        assert particle in valid_options
        assert len(pos) == 3 and len(dire) == 3
        
        self.N = N
        self.pos = np.array(pos)
        self.dire = np.array(dire)
        self.dire = self.dire/np.sqrt(sum(self.dire**2))

        self.particle = particle
        self.container = deque()
        cdef long int SEED = randint(100_000)
        self.gen = mixmax_engine(0, 0, 0, SEED);
        self.genPTR = &gen;

        cdef STATE state
        state.pos.x = self.pos[0]
        state.pos.y = self.pos[1]
        state.pos.z = self.pos[2]

        if not region.is_inside(state.pos):
            raise RuntimeError("Particles will not be inside the starting region.")
            import time
            time.sleep(1000)

        state.dire.x = self.dire[0]
        state.dire.y = self.dire[1]
        state.dire.z = self.dire[2]

        state.genPTR = self.genPTR
        state.current_region = <void*> region



        state.E = E

        cdef int i
        cdef Particle p

        if particle == "photon":
            for i in range(self.N):
                p = Photon._new(state)
                p.throwAZIMUTH()
                self.container.append(p)
        elif particle == "electron":

            for i in range(self.N):
                p = Electron._new(state)
                p.throwAZIMUTH()
                self.container.append(p)
        else:
            for i in range(self.N):
                p = Positron._new(state)
                p.throwAZIMUTH()
                self.container.append(p)



cdef SquareSurface(Source):
    pass
