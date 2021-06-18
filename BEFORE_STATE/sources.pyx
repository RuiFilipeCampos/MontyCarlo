# cython: profile=True

print(">>>>>   IMPORTING SOURCES")


from .particles.particle cimport STATE
from particles.photons cimport Photon
from particles.electrons cimport Electron
from particles.particle cimport Particle
from tools.vectors cimport Vector
from geometry.main cimport Volume
from collections import deque


from random.mixmax.interface cimport mixmax_engine


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
    cdef object particles, photon_mesh, electron_mesh, positron_mesh
    
    
    
    def __str__(self):
        
        to_print = ">>>>>> SOURCE \n"
        
        for particle in self.particles:
            if isinstance(particle, Photon): to_print += ">>> PHOTON \n"
            elif isinstance(particle, Electron): to_print += ">>> Electron \n"
            elif isinstance(particle, Positron): to_print += ">>> Positron \n"
            
            X, Y, Z = particle.track()
            
            for x, y, z in zip(X, Y, Z):
                to_print += f"{x}         | {y}      |  {z}  \n"
                
        return to_print
            

    cpdef void set_mesh(self):
    

        
        electron_cell = deque()
        electron_points = deque()
        cdef int electron_Npoints = 0
        
        photon_cell = deque()
        photon_points = deque()
        cdef int photon_Npoints = 0
        
        positron_cell = deque()
        positron_points = deque()
        cdef int positron_Npoints = 0

        print("getting particles")
        for particle in self.particles:
            if isinstance(particle, Photon):
                photon_Npoints += particle.add_to_cell(photon_cell, photon_points, photon_Npoints)
                #
            
            elif isinstance(particle, Electron):
                electron_Npoints += particle.add_to_cell(electron_cell, electron_points, electron_Npoints)

                #
                
            elif isinstance(particle, Positron):
                positron_Npoints += particle.add_to_cell(positron_cell, positron_points, positron_Npoints)
                
        print("adding meshes")


        
        self.photon_mesh = pv.PolyData(np.array(photon_points),  lines = np.array(photon_cell), n_lines = photon_Npoints)
        self.electron_mesh = pv.PolyData(np.array(electron_points) , lines = np.array(electron_cell), n_lines = electron_Npoints)
        self.positron_mesh = pv.PolyData( np.array(positron_points),  lines = np.array(positron_cell), n_lines = positron_Npoints)
                                         

                                         
                                         
        
    def new_plot(self, ph_opacity = .1,el_opacity = 1, po_opacity = 1):    
        import pyvista as pv
        pv.set_plot_theme("night")
        plotter = pv.Plotter()
        
        
        
        plotter.add_mesh(self.photon_mesh, line_width = 0.001,    color='white', opacity=ph_opacity)
        plotter.add_mesh(self.electron_mesh, line_width = 0.001,  color='blue', opacity=el_opacity)
        plotter.add_mesh(self.positron_mesh,line_width = 0.001,   color='red', opacity=po_opacity)
                
        return plotter


cdef class DirectionalElectrons(Source):
    cdef public object photon_spec
    cdef int N
    cdef int Ntot
    cdef Volume region 
    cdef double E
    
    
    

        
        
        
    
    def __init__(self, int N, double E, Volume current_region, Volume space):
        self.photon_spec = deque()
        self.N = N
        self.E = E
        cdef Vector pos = Vector(0, 0, 0)
        cdef Vector ex, ey, ez
        self.region = current_region
        
        ex = Vector(1, 0, 0)
        ey = Vector(0, 1, 0)
        ez = Vector(0, 0, 1)
        cdef Electron p
        
        self.particles = deque()
        for _ in range(N):
            

            p = Electron._new(E, 
                              0, 0, 0,
                              0, 1, 0,
                              0, 0, 1,
                              current_region)
            
            p.deposit()
            self.particles.append(p)
        
    
    
    cpdef run(self, output = None,secondary = True,  mod = 100):
        from time import perf_counter_ns
        
        cdef double t0, tf, dt
        
        
        if not secondary:
            t0 = perf_counter_ns()
            self._runNS()
            tf = perf_counter_ns()
        elif output is None:
            t0 = perf_counter_ns()
            self._run()
            tf = perf_counter_ns()
        else:
            t0 = perf_counter_ns()
            self._runWITHPROGRESS(mod)
            tf = perf_counter_ns()
            
        dt = tf - t0
        print("")
        print(f"INITIAL ENERGY: {self.E}")
        print(f"ELAPSED: {dt/60*1e-9}min")
        print(f"PER PARTICLE: {1e-3*dt/self.N}microseconds")
        print(f"PER TRACK: {1e-3*dt/self.Ntot}microseconds")
        print(f"# OF SECONDARY PARTICLES PRODUCED: {self.Ntot - self.N}")
        print(self.region.material.__repr__())
        return self.photon_spec
    
    
    cpdef spectrum(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        while nCURRENT < nSIMULATE:
            #print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            if isinstance(p, Photon):
                self.photon_spec.append(p.E)
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
            #print(p)
        print(len(self.particles))
        self.Ntot =  nSIMULATE
        return self.photon_spec

    
    
    
    cdef void _runWITHPROGRESS(self, int mod):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        cdef object simulated = deque()
        while nCURRENT < nSIMULATE:
            
            print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            # if isinstance(p, Photon):
            #     self.photon_spec.append(p.E)
            
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            simulated.append(p)
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
        self.particles = simulated
        self.Ntot =  nSIMULATE

            #print(p)
            #if nCURRENT % mod == 0:
                #print(f"{nCURRENT}/{nSIMULATE}")
            
    cdef void _run(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        while nCURRENT < nSIMULATE:
            #print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
            #print(p)
        self.Ntot =  nSIMULATE
    cdef void _runNS(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        while nCURRENT < nSIMULATE:
            #print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            #self.particles.extend(p.secondary)
            
            #nSIMULATE += p.nSECONDARY
            nCURRENT += 1
            #print(p)
        self.Ntot =  nSIMULATE
               
            
    def get_particles(self):
        return self.particles



from .particles.positrons cimport Positron
cdef class DirectionalPositrons(Source):
    cdef public object photon_spec
    cdef int N
    cdef int Ntot
    cdef Volume region 
    cdef double E
    
    def __init__(self, int N, double E, Volume current_region, Volume space):
        self.photon_spec = deque()
        self.N = N
        self.E = E
        cdef Vector pos = Vector(0, 0, 0)
        cdef Vector ex, ey, ez
        self.region = current_region
        
        ex = Vector(1, 0, 0)
        ey = Vector(0, 1, 0)
        ez = Vector(0, 0, 1)
        cdef Positron p
        
        self.particles = deque()
        for _ in range(N):
            

            p = Positron._new(E, 
                              0, 0, 0,
                              0, 1, 0,
                              0, 0, 1,
                              current_region)
            
            p.deposit()
            self.particles.append(p)
        
    
    
    cpdef run(self, output = None,secondary = True,  mod = 100):
        from time import perf_counter_ns
        
        cdef double t0, tf, dt
        
        
        if not secondary:
            t0 = perf_counter_ns()
            self._runNS()
            tf = perf_counter_ns()
        elif output is None:
            t0 = perf_counter_ns()
            self._run()
            tf = perf_counter_ns()
        else:
            t0 = perf_counter_ns()
            self._runWITHPROGRESS(mod)
            tf = perf_counter_ns()
            
        dt = tf - t0
        print("")
        print(f"INITIAL ENERGY: {self.E}")
        print(f"ELAPSED: {dt/60*1e-9}min")
        print(f"PER PARTICLE: {1e-3*dt/self.N}microseconds")
        print(f"PER TRACK: {1e-3*dt/self.Ntot}microseconds")
        print(f"# OF SECONDARY PARTICLES PRODUCED: {self.Ntot - self.N}")
        print(self.region.material.__repr__())
        return self.photon_spec
    
    
    cpdef spectrum(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        while nCURRENT < nSIMULATE:
            #print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            if isinstance(p, Photon):
                self.photon_spec.append(p.E)
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
            #print(p)
        print(len(self.particles))
        self.Ntot =  nSIMULATE
        return self.photon_spec

    
    
    
    cdef void _runWITHPROGRESS(self, int mod):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        cdef object simulated = deque()
        while nCURRENT < nSIMULATE:
            
            print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            # if isinstance(p, Photon):
            #     self.photon_spec.append(p.E)
            
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            simulated.append(p)
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
        self.particles = simulated
        self.Ntot =  nSIMULATE

            #print(p)
            #if nCURRENT % mod == 0:
                #print(f"{nCURRENT}/{nSIMULATE}")
            
    cdef void _run(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        while nCURRENT < nSIMULATE:
            #print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
            #print(p)
        self.Ntot =  nSIMULATE
    cdef void _runNS(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        while nCURRENT < nSIMULATE:
            #print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            #self.particles.extend(p.secondary)
            
            #nSIMULATE += p.nSECONDARY
            nCURRENT += 1
            #print(p)
        self.Ntot =  nSIMULATE
               
            
    def get_particles(self):
        return self.particles




cdef class IsotropicSourceEl(Source):
    cdef int N, Ntot
    cdef double E
    cdef Volume region
    
    def __init__(self, int N, double E, Volume current_region, Volume space):
        self.N = N
        self.E = E
        self.region = current_region
        cdef Vector pos = Vector(0, 0, 0)
        cdef Vector ey, ez
        
        ey = Vector(0, 1, 0)
        ez = Vector(0, 0, 1)
        
        
        self.particles = deque()
        cdef int _
        for _ in range(N):
            
            self.particles.append(Electron._newISOTROPIC(E, 
                              0, 0, 0,
                              current_region, genPTR))
            
    def run(self):
        self._runWITHPROGRESS(500)
    def get_particles(self):
        return self.particles 
    
    cdef void _runWITHPROGRESS(self, int mod):
        import time
        
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        
        simulated = deque()
        while nCURRENT < nSIMULATE:
            
            print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            simulated.append(p)
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
        self.particles = simulated
        self.Ntot = nSIMULATE



cdef class IsotropicSourcePh(Source):
    cdef int N, Ntot
    cdef double E
    cdef Volume region
    
    def __init__(self, int N, double E, Volume current_region, Volume space):
        self.N = N
        self.E = E
        self.region = current_region
        cdef Vector pos = Vector(0, 0, 0)
        cdef Vector ey, ez
        
        ey = Vector(0, 1, 0)
        ez = Vector(0, 0, 1)
        
        
        self.particles = deque()
        cdef int _
        for _ in range(N):
            
            self.particles.append(Photon._newISOTROPIC(E, 
                              0, 0, 0,
                              current_region, genPTR))
            
    def run(self):
        self._runWITHPROGRESS(500)
    def get_particles(self):
        return self.particles 
    
    cdef void _runWITHPROGRESS(self, int mod):
        import time
        
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        
        simulated = deque()
        while nCURRENT < nSIMULATE:
            
            print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            simulated.append(p)
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
        self.particles = simulated
        self.Ntot = nSIMULATE







cdef class DirectionalPhotons(Source):
    cdef int N, Ntot
    cdef double E
    cdef Volume region
    
    def __init__(self, int N, double E, Volume current_region):
        self.N = N
        self.E = E
        self.region = current_region
        cdef Vector pos = Vector(0, 0, 0)
        cdef Vector ey, ez
        
        ey = Vector(0, 1, 0)
        ez = Vector(0, 0, 1)
        
        cdef Photon p
        self.particles = deque()
        cdef int _
        for _ in range(N):
            
            p = Photon._new(E, 
                              0, 0, 0,
                              0, 1, 0,
                              0, 0, 1,
                              current_region)
            p.deposit()
            self.particles.append(p)

            


        
    
    
    cpdef run(self, output = None, mod = 100, secondary = True):
        from time import perf_counter_ns
        cdef int i
        cdef double t0, tf
        if not secondary:
            t0 = perf_counter_ns()
            self._runNS()
            tf = perf_counter_ns()
        elif output is None:
            t0 = perf_counter_ns()
            i = self._run()
            tf = perf_counter_ns()

        else:
            t0 = perf_counter_ns()
            self._runWITHPROGRESS(mod)
            tf = perf_counter_ns()
            
        cdef double dt = tf - t0
        print("")
        print(f"INITIAL ENERGY: {self.E}")
        print(f"ELAPSED: {dt/60*1e-9}min")
        print(f"PER PARTICLE: {1e-3*dt/self.N}microseconds")
        print(f"PER TRACK: {1e-3*dt/self.Ntot}microseconds")
        print(f"# OF SECONDARY PARTICLES PRODUCED: {self.Ntot - self.N}")
        print(self.region.material.__repr__())
    
    
    cdef void _runWITHPROGRESS(self, int mod):
        import time
        
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        
        simulated = deque()
        while nCURRENT < nSIMULATE:
            
            print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            simulated.append(p)
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
        self.particles = simulated
        self.Ntot = nSIMULATE

            
    cdef int _run(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        while nCURRENT < nSIMULATE:
            #print(f"{nCURRENT}/{nSIMULATE}")
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
            #print(p)
        self.Ntot = nSIMULATE
        return nSIMULATE

    def run1(self):
        self._run1()
        
    cdef void _run1(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        
        
        
        from time import perf_counter
        cdef double t0 = perf_counter()
        cdef int i
        while nCURRENT < nSIMULATE:
            print(nCURRENT, nSIMULATE)
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1

            
        self.Ntot = nSIMULATE 
        print(nSIMULATE)           
  
    cdef void _run0(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        while nCURRENT < self.N:
            #print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            self.particles.extend(p.secondary)
            
            nSIMULATE += p.nSECONDARY
            nCURRENT += 1
            #print(p)
        self.Ntot = nSIMULATE    


    cdef void _runNS(self):
        cdef Particle p
        cdef int nSIMULATE = self.N
        cdef int nCURRENT = 0
        
        while nCURRENT < nSIMULATE:
            #print(f"{nCURRENT}/{nSIMULATE}")
            
            p = self.particles.popleft()
            
            #print(p)
            p._run(genPTR)
            #self.particles.extend(p.secondary)
            
            #nSIMULATE += p.nSECONDARY
            nCURRENT += 1
            #print(p)
        self.Ntot =  nSIMULATE
           
            
    def get_particles(self):
        return self.particles 
    
    
    
    # cdef void _run(self):
    #     cdef Particle p
    #     cdef int nSIMULATE = self.N
    #     cdef int nCURRENT = 0
        
        
    #     secondary = []
        
        
    #     while nCURRENT < self.N:
    #         #print(f"{nCURRENT}/{nSIMULATE}")
            
    #         p = self.particles[nCURRENT]
            
    #         #print(p)
    #         p._run(genPTR)
    #         self.particles += p.secondary
            
    #         nSIMULATE += p.nSECONDARY
            
    #         self.particles[nCURRENT] = 0
    #         del p
            
    #         nCURRENT += 1
    #         #print(p)
    #     self.Ntot = nSIMULATE
        





# 
# import mayavi.mlab as mlab
# import numpy as np

# from geometry.primitives cimport Volume, UnboundedVolume

# cdef class DirectionalPointSource:
#     cdef double[:] xx, yy, zz, ss
#     cdef double[:, :] CONNECT
    
    
#     @staticmethod
#     def new(SPACE, CURRENT, N):
#         return DirectionalPointSource._new(SPACE, CURRENT, N)
    
#     @staticmethod
#     cdef DirectionalPointSource _new(UnboundedVolume SPACE, 
#                                      Volume CURRENT,
#                                      int N):
        
#         cdef DirectionalPointSource self
#         self = <DirectionalPointSource>DirectionalPointSource.__new__(DirectionalPointSource)
        
#         cdef int index = 0
#         cdef list x, y, z, connections, s
#         x = list()
#         y = list()
#         z = list()
#         connections = list()
#         s = list()
#         cdef list X, Y, Z
        
#         cdef int i
#         for i in range(N):
        
#             p = Photon.new(SPACE, CURRENT)
#             p.run()
#             X, Y, Z = p.getTrack()
#             N = len(X)
#             S = p.getEnergy()
#             S.append(10)
            
#             x.append(np.array(X))
#             y.append(np.array(Y))
#             z.append(np.array(Z))
#             s.append(np.array(S))
            
#             connections.append(np.vstack(
#                             [np.arange(index,   index + N - 1.5),
#                             np.arange(index + 1, index + N - .5)]
#                                 ).T)
        
        
#             index += N
    
#         print("FINISHED SIMULATION")
    
#         # Now collapse all positions, scalars and connections in big arrays
        
#         self.xx = np.hstack(x)
#         self.yy = np.hstack(y)
#         self.zz = np.hstack(z)
#         self.ss = np.hstack(s)
#         self.CONNECT = np.vstack(connections)
#         return self
    
#         # Create the points

    
    
#     def plot(self):
#         src = mlab.pipeline.scalar_scatter(self.xx, self.yy, self.zz)
#         src.mlab_source.dataset.lines = np.array(self.CONNECT)
#         src.update()
        
#         mlab.pipeline.surface(src, colormap='jet', line_width=1, opacity=.04)
#         mlab.show()
        


            
            
            
        