#internal imports

from .particles import electrons as el
from .particles import photons   as ph

#import .particles.electrons as el
#import .particles.photons   as ph

#import montecarlo.electrons as el
#import montecarlo.photons as ph


#external imports
from numba import *
import time

from pyquaternion import Quaternion
rotate = lambda vec,axis,theta: Quaternion(axis=axis, angle=theta).rotate(vec)


#hopefully I can remove this
#from mayavi.mlab import *

from numpy import *
from numpy.random import rand

from .tools.vectors import Vector
from .particles.particle import choose # overriding choose function of numpy
#from montecarlo.particle import choose 
from .tools.performance import timer



class Source:
    @timer
    def __init__(self,
                 space,
                 particle = "p",
                 pos      = (0,0,0),
                 N        = 100,
                 E        = 10,                  #to do: suport for lists
                 axis     = (0, 0),
                 theta    = 0,                   #to do: suport for lists
                 phi      = 0,
                 simulate_secondary = False):    #to do: suport for lists
        
        """
        This is a doc file asking  myself to write the doc file.....
        """
        
        if space == []:              raise ValueError("> Space can't be empty!")
        if type(space) is not list:  raise ValueError("> Space must be a list!")
        
        if particle == "p":
            self.ParticleType = ph.Photon
            print("Hi. I'm a photon source. ^.^")

        if particle == "e":
            print("Hi. I'm an electron source. ~.~")
            self.ParticleType = el.Electron
        
        self.space     = space
        self.simulated = False
       
        #position
        if callable(pos):
            _pos = pos

        if type(pos) is tuple:
            self.pos = Vector(pos[0], pos[1], pos[2])
            def _pos():
                return self.pos

        if type(pos) is dict:
            pos = list(pos.items())
            self.pos_args = [item for tup in pos for item in tup]
            def _pos():
                return choose(*self.pos_args)
        
        #PROB:E
        #spectrum
        if type(E) is dict:
            E  = list(E.items())
            self.args = [item for tup in E for item in tup]
            print(self.args) 
            def _E():
                return choose(*self.args)

        if type(E) is int or type(E) is float:
            def _E():
                return E

        if callable(E):
            _E = E

        if callable(phi):
            _phi = phi

        if type(phi) is int or type(phi) is float:
            def _phi():
                return phi 

        if callable(theta):
            _theta = theta

        if type(theta) is int or type(theta) is float:
            def _theta():
                return theta




        self.ex = Vector(1., 0., 0.)
        self.ey = Vector(0., 1., 0.)
        self.ez = Vector(0., 0., 1.)
        #self.change_direction(axis[0], axis[1])

        print("> Creating source...")
        #Populating Source
        t0 = time.perf_counter()
        #from tqdm import trange
        for _ in range(N):
            #progress = 100*_/N
            #if progress % 5 == 0:
                #tf = time.perf_counter()
                #print('Progress: ', progress, '%', ' - ', (tf-t0)/60, 'min')

            self.particle = self.ParticleType(pos   = _pos(),
                                              theta = _theta(),
                                              phi   = _phi(),
                                              E     = _E(),
                                              space = self.space,
                                              ex    = self.ex,
                                              ey    = self.ey,
                                              ez    = self.ez,
                                              simulate_secondary = simulate_secondary)
            #self.particle.run()
            #del self.particle

        print(f"""
                ------------------------------------------------
                Source Paramaters (these might be just samples)
                ------------------------------------------------
                pos   = {_pos()}
                E     = {_E()}
                theta = {_theta()}
                phi   = {_phi()}

                simulate_secondary = {simulate_secondary}
                -----------------------------------------------
                """)
        print(f"> Source has been populated with {N} particles!")
        print(f">> note: simulate_secondary = {simulate_secondary}")

 






    @timer
    def simulate(self):
        """
        Run simulation.
        """

        if self.simulated is True:
            raise RuntimeError("> Source has already been simulated.")

        print("> Simulating source...")
        t0 = time.perf_counter()
        N = len(self.particle_)

        
        #from tqdm import tqdm
        for p in self.particle_:
            p.run()

            
        print("> Done.")
        self.simulated = True
        
        #print("Deleting particles...")
        #for n, _ in enumerate(self.particle_):
        #    del self.particle_[n]
        #self.name = ""
        #print('Saving results to {self.name}.txt')


    def change_direction(self, theta, phi):
        '''Rotate the frame.'''
        axis = self.ez
        self.ey = rotate(self.ey, axis, phi)
        self.ex = rotate(self.ex, axis, phi)

        axis = self.ey
        self.ez = rotate(self.ez, axis, theta)
        self.ex = rotate(self.ex, axis, theta)

    def __setattr__(self, name, value):
        """
        Overloaded __setattr__: keep a record of all previous values of a given variable
        """
        if name not in self.__dict__:
            self.__dict__[name + "_"]  = [value]
            self.__dict__[name]        = value
        else:
            self.__dict__[name + "_"] += [value]
            self.__dict__[name]        = value

    ## Plot
    def plot(self, fig = None, plot_secondary = False):
        if fig is None:
            fig = figure()

        print(fig)
        for particle in self.particle_:
            particle.add_plot(fig, plot_secondary = plot_secondary)
        return fig

    def _plot(self):
        from vispy import app, scene

        canvas = scene.SceneCanvas(keys = 'interactive',  title='plot3d', show=True)
        view   = canvas.central_widget.add_view()
        view.camera = 'turntable'
        view.camera.fov = 45
        view.camera.distance = 6

        print("> Plotting particles...")
        for particle in self.particle_:
            particle._add_plot(view.scene)

        print("> Running app...")
        app.run()

