"""
Code that is specific to a photon type.

To do's here:
- would be cool to have a cut off for photons
- need to update the models
"""


#External Imports
from numpy import *
from numpy.random import rand, randint
#from numba import *
#import pickle

#Local Imports
#import montecarlo.particle as pa
#import montecarlo.electrons as e
#from montecarlo.particle import choose, rotate

from . import particle as pa
from .particle import choose, rotateAngle, rotateCos
from . import electrons as e


class Photon(pa.Particle):
    def __init__(self,
                 space = [],
                 E     = 6,
                 pos   = [0, 0, 0],
                 theta = 0,
                 phi   = 0,
                 ex = [1,0,0], ey = [0,1,0], ez = [0,0,1],
                 simulate_secondary = False):
        
        super().__init__(space   = space,
                         E       = E,
                         pos     = pos,
                         theta   = theta,
                         phi     = phi,
                         ex = ex, ey = ey, ez = ez,
                         simulate_secondary = simulate_secondary)

        self.k = E/0.511 

        #self.interactions = 
        



    #@timer
    def run(self):
        try:
            i = 0
            while 1:
                self.propagate()
                self.interaction = choose(self._coherent,       self.coherent.CS(self.E))
                                          #self._incoherent,     self.incoherent(self.E))
                                          #self._photoelectric,  self.photoelectric(self.E),
                                          #self._pairproduction, self.pairproduction(self.E))
                self.interaction()
                i += 1
                
        except pa.StopSimulation:
            pass
            #print("Number of interactions:", i)
        
        if self.simulate_secondary is True:
            if self.children is not None:
                for particle in self.children_[1:]:
                    particle.simulate_secondary = True
                    particle.run()

    def update_mu(self):
        """ I think this method is taking longer than it should!"""
        self.mu_tot = self.coherent.CS(self.E) #+ #self.incoherent(self.E) #+ self.photoelectric(self.E) + self.pairproduction(self.E)
        

    def update_coefs(self):
        self.current_material = self.current_region.material.photon
        
        self.coherent         = self.current_material.coherent
        #self.compton          = self.current_material.incoherent
        #self.photoelectric    = self.current_material.photoelectric
        #self.pairproduction   = self.current_material.pairproduction

        self.update_mu()
        
        
    
    def _coherent(self):
        """
        Samples angular deflection for the coherent scattering.
        To do: Low energy stuff.
        """
        
        #k = self.E/0.511 #probably should change k to be a state variable


        
        x_max = 20.6074*2*self.k
        r_max = self.coherent.FF.cumul(x_max) #internals of this needs work
        
        while 1:
            r  = rand()*r_max
            x2 = self.coherent.FF.invCum(r)
            cos  = 1 - 0.5 * x2 / (20.6074*self.k)**2
            if rand() < self.coherent.g(cos):
                self.change_direction(cos, 2*pi*rand())
                break
            



    def _compton(self):
        """
        There might be better ways to do this sampling.
        """

        k = self.E/0.511
        t1 = log(1+2*k)
        t2 = 2*k*(1+k)/(1+2*k)**2
        tau_min = 1/(1 + 2*k)
        
        for _ in range(10_000):
            r = 2*rand()
            i = 1 if r < t1/(t1 + t2) else 2
            r = rand()
            tau = tau_min**r if i == 1 else (tau_min**2 + r*(1 - tau_min**2))**.5

            N = (1 - tau)*((2*k-1)*tau - 1)
            D = k**2 * tau * (1 + tau**2)

            
            T = (1 - N/D)*self.coherent.S(x)
        


                #num = cos(theta/2)
                #den = sin(theta/2)*(1+self.E/0.511)
                #theta_el = arctan(num/den)
                
##                self.children = e.Electron(pos   = self.pos,
##                                           theta = theta_el,#                 confirm
##                                           phi   = phi + pi,#                 confirm
##                                           E     = self.E*(1-fraction),
##                                           space = self.space, 
##                                           ex    = self.ex,
##                                           ey    = self.ey,
##                                           ez    = self.ez,
##                                           simulate_secondary = self.simulate_secondary)
                
                #self.E = self.E*fraction
                #self.update_mu()
                #break
        else:
            raise RuntimeError("Rejection sampling exceeded the maximum number of iterations.")

    def _pairproduction(self):
        #assert False
        raise pa.StopSimulation

    def _photoelectric(self):
        if self.E < 0.511: theta = pi
        else:              theta = 0
        
##        self.children = e.Electron(pos   = self.pos,
##                                   theta = theta,#       to improve  
##                                   phi   = 2*pi*rand(),
##                                   E     = self.E,
##                                   space = self.space,
##                                   ex    = self.ex,
##                                   ey    = self.ey,
##                                   ez    = self.ez,
##                                   simulate_secondary = self.simulate_secondary)
        #assert False
        raise pa.StopSimulation






        
