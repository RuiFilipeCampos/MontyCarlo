


#External Imports
from numpy import *
from numpy.random import *
#from numpy.linalg import inv
#from matplotlib.pyplot import *
from pyquaternion import Quaternion
from scipy.integrate import quad


#Internal Imports
from . import particle as pa
from . import photons  as ph
from .particle import choose, dist, rotateAngle, rotateCos





########################################################################3
from functools import wraps
import time

def timer(func):
    @wraps(func) #copies metadata, id, docstring etc
    def wrapper(*args, **kwargs):
        
        t0 = time.time()
        result = func(*args, **kwargs)
        tf = time.time()
        
        print(func.__name__, tf-t0)
        
        return result
    return wrapper
########################################################################3










class Electron(pa.Particle):
    """Electron object."""

    def __init__(self,
                 space = [],
                 E     = 6,
                 pos   = [0,0,0],
                 theta = 0,
                 phi   = 0,
                 ex = [1,0,0], ey = [0,1,0], ez = [0,0,1],
                 simulate_secondary = False):
        
        super().__init__(space   = space,
                         E       = E*1e6,
                         pos     = pos,
                         theta   = theta,
                         phi     = phi,
                         ex = ex, ey = ey, ez = ez,
                         simulate_secondary = simulate_secondary)
        

    #@timer        
    def run(self):
        try:
            assert self.current_region.imp is 1
            while 1:
                self.propagate()
                self.interaction = choose(self._elastic,  self.elastic(self.E), 
                                          self._brem,     self.brem(self.E)  )
                self.interaction()
                self._inelastic()

                assert self.E > 50e3
                
        except AssertionError:
            pass
        
        if self.simulate_secondary is True:
            if self.children is not None:
                for particle in self.children_[1:]:
                    particle.simulate_secondary = True
                    particle.run()


    def update_coefs(self):
        self.current_material = self.current_region.material.electron
        self.brem             = self.current_material.brem
        self.elastic          = self.current_material.elastic
        self.sp               = self.current_material.sp

        self.update_mu()

        
    def update_mu(self):
        self.mu_tot = self.elastic(self.E) + self.brem(self.E)

    def _doHinge(self):
        
        self.elastic.update_dist(self.E)
        mu_c = self.elastic.mu_c()
        
        if mu_c == 0: return 0
        
        self.mu_c = mu_c
        last_pos  = self.__dict__["pos_"].pop()
        
        t = rand()*self.dL
        self.pos = self.pos_[-1] + self.ez*t
        
        if self.pos not in self.current_region:
            self.current_region = self.whereAmI()
            assert self.current_region.imp == 1
            self.update_coefs()

        try:
            a, b = self.elastic.F(mu_c, self.dL)
            mu   = b*rand() if rand() < a else b + (1-b)*rand()
        except AssertionError:
            mu = rand()
            
        self.mu    = mu
        self.theta = arccos(1-2*self.mu)
        self.phi   = 2*pi*rand()
        
        self.change_direction(self.theta, self.phi)
        self.pos = self.pos + self.ez*(self.dL - t)
        return mu_c


        

    def _elastic(self):
        #prob = self.elastic.prob(self.E)
        
        mu_c   = self._doHinge()
        invCum = self.elastic.invCum

        r_c = mu_c + (1-mu_c)*rand() #get angles above cut off
        mu  = invCum(r_c)
        
        self.theta = arccos(1-2*mu)
        self.phi   = 2*pi*rand()
        
        self.change_direction(self.theta, self.phi)
        
    def _inelastic(self):
        self.E -= self.dL*self.sp(self.E)
        self.update_mu()



    def _brem(self):
        self.brem.X.update_dist(self.E)
        self.k      = self.brem.X.invCum(rand())
        
        self.W = self.E*self.k
        
        self.children = ph.Photon(pos    = self.pos,
                                  theta  = 0,              # to do
                                  phi    = 2*pi*rand(),
                                  E      = self.W/1e6,
                                  space  = self.space,
                                  ex     = self.ex,
                                  ey     = self.ey,
                                  ez     = self.ez,
                                  simulate_secondary = self.simulate_secondary)
        self.E -= self.W
        self.update_mu()


