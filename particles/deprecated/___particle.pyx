
# cython: profile=True
"""
particle.py
----------------------------------------------------------------------------
Code that is inherited by all particle types.
Handles reference frame stuff, propagation of particle, boundary crossing.

It also contains some functions that turn out to be useful just about anywhere:
- choose -> chooses elements with provided probabilities
- dist   -> calculates distance between two points
----------------------------------------------------------------------------


There is some work to be done here.

Problems:
. I am not sure if the self.__setattri__ method should be overloaded here.
I do like to have it, it makes the code cleaner in many situations. However, I'm pretty
sure it is slowing down the simulation. I am considering deprecating it.
Also considering putting it as a debug option.

. The use of pyquaternion - I believe that I'm only using one equation from quaternion
theory, probably don't need to import an entire module for that..

. the simulation parameters shouldn't be a variable in this namespace. Need to find a way
to make them user definable. (maybe move them to the sources module)

"""


#External Imports
from numpy import *



from numpy.random import *
from pyquaternion import Quaternion
import bisect

#Internal Imports
from ..tools.performance import timer

 

###### SIMULATION PARAMETERS

cdef double eps   = 10  #size of the smallest pathlenght
cdef double error = 0.001 #fraction of photons that would miss eps





def rotateAngle(vec, axis, theta):
    """ Rotates vec by an angle theta along the provided axis """
    return Quaternion(axis = axis, angle = theta).rotate(vec)

#rotate = lambda vec, axis, theta: Quaternion(axis=axis, angle=theta).rotate(vec)


def rotateCos(vector, axis, c):
    """Rotate vector along axis with angle defined by c = cos(theta)"""
    vector = Quaternion(0, *vector)
    
    c2 = sqrt(.5*(1 + c))
    s2 = sqrt(.5*(1 - c))
    
    q = Quaternion(c2, *(s2*axis))
    return (q * vector * q.conjugate).imaginary


def choose(*args):
    """
    Return elements with provided probabilities.
    """
    n = len(args)

    if n % 2 != 0:
        raise Exception("Nummber of elements not equal to number of coefs.")

    elements = [args[2*i]   for i in range(n//2)]
    coefs    = [args[2*i+1] for i in range(n//2)]

    
    #elements = [args[i] for i in range(0, n//2)]
    #coefs    = [args[i] for i in range(n//2, n)]

    coefs = array(coefs)
    probs = coefs/sum(coefs)
    cumul = [sum(probs[0:i]) for i in range(n//2)] + [1]
    #print(probs)
    r = rand()
    
    i = searchsorted(cumul, r, side="left")
    return elements[i - 1]    


def dist(P, Q):
    """Find distance between two points."""
    return sqrt(sum((p-q)**2 for p, q in zip(P, Q)))


class StopSimulation(Exception):
    """ Custom Type Error to stop the simulation. """
    pass

cimport numpy as np


cdef class CParticle:
    cdef public double E, theta, phi
    cdef public bint simulate_secondary
    #cdef np.ndarray[np.float64_t, ndim=3] pos, ex, ey, ez
    cdef object pos, ex, ey, ez
    cdef public list space

    def __init__(self, 
                 list space, 
                 double E, 
                 object pos, 
                 double theta, 
                 double phi,
                 object ex,
                 object ey,
                 object ez,
                 bint simulate_secondary):


        self.simulate_secondary = simulate_secondary
        
        self.space = space
        self.pos   = pos

        self.ex, self.ey, self.ez = ex, ey, ez

        self.theta = theta
        self.phi   = phi
        self.change_direction(cos(self.theta), self.phi)
        
            
        self.E = E
        self.current_region = self.whereAmI()

        if self.__class__.__name__ == 'Photon':
            self.current_region.E0 += self.E
            self.current_region.N0 += 1
            self.current_region.E += self.E
            self.current_region.Nin += 1
            
        self.update_coefs()
        #self.children = None




        
##    def __setattr__(self, name, value):
##        """Overloaded __setattr__: keep a record of all previous values of a given variable"""
##        if name not in self.__dict__:
##            self.__dict__[name + "_"] =  [value]
##            self.__dict__[name] = value
##        else:
##            self.__dict__[name + "_"] += [value]
##            self.__dict__[name] = value
            
    

            
    def whereAmI(self):
        ''' Returns the region object that contains the position of the particle'''
        
        for region in self.space:
            if self.pos in region:
                return region

    #@timer
            #testar contra força bruta
    cdef bint propagate(self):
        
        cdef double r    = rand()
        cdef double L    = -log(r)/self.mu_tot
        cdef double prob = 1 - eps/L
            
        cdef list samples = [0]

        cdef int k
        cdef double dL

        while not prob < error:
            dL = rand()*L
            if self.pos + self.ez*dL in self.current_region:
                bisect.insort_left(samples, dL)
                prob *= prob
            else:
                L = dL
                k = searchsorted(samples, L, side = 'left')
                samples = samples[0:k]

                prob = (1 - eps/L)**(len(samples))

        proposal = self.pos + self.ez * L
        
        if proposal in self.current_region:
            self.pos = proposal
            return True

        
        #cdef np.ndarray[np.float_t, ndim=3] initial = self.pos + self.ez*samples[-1]
        #cdef np.ndarray[np.float_t, ndim=3] final   = proposal
        #cdef np.ndarray[np.float_t, ndim=3] mid_point

        cdef object initial = self.pos + self.ez*samples[-1]
        cdef object final   = proposal
        cdef object mid_point

        cdef float distance = dist(initial, final)

        while distance > 0.00001:
            mid_point = initial + self.ez * distance / 2
            if mid_point in self.current_region: initial = mid_point
            else:                                final   = mid_point
            distance = dist(initial, final)
                


        
        self.pos = final

        self.current_region.E -= self.E
        self.current_region.Nout += 1
        
        self.current_region = self.whereAmI()

        self.current_region.Nin += 1
        self.current_region.E += self.E

        if not self.current_region.imp == 1:
            raise StopSimulation
        #bint self.current_region.imp == 1
        
        self.update_coefs()
        self.propagate()
        return False
        
        
                
            

        
        




    
    def change_direction(self, float cos, float phi):
        '''Rotates the local frame.'''
        cdef object axis
        axis = self.ez
        self.ey = rotateAngle(self.ey, axis, phi)
        #self.ex = rotate(self.ex, axis, phi)

        axis = self.ey
        self.ez = rotateCos(self.ez, axis, cos)
        #self.ex = rotate(self.ex, axis, theta)



















    #maybe use this to store theta and phi
##    def change_direction(self, theta, phi):
##        '''Rotates the local frame.'''
##        axis = self.ez
##        self.ey = rotate(self.ey, axis, phi)
##        #self.ex = rotate(self.ex, axis, phi)
##
##        axis = self.ey
##        self.ez = rotate(self.ez, axis, theta)
##        #self.ex = rotate(self.ex, axis, theta



class Particle(CParticle):
    """  Contains logic for all particle types. """
                
    
    def __init__(self,
                 space  = [],
                 E      = 6.,
                 pos    = [0., 0., 0.],
                 theta  = 0.,
                 phi    = 0.,
                 ex = [1.,0.,0.], ey = [0.,1.,0.], ez = [0.,0.,1.],
                 simulate_secondary = False):
        
        ex, ey, ez = map(array, [ex, ey, ez])
        pos = array(pos)

        super().__init__(space, E, pos, 
                         theta, phi, 
                         ex, ey, ez, 
                         simulate_secondary)



    def summary(self):
        return {'pos': self.pos,
                'E'  : self.E,
                'ez' : self.ez}


    
    def plot(self, tube_radius = 1, line_width = 1):
        """Plot the path of the particle. (uses mayavi)"""
        history = array(self.pos_)
        X, Y, Z = history[:,0], history[:,1], history[:,2]
        plot3d(X, Y, Z, tube_radius = tube_radius, line_width = line_width)
        show()



    def _add_plot(self, fig):
        """
        Don't use this. Just an attempt of using other modules for visualizing.
        Keeping it here to remember later that there are other options.
        (there is also paraview, but it's very time consuming to learn that)s
        """
        from vispy import scene, visuals
        Plot3D = scene.visuals.create_visual_node(visuals.LinePlotVisual)


        history = array(self.pos_)
        X,Y,Z = history[:,0], history[:,1], history[:,2]
        pos = c_[X, Y, Z]
        
        Plot3D(pos, width = 2.0, color='red', edge_color='w', symbol='o',
              face_color = (0.2, 0.2, 1, 0.8),
              parent     = fig)
        
    def add_plot(self, fig, plot_secondary = False):
        """
        Add trajectory to a mayavi figure. Used by a source instance to plot the
        results of the simulation.
        """
        
        name = self.__class__.__name__
        if name == 'Photon': color = (0, 0, 1)
        else:                color = (1, 0, 0)
        
        history = array(self.pos_)
        X,Y,Z   = history[:,0], history[:,1], history[:,2]
        plot3d(X, Y, Z, figure=fig, tube_radius=.8, color=color)
        
        if plot_secondary is True:
            for p in self.children_[1:]:
                p.add_plot(fig, plot_secondary = True)