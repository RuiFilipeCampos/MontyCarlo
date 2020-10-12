
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
#from numba import *

#from mayavi.mlab import *

#import pickle #probly not needed anymore
#import copy   #probly not needed anymore
import bisect

#searchsorted = njit(lambda a, x: searchsorted(a, x))
#insort_left = njit(bisect.insort_left)



###### SIMULATION PARAMETERS

eps   = 0.1  #size of the smallest pathlenght
error = 0.001 #fraction of photons that would miss eps






###################DEBUG##############################################3
#should be its own module
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

@njit
def dist(P, Q):
    """Find distance between two points."""
    return sqrt(sum((P-Q)**2))
    #return sqrt(sum((p-q)**2 for p, q in zip(P, Q)))



class StopSimulation(Exception):
    """
    Custom Type Error to stop the simulation.
    """
    pass

class Particle:
    """
    Contains logic for all particle types.
    """
                   
    
    def __init__(self,
                 space  = [],
                 E      = 6,
                 pos    = [0, 0, 0],
                 theta  = 0,
                 phi    = 0,
                 ex = [1,0,0], ey = [0,1,0], ez = [0,0,1],
                 simulate_secondary = False):
        
        self.simulate_secondary = simulate_secondary
        
        self.space = space
        self.pos   = array(pos)

        self.ex, self.ey, self.ez = map(array, [ex, ey, ez])

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

        self.total_distance = 0

        self.children = None


        
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
    def propagate(self):
        #this thing is gonna have to be an interface to C or Fortran
        #I can't see any other way
        #for now, I believe I can refactor all this so that I can use JIT
        """
        Propagates particle.
        - sample from exp(-mu*L)
        - randomly search for any boundary crossing - search is complete once boundary is found, or the probability
        of a region of size eps existing within the trajectory is smaller than 'error'
        - construct proposal
        - return proposal if in region
        - if not, commence crossing algorithm (it assumes that there is no extra boundaries)
        """

        
        self.pos = Particle._propagate(self.mu_tot,
                                       self.pos,
                                       self.ez,
                                       self.current_region.condition)

        if self.pos not in self.current_region:
            self.current_region = self.whereAmI()

            if not self.current_region.imp == 1:
                raise StopSimulation
        #bint self.current_region.imp == 1
        
            self.update_coefs()
            self.propagate()
            
        #return None

     
    #@njit
    def _propagate(mu_tot, pos, ez, cond):


        
        r = rand()
        L = -log(r)/mu_tot
        
        dL = 0.

        prob = 1 - eps/L
        
        samples = [0]
                        
        while not prob < error:
            dL = rand()*L
            if cond(pos + ez*dL):
                i = 0
                for element in samples:
                    if samples[i] > L:
                        samples.insert(i, L)
                        break
                    else: i += 1

                prob *= prob
            else:
                L = dL
                i = 0
                for element in samples:
                    if samples[i] > L:
                        samples = samples[0:i]
                        prob = (1 - eps/L)**(len(samples))
                        break
                    else:
                        i += 1
                    
                
                
                

        proposal = pos + ez * L

        if cond(proposal):
            return proposal
        
##
##        while dL < L:
##            if not cond(pos + dL*ez):
##                proposal = pos + dL*ez
##                break
##            dL += 0.1
##        else: return pos + ez*L


        
        initial = pos + ez*samples[-1]
        final   = proposal

        distance = dist(initial, final)

        while distance > 0.00001:
            mid_point = initial + ez * distance *0.5
            if cond(mid_point): initial = mid_point
            else:               final   = mid_point
            
            distance = dist(initial, final)

        return final


        



        
        
                
                
            
                
        
##        if proposed_pos not in self.current_region:
##            self.walk_to_surface(self.pos, proposed_pos)
##            self.current_region = self.whereAmI()
##
##            state = self.summary()
##            try:    self.current_region.tally += [state]
##            except: self.current_region.tally  = [state]
##                    
##            assert self.current_region.imp == 1
##            
##            self.update_coefs()
##            self.propagate()
##        else:
##            answer, pos = self.any(self.pos, proposed_pos)
##            if answer:
##                self.walk_to_surface(self.pos, pos)
##                self.update_coefs()
##                self.propagate()
##            else:
##                self.pos = proposed_pos
##
##
##
##            
##           
##
##
##
##
##
##    def any(self, initial_pos, final_pos):
##        dL = sqrt(sum( (final_pos - initial_pos)**2 ))
##        max_size = 2
##        min_size = max_size/2
##        
##        if dL < min_size:
##            return (False, final_pos)
##        
##        size = rand()*(max_size - min_size) + max_size/2
##        prob = 1 - size/dL
##        
##        while prob > 0.01:
##            probe_pos = self.pos + self.ez*rand()*dL
##            
##            if probe_pos not in self.current_region:
##                return (True, probe_pos)
##            
##            prob *= prob
##            
##        return (False, final_pos)
##
##
##
##
##    def walk_to_surface(self, initial_pos, final_pos):
##        distance = dist(initial_pos, final_pos)
##
##        while distance > 0.01:
##            
##            new_pos = initial_pos + 0.5*distance*self.ez
##            if new_pos in self.current_region:
##                answer, pos = self.any(initial_pos, new_pos)
##                if answer: final_pos   = pos
##                else:      initial_pos = new_pos
##            else:
##                final_pos = new_pos
##
##            distance = dist(initial_pos, final_pos)
##
##        self.pos = final_pos

        
        

    #maybe use this to store theta and phi
##    def change_direction(self, theta, phi):
##        '''Rotates the local frame.'''
##        axis = self.ez
##        self.ey = rotate(self.ey, axis, phi)
##        #self.ex = rotate(self.ex, axis, phi)
##
##        axis = self.ey
##        self.ez = rotate(self.ez, axis, theta)
##        #self.ex = rotate(self.ex, axis, theta)


    
    def change_direction(self, cos, phi):
        '''Rotates the local frame.'''
        axis = self.ez
        self.ey = rotateAngle(self.ey, axis, phi)
        #self.ex = rotate(self.ex, axis, phi)

        axis = self.ey
        self.ez = rotateCos(self.ez, axis, cos)
        #self.ex = rotate(self.ex, axis, theta)















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



