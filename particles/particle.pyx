# cython: profile=True

#External Imports
from numpy import *
from numpy.random import *
import numpy as np
#from pyquaternion import Quaternion
#import bisect

#Internal Imports
from ..tools.performance import timer
from ..tools.vectors cimport Vector
from ..geometry.primitives cimport Volume

cdef object choose(list cumul, list items):

	cdef double r = rand()*cumul[-1]
	cdef int i = 0


	for i in range(len(cumul)):
		if cumul[i] > r:
			return items[i-1]



class StopSimulation(Exception):
	""" Custom Type Error to stop the simulation. """
	pass

cdef class Particle:
	def __init__(self, 
				 Volume space, 
				 double E, 
				 Vector pos, 
				 double theta, 
				 double phi,
				 Vector ex,
				 Vector ey,
				 Vector ez,
				 bint simulate_secondary,
				 Volume current_region):


		self.simulate_secondary = simulate_secondary
		
		self.space = space
		self.pos   = pos
		self.current_region = current_region

		self.ex, self.ey, self.ez = ex, ey, ez

		self.theta = theta
		self.phi   = phi

		self.change_direction(cos(self.theta), 
							      self.phi)
		
		self.E = E
		self.update_coefs() #probably should be deffered to derived init


		#self.children = None




		
#    def __setattr__(self, name, value):
#        """Overloaded __setattr__: keep a record of all previous values of a given variable"""
#        if name not in self.__dict__:
#            self.__dict__[name + "_"] =  [value]
#            self.__dict__[name] = value
#        else:
#            self.__dict__[name + "_"] += [value]
##            self.__dict__[name] = value
#


	cdef int propagate(self) except -1:
 		#proposal
		cdef double L = -log(rand())/self.mu_tot


		#ask geometry for position and region
		cdef Volume region
		self.pos, region = self.current_region.getCrossings(self.pos, 
															self.ez,
															L)

		if region is self.current_region:
			return 0



		#change recrod before leaving region
		self.current_region = region
		#####################################

		if self.current_region.imp is 0:
			raise StopSimulation

		#change record before continuing with simulation
		self.current_region.N += 1
		#################################################


		self.update_coefs()
		self.propagate()
		return 1

	cdef public update_coefs(self):
		print("> Some particle tried to call update_coefs from particle.pyx!")
		raise StopSimulation

	cdef change_direction(self, double cos, double phi):
		'''Rotates the local frame.'''


		cdef Vector axis

		axis = self.ez
		self.ey = self.ey.rotateCos(axis, np.cos(phi))
		#self.ex = rotate(self.ex, axis, phi)

		axis = self.ey
		self.ez = self.ez.rotateCos(axis, cos)
		#self.ex = rotate(self.ex, axis, theta)



















































































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
	***NOTE TO SELF: Avoid using this plz.

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


def dist(object P, object Q):
	"""Find distance between two points."""
	return sum((P-Q)**2)**.5


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







# class Particle(CParticle):
#     """Python wrapper.
#     Cython does not allow keywords.
#     Might need to change this behaviour later though, I need the inits
#     to be as quick as possible.


#     Nevermind, I'm already in the process of changing this stuff.
#     Keyword args should be reserved to sources.


#     Maybe keep this for a user front end, idk"""
				
	
#     def __init__(self,
#                  space  = [],
#                  E      = 6.,
#                  pos    = [0., 0., 0.],
#                  theta  = 0.,
#                  phi    = 0.,
#                  ex = [1., 0., 0.], ey = [0., 1., 0.], ez = [0., 0., 1.],
#                  current_region = None,
#                  simulate_secondary = False):
		

#         pos = Vector(pos[0], pos[1], pos[2])
#         ex = Vector(ex[0], ex[1], ex[2])
#         ey = Vector(ey[0], ey[1], ey[2])
#         ez = Vector(ez[0], ez[1], ez[2])

		

#         super().__init__(space, E, pos, 
#                          theta, phi, 
#                          ex, ey, ez, 
#                          simulate_secondary, current_region)



#     def summary(self):
#         return {'pos': self.pos,
#                 'E'  : self.E,
#                 'ez' : self.ez}


	
#     def plot(self, tube_radius = 1, line_width = 1):
#         """Plot the path of the particle. (uses mayavi)"""
#         history = array(self.pos_)
#         X, Y, Z = history[:,0], history[:,1], history[:,2]
#         plot3d(X, Y, Z, tube_radius = tube_radius, line_width = line_width)
#         show()



#     def _add_plot(self, fig):
#         """
#         Don't use this. Just an attempt of using other modules for visualizing.
#         Keeping it here to remember later that there are other options.
#         (there is also paraview, but it's very time consuming to learn that)s
#         """
#         from vispy import scene, visuals
#         Plot3D = scene.visuals.create_visual_node(visuals.LinePlotVisual)


#         history = array(self.pos_)
#         X,Y,Z = history[:,0], history[:,1], history[:,2]
#         pos = c_[X, Y, Z]
		
#         Plot3D(pos, width = 2.0, color='red', edge_color='w', symbol='o',
#               face_color = (0.2, 0.2, 1, 0.8),
#               parent     = fig)
		
#     def add_plot(self, fig, plot_secondary = False):
#         """
#         Add trajectory to a mayavi figure. Used by a source instance to plot the
#         results of the simulation.
#         """
		
#         name = self.__class__.__name__
#         if name == 'Photon': color = (0, 0, 1)
#         else:                color = (1, 0, 0)
		
#         history = array(self.pos_)
#         X,Y,Z   = history[:,0], history[:,1], history[:,2]
#         plot3d(X, Y, Z, figure=fig, tube_radius=.8, color=color)
		
#         if plot_secondary is True:
#             for p in self.children_[1:]:
#                 p.add_plot(fig, plot_secondary = True)