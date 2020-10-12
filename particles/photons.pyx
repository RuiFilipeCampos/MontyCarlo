# cython: profile=True



#Error messages (to be moved to its own module)
errorMSG1 = "Exhausted allowed number of iterations for rejection sampling."

#External Imports
from numpy import *
from numpy.random import rand, randint
#import pickle -> probly not needed any more?


#Local Imports
from .particle import StopSimulation
from .particle cimport Particle
from ..tools.vectors import *
# --  -- from . import electrons as e


#settings
from ..settings import __photonCUTOFF__



cdef object choose(list cumul, list items):

	cdef double r = rand()*cumul[-1]
	cdef int i = 0

	for i in range(len(cumul)):
		if cumul[i] > r:
			return items[i-1]



cdef class Photon(Particle):
	cdef public double k

	cdef public int N_coh, N_incoh, N_photo, N_pair, N_trip

	cdef public list  _interactions
	cdef public list  interactions

	cdef public object current_material

	cdef public object coherent, incoherent, photoelectric
	cdef public object pairproduction, tripletproduction, inter

	cdef public double mu_tot
	cdef public list   mu_cumul

	cdef public list X, Y, Z


	def __init__(self, space, current_region,
				 E     = 6.,
				 pos   = Vector(0., 0., 0.),
				 theta = 0.,
				 phi   = 0.,
				 ex    = Vector(1., 0., 0.), 
				 ey    = Vector(0., 1., 0.), 
				 ez    = Vector(0., 0., 1.),
				 simulate_secondary = False):

		self.k = E/0.511

		self._interactions = [self._coherent,
							  self._incoherent,
							  self._photoelectric,
							  self._pairproduction,
							  self._tripletproduction]

		self.X = []
		self.Y = []
		self.Z = []



		super().__init__(space, E, pos, 
						 theta, phi, 
						 ex, ey, ez, 
						 simulate_secondary,
						 current_region) 

	#@timer

	def run(self):
		self._run()



	cdef public _run(self):
		try:
			while 1:
				self.X += [self.pos.x]
				self.Y += [self.pos.y]
				self.Z += [self.pos.z]

				self.propagate()
				#print(self.mu_cumul)
				#self.mu_tot = 0.

				self.inter = choose(self.mu_cumul, self._interactions)

				self.inter(self)


		except StopSimulation:
			self.X += [self.pos.x]
			self.Y += [self.pos.y]
			self.Z += [self.pos.z]
			self.E = 0.511*self.k
			if self.simulate_secondary is True:
				if self.children is not None:
					for particle in self.children_[1:]:
						particle.simulate_secondary = True
						particle.run()



	cdef public update_coefs(self):
		"""
		Updates all references. Usually called when there is a region crossing.
		"""
		#print(self.current_region.material)
		self.current_material = self.current_region.material.photon
		
		#these references are used by their corresponding _fooInteraction method
		self.coherent          = self.current_material.coherent
		self.incoherent        = self.current_material.incoherent
		self.photoelectric     = self.current_material.photoelectric
		self.pairproduction    = self.current_material.pairproduction
		self.tripletproduction = self.current_material.tripletproduction


		#this list is used so I can iterate over all interactions
		#to create the cumulutative list
		self.interactions = [self.coherent, 
							 self.incoherent,
							 self.photoelectric,
							 self.pairproduction,
							 self.tripletproduction]
		
		self.update_mu()



	cdef public update_mu(self):
		"""
		Updates cross sections. Called when there is a region crossing or
		when the value of energy has changed.
		"""

		self.mu_tot   = 0.
		self.mu_cumul = [0.]

		self.E = self.k * 0.511

		for interaction in self.interactions:
			self.mu_tot   += interaction.CS(self.E)
			self.mu_cumul += [self.mu_tot]


	def summary(self):
		print(f"coh: {self.N_coh} \n incoh: {self.N_incoh}")
		print(f"pair: {self.N_pair} \n trip: {self.N_trip}")








































	cdef _coherent(self):
		"""
		- Samples angular deflection for the coherent scattering.
		- Changes direction accordingly.

		To do: Low energy stuff.
		"""

		self.N_coh += 1
		
		cdef double r, x2, cos
		cdef double x_max = 20.6074*2*self.k
		cdef double r_max = self.coherent.FF.cumul(x_max) #internals of this needs work
		
		while 1:
			#Sample x**2 from squared form factor (limited in (0, x_max**2))
			r  = rand()*r_max

			x2 = self.coherent.FF.invCum(r)

			#Get cos(theta) from x**2 and k = E/0.511MeV
			cos  = 1 - 0.5 * x2 / (20.6074*self.k)**2

			#Sample from thomson scattering.
			if rand() < (1+cos**2)*.5:
				self.change_direction(cos, 2*pi*rand())
				break


	cdef _incoherent(self):
		"""
		Incoherent sampling.
		"""

		self.N_incoh += 1
		
		cdef double t1 = log(1 + 2*self.k)
		cdef double t2 = 2*self.k*(1 + self.k)/(1 + 2*self.k)**2
		cdef double tau_min = 1/(1 + 2*self.k)
		cdef double tau, cos, N, D, sin2, x , T
		
		cdef int _
		for _ in range(10_000):

			#Sample fractional energy loss(tau)
			if rand() < t1/(t1 + t2): tau = tau_min**rand()
			else:                     tau = (tau_min**2 + rand()*(1 - tau_min**2))**.5

			#Calculate cosine from current energy and tau
			cos = (self.k + 1 - 1/tau)/self.k

			#constructing eqn in braces
			N = (1 - tau)*( (2 * self.k + 1) * tau - 1)
			D = self.k**2 * tau * (1 + tau**2)

			#constructing argument for Incoherent Form Factor
			sin2 = (.5 * (1 - cos))**.5
			x = self.k*sin2*2*20.6074
			
			#calculating proposal
			T = (1 - N/D)*self.incoherent.S(x)

			if rand() < T:
				self.change_direction(cos, 2*pi*rand())
				self.k = tau*self.k
				self.update_coefs()
				if self.k < __photonCUTOFF__: #note to self: make cut off adjustable
					raise StopSimulation
				
				break
		else: 
			print("k =", self.k)
			raise RuntimeError("Incoherent Scattering Sampling: " + errorMSG1)


	cdef _pairproduction(self):
		self.N_pair += 1
		
		#assert False
		raise StopSimulation

	cdef _photoelectric(self):
		self.N_photo += 1

		raise StopSimulation
		#first select the element that will be ionized
		cdef list P = [0.]
		cdef double Ptot = 0.
		cdef list elements = []

		for element in self.photoelectric:
			elements += [element]
			Ptot     += element(self.k*0.511)
			P        += [Ptot]

		element = choose(P, elements) #element that will be ionized

		#then select the active shell
		P = [0.]
		Ptot = 0.
		cdef shells = []
		for shell in element:
			shells += [shell]
			Ptot += shell(self.k*0.511)
			P    += [Ptot]

		shell = choose(P, shells)		
		raise StopSimulation



	def __repr__(self):
		return f"<Photon: pos = {self.pos}, ez = {self.ez}, k = {self.k}>"


	cdef _tripletproduction(self):
		self.N_trip += 1
		#raise StopSimulation
		raise StopSimulation


		
