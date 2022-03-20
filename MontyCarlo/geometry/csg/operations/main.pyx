
cdef double delta = 1e-10


cdef class CSGop(CSGvol):
	cdef CSGvol R, L
	cdef double (*rule)(double, double)

	def __init__(self, CSGvol L, CSGvol R):
		self.L = L
		self.R = R

	cdef bint is_inside(self, double3 pos):
		raise RuntimeError("'is_inside' called from virtual Volume.BVH.CSGvol.CSGop")

	def translate(self, dx, dy, dz):
		if isinstance(self.L, Transform):
			self.L.primitive.translate(dx, dy, dz)
		else: self.L.translate(dx, dy, dz)

		if isinstance(self.R, Transform):
			self.R.primitive.translate(dx, dy, dz)
		else: self.R.translate(dx, dy, dz)


	def rotate(self, axis, angle):
		self.L.rotate(axis, angle)
		self.R.rotate(axis, angle)

	cdef intLIST intersect(self, double3& pos, double3& dire):
		raise RuntimeError("Called from virtual;")
		
		#return self.rule(self.L.intersect(pos, dire), self.R.intersect(pos, dire))



cdef class Subtraction(CSGop):

	"""
	SDF:  max(-SDF_L, SDF_R)
		
	RAY TRACE:
		
		   L/R     IN    OUT     BORDER  
		 -------- ---- -------- -------- 
		  IN       OUT      IN       BORDER      
		  OUT      OUT     OUT      OUT  
		  BORDER   OUT  BORDER      
		  
		  1 = out
		  -1 = in
		  0 = border
			
	"""
	def __init__(self, CSGvol L, CSGvol R):
		# child nodes
		super(Subtraction, self).__init__(L, R)
		#self.mesh = L.mesh - R.mesh




	cdef double SDF(self, double3 pos):
		return fmax(self.L.SDF(pos), -self.R.SDF(pos))

	cdef bint is_inside(self, double3& pos):
		return self.L.is_inside(pos) and (not self.R.is_inside(pos))

	def __repr__(self):
		return "<Subtraction>"

	cdef intLIST intersect(self, double3& pos, double3& dire):
		IF VERBOSE: print("SUBTRACTING: \n --left-- \n")

		cdef intLIST L = self.L.intersect(pos, dire)

		if L.size() == 0:
			return L

		IF VERBOSE: print("\n --right-- \n")

		cdef intLIST R = self.R.intersect(pos, dire)

		if R.size() == 0:
			return L

		return intMinus(L, R)








			   
cdef class Union(CSGop):

	"""
	SDF:
		min(SDF_a, SDF_b)
		
	RAY TRACE:
		
		   L/R     IN    OUT     BORDER  
		 -------- ---- -------- -------- 
		  IN       IN   IN       IN      
		  OUT      IN   OUT      BORDER  
		  BORDER   IN   BORDER      
		  
		  1 = out
		  -1 = in
		  0 = border
			
	"""
	def __init__(self, CSGvol L, CSGvol R):
		# child nodes
		super(Union, self).__init__(L, R)
		#self.mesh = L.mesh + R.mesh

	def __repr__(self):
		return "<Union>"

	cdef double SDF(self, double3 pos):
		return fmin(self.L.SDF(pos), self.R.SDF(pos))

	cdef intLIST intersect(self, double3& pos, double3& dire):
		IF VERBOSE: print("UNION: \n --left-- \n")

		cdef intLIST L = self.L.intersect(pos, dire)

		if L.size() == 0:
			return self.R.intersect(pos, dire)

		IF VERBOSE: print("\n --right-- \n")

		cdef intLIST R = self.R.intersect(pos, dire)

		if R.size() == 0:
			return L

		return intPlus(L, R)

	cdef bint is_inside(self, double3& pos):
		return self.L.is_inside(pos) or self.R.is_inside(pos)


cdef class Intersection(CSGop):

	"""
	SDF:
		max(SDF_a, SDF_b)
		
	RAY TRACE:
		
		   L/R     IN    OUT     BORDER  
		 -------- ---- -------- -------- 
		  IN       IN      OUT       BORDER      
		  OUT      OUT     OUT      OUT  
		  BORDER   BORDER  OUT      
		  
		  1 = out
		  -1 = in
		  0 = border
			
	"""
	def __init__(self, CSGvol L, CSGvol R):
		# child nodes
		super(Intersection, self).__init__(L, R)
		#self.mesh = L.mesh.boolean_cut(R.mesh)


	def __repr__(self):
		return "<Intersection>"
 
	cdef double SDF(self, double3 pos):
		return fmax(self.L.SDF(pos), self.R.SDF(pos))

	cdef intLIST intersect(self, double3& pos, double3& dire):
		IF VERBOSE: print("INTERSECTING: ")
		cdef intLIST L = self.L.intersect(pos, dire)

		if L.size() == 0:
			return L

		cdef intLIST R = self.R.intersect(pos, dire)
		if R.size() == 0:
			return R

		return intIntersect(L, R)

	cdef bint is_inside(self, double3& pos):
		return self.L.is_inside(pos) and self.R.is_inside(pos)