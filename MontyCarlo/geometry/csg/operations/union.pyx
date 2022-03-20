
			   
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
