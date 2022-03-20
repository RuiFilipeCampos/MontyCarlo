

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