



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