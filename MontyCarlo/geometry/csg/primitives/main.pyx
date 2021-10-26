




cdef class Primitive(CSGvol):
	cdef Transform tr


	def __init__(self):
		super(Primitive, self).__init__()
		self.tr = Identity(self)

	def translate(self, dx, dy, dz):
		self.tr = self.tr.translate(dx, dy, dz)
		#self.mesh.translate([dx, dy, dz])

	def rotate(self, axis, angle):
		self.tr = self.tr.rotate(axis, angle)

	@property
	def matrix(self):
		return self.tr.matrix

	@property
	def inv_matrix(self):
		return self.tr.inv_matrix