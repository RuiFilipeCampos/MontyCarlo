





cdef class Primitive(CSGvol):
    cdef (*transform)(double3 pos, double[:] matrix)
    cdef double[:] matrix
    cdef double[:] inverse_matrix

	def __init__(self, *args, **kwargs):

        self.matrix = identity()

        translated = False
        rotated = False

        if kwargs['transforms']:
            for transformation in kwargs['transforms']:

                if transformation[0] == "translate":
                    translated = True
                    pass

                if transformation[0] == "rotate":
                    rotated = True
                    pass

        if translated and rotated:
            self.transform = general_transform
        
        elif translated:
            self.transform = translate
        
        elif rotated:
            self.transform = rotate


		super(Primitive, self).__init__(*args, **kwargs)

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