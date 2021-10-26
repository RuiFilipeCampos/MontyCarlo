


ctypedef (*map_t)(double3& pos, double[:] transformation)


cdef class Primitive(CSGvol):
    cdef map_t transform

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
