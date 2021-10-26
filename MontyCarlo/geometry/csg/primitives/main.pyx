


ctypedef (*map_t)(double3& pos, double[:] transformation)


cdef class Primitive(CSGvol):

    cdef map_t apply_transform
	cdef double[16] direct_transform
    cdef double[16] inverse_transform


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
            self.apply_transform = &apply_general_transform
        
        elif translated:
            self.apply_transform = &apply_translation
        
        elif rotated:
            self.apply_transform = &apply_rotation

		super(Primitive, self).__init__(*args, **kwargs)
