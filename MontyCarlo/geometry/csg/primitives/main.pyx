

import numpy as np



cdef void apply_general_transform(double3& pos, double[:] transform):
    cdef double3 _pos = pos
    pos.x = transform[0]*_pos.x + transform[1]*_pos.y + transform[2] *_pos.z  + transform[3]
    pos.y = transform[4]*_pos.x + transform[5]*_pos.y + transform[6] *_pos.z  + transform[7]
    pos.z = transform[8]*_pos.x + transform[9]*_pos.y + transform[10]*_pos.z  + transform[11]

cdef void apply_translation(double3& pos, double[:] transform):
    pos.x += transform[3]
    pos.y += transform[7]
    pos.z += transform[11]

cdef void apply_rotation(double3& pos, double[:] transform):
    cdef double3 _pos = pos
    pos.x = transform[0]*_pos.x + transform[1]*_pos.y + transform[2] *_pos.z
    pos.y = transform[4]*_pos.x + transform[5]*_pos.y + transform[6] *_pos.z
    pos.z = transform[8]*_pos.x + transform[9]*_pos.y + transform[10]*_pos.z




cdef class Primitive(CSGvol):
	def __init__(self, *args, **kwargs):

        tmp_direct_transform = np.array(
            [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
            ], dtype = np.float64,
        )


        self.translated = False
        self.rotated = False

        if kwargs['transforms']:
            for transformation in kwargs['transforms']:

                if transformation[0] == "translate":
                    self.translated = True
                    tmp_direct_transform[3]  += transformation[1]
                    tmp_direct_transform[7]  += transformation[2]
                    tmp_direct_transform[11] += transformation[3]

                elif transformation[0] == "rotate":
                    self.rotated = True
                    axis = np.array([
                        transformation[1][0], 
                        transformation[1][1], 
                        transformation[1][2], 
                    ], 
                    dtype = np.float64 )

                    axis = axis/np.sqrt(np.sum(axis**2))

                    ux = axis[0]
                    uy = axis[1]
                    uz = axis[2]

                    angle = transformation[2]

                    rotation_matrix = np.zeros((4,4))

                    _cos = cos(angle)
                    oneMcos = 1 - _cos
                    _sin = sqrt(1 - _cos*_cos)

                    rotation_matrix[0, 0] = _cos + ux*ux*oneMcos
                    rotation_matrix[0, 1] = ux*uy*oneMcos - uz*_sin
                    rotation_matrix[0, 2] = ux*uz*oneMcos + uy*_sin

                    rotation_matrix[1, 0] = uy*ux*oneMcos + uz*_sin
                    rotation_matrix[1, 1] = _cos + uy*uy*oneMcos
                    rotation_matrix[1, 2] = uy*uz*oneMcos - ux*_sin

                    rotation_matrix[2, 0] = uz*ux*oneMcos-uy*_sin
                    rotation_matrix[2, 1] = uz*uy*oneMcos + ux*_sin
                    rotation_matrix[2, 2] = _cos + uz*uz*oneMcos

                    rotation_matrix[3, 3] = 1

                    tmp_direct_transform = rotation_matrix @ tmp_direct_transform

        tmp_inverse_transform = np.linalg.inv(tmp_direct_transform)
        
        cdef int i, y, x

        i = 0
        for y in range(4):
            for x in range(4):
                self.direct_transform[i]  = tmp_direct_transform[x, y]
                self.inverse_transform[i] = tmp_inverse_transform[x, y]
                i += 1

        self.set_map()


            

		super(Primitive, self).__init__(*args, **kwargs)

    cdef void set_map(Primitive self):
        
        if self.translated and self.rotated:
            self.apply_transform = &apply_general_transform
            self.has_transform = True
        
        elif self.translated:
            self.apply_transform = &apply_translation
            self.has_transform = True
        
        elif self.rotated:
            self.apply_transform = &apply_rotation
            self.has_transform = True
        else:
            self.has_transform = False



def Translate(self, *args, **kwargs):

    cdef double x, y, z
    x = kwargs['x']
    y = kwargs['y']
    z = kwargs['z']

    for volume in args:
        (<Primitive> volume).direct_transform[3] += x
        (<Primitive> volume).inverse_transform[3] -= x

        (<Primitive> volume).direct_transform[3] += y
        (<Primitive> volume).inverse_transform[3] -= y

        (<Primitive> volume).direct_transform[3] += z
        (<Primitive> volume).inverse_transform[3] -= z


        (<Primitive> volume).translated = True
        (<Primitive> volume).set_map()         

        Translate(*list(volume))


def Rotate(self, *args, axis = (0, 0, 1), angle=0):
    axis = np.array([
        transformation[1][0], 
        transformation[1][1], 
        transformation[1][2], 
    ], 
    dtype = np.float64 )

    axis = axis/np.sqrt(np.sum(axis**2))

    ux = axis[0]
    uy = axis[1]
    uz = axis[2]

    angle = transformation[2]

    rotation_matrix = np.zeros((4,4))

    _cos = cos(angle)
    oneMcos = 1 - _cos
    _sin = sqrt(1 - _cos*_cos)

    rotation_matrix[0, 0] = _cos + ux*ux*oneMcos
    rotation_matrix[0, 1] = ux*uy*oneMcos - uz*_sin
    rotation_matrix[0, 2] = ux*uz*oneMcos + uy*_sin

    rotation_matrix[1, 0] = uy*ux*oneMcos + uz*_sin
    rotation_matrix[1, 1] = _cos + uy*uy*oneMcos
    rotation_matrix[1, 2] = uy*uz*oneMcos - ux*_sin

    rotation_matrix[2, 0] = uz*ux*oneMcos-uy*_sin
    rotation_matrix[2, 1] = uz*uy*oneMcos + ux*_sin
    rotation_matrix[2, 2] = _cos + uz*uz*oneMcos

    rotation_matrix[3, 3] = 1

    tmp_direct_transform = rotation_matrix @ tmp_direct_transform

