cdef class Vector:
    cdef double x, y, z, real
    cdef Vector result
    cdef double norm(self)
    cdef Vector normalize(self)
    cdef Vector _mult(self, Vector q)
    cdef Vector _conj(self)
    cdef Vector rotateCos(self, Vector axis, double cos)
    cdef Vector rotateAngle(self, Vector axis, double theta)



    @staticmethod
    cdef Vector _new(double x, double y, double z)

    cdef void fastNORMALIZE(self)
    



    cdef Vector ADD(Vector self, Vector v2)
    
    cdef Vector SUB(Vector self, Vector v2)
    cdef Vector MUL(Vector self, double other)    
    
    cdef double DOT(Vector self, Vector v2)
    
    
    cdef double DIST(Vector self, Vector other)
            
    cdef Vector DIV(Vector self, double s)    
    cdef Vector XOR(Vector self, Vector q)