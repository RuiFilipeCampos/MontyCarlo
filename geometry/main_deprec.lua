# distutils: language = c++
# cython: profile=False

print(">>>>>   IMPORTING GEOMETRY")

from libc.math cimport fmin, fmax, sqrt , cos, pi
from libcpp.vector cimport vector

from numpy.math cimport INFINITY

#from ..tools.vectors cimport Vector

import numpy as np 

#from ..random.mixmax.interface cimport mixmax_engine

from numpy.random import randint


from numpy cimport ndarray


from libc.stdlib cimport malloc, free
cdef double MINdr = 0

from libcpp.list cimport list as cpplist
    
from cython.operator import dereference, preincrement



cdef struct double3:
    double x, y, z
    
cdef struct Ray:
    double3 O
    double3 D

cdef struct intersection:
    double3 pos
    double d



cdef class Volume:

    
    
    
    def __gt__(self, other):
        self.tr.b.x = -other[0]
        self.tr.b.y = -other[1]
        self.tr.b.z = -other[2]
        
        self.mesh.translate(other)
        return True
        
    def __lshift__(self, material):
        self.fill(material)
    
    def __enter__(self): return self
    def __exit__(self, *args, **kwargs):
        cdef Volume region
        cdef int i
        if self.Ninner > 0:
            self.inner = <void**>malloc(self.Ninner * sizeof(void*))
            for i, region in enumerate(self.tmp_inner):
                self.inner[i] = <void*> region
                

            
        
        if not (isinstance(self.outer, InfiniteVolume) or isinstance(self.outer, bInfiniteVolume) or isinstance(self, InfiniteVolume) or isinstance(self, bInfiniteVolume)):
            print("here:")
            proposed_mesh = self.mesh.boolean_cut(self.outer.mesh)
            print(">>>>!")
            if not proposed_mesh.number_of_cells == 0:
                self.mesh = proposed_mesh
       # self.generate_mesh(10_000)


    def set_array(self, arr):
        self.voxels = arr
    
    def get_arr(self):
        return self.voxels
    
    def setOuter(self, other):
        self.outer = other
    
    def fill(self, mat):
        self.material = mat 
        
    cdef bint move(self, STATE& state):
        raise RuntimeError("move called from virtual ")
    cdef double SDF(self, double3& pos):
        raise RuntimeError("SDF FROM VOL WAS CALLED")
        

    
    
cimport cython
@cython.boundscheck(False)
@cython.wraparound(False) 
@cython.initializedcheck(False)
@cython.cdivision(True)
cdef class BoundingVolume(Volume):
    
    def __init__(self):
        super(BoundingVolume, self).__init__()
    
    def __contains__(self, other):

        if isinstance(other, Volume):
            print(self, "in", other)
            self.eat(other)
            return True
        
        
        try:
            if len(other) == 3:
                for x in other:
                    try: int(x)
                    except TypeError:
                        raise TypeError("One of the elements in container is not a number.")
                        
                distance = self.SDF(other[0], other[1], other[2])
                if distance < 0: return True
                return False
    
            return RuntimeError("Container is not 3d vector. (lenght 3)")
        except TypeError:
            raise TypeError(f"Undefined behaviour for object of type '{type(other)}'.")
            
        



    
    #PYTHON INTERFACE - CONSTRUCTION OF BVH
    def eat(self, Volume other):
        print("1")
        if self.Ninner == 0:
            print("2")
            self.Ninner = 1
            print("3") 
            self.tmp_inner = [other]
        else:
            print("else")
            self.tmp_inner += [other]
            self.Ninner += 1
        print("other")
        other.setOuter(self)
        
    cdef Volume search_inner_asOUTER(self, STATE *state, void *to_ignore):
        
        cdef int i
        
        for i in range(self.Ninner):
            if self.inner[i] == to_ignore: 
                continue
        
            if (<Volume> self.inner[i]).SDF(state) < 0:
                return (<Volume> self.inner[i]).search_inner(state)
        
        return self
        
    
    cdef Volume search_inner(self, double x, double y, double z):
        cdef int i
        
        for i in range(self.Ninner):
            if (<Volume> self.inner[i]).SDF(x, y, z) < 0:
                return (<Volume> self.inner[i]).search_inner(x, y, z)
        
        return self
        
        
   
    cdef double find_safe_dr(self, STATE* state):
        cdef double  newD
        cdef double currentD = (<Volume> self.inner[0]).SDF(x, y, z)
        self.closest = 0
        cdef int i
       
        for i in range(1, self.Ninner):
           newD = fmax( (<Volume> self.inner[i]).SDF(x, y, z),  (<Volume> self.inner[i]).last )
           if newD < currentD:
               currentD = newD
               self.closest = i

        return currentD
    


            
        
        
    
    cdef bint move(self, double L, STATE *state):
        # cdef double dr =  fmin(self.find_safe_dr(self.p.x, self.p.y, self.p.z), L)

        # self.p.x += dr*self.p.ezx
        # self.p.y += dr*self.p.ezy
        # self.p.z += dr*self.p.ezz
        

        
        cdef double dr 
       
        cdef double minDIST0 = INFINITY
        cdef double minDIST1 = INFINITY
        cdef double testDIST
        
        
        while 1:
            
            # do ray trace if marked in previous step
            if state.ray_trace:
            
                
                # get intersections with the closest volume
                (<Volume> self.inner[self.closest]).intersect(state.pos, state.dire, self.intPTR)
                
                # no intersections, mark closest volume to be ignored, turn off ray trace, clear intersections
                if self.intPTR.size() == 0:
                    self.to_test[self.closest] = False
                    state.ray_trace = False
                    self.intPTR.clear()
                    continue
                
                # an intersection has ocurred, find the closest
                
                self.it = self.intPTR.begin()
          
                cdef int i

                for i in range(self.intpTR.size()):
                
                    if self.it.d < current_inters.d:
                        int0 = <intersection*> self.it
                        continue
                    elif self.it.d < int1.d:
                        int1 = <intersection*> self.it
                
                if self.closest == 0:
                    self.outer.find_inner_asOuter(state, self.inner[self.closest])
                else:
                    (<Volume> self.inner[self.closest]).find_inner(state)
                
                
                if state.current_region == self.inner[self.closest]:
                    ...
                
                    
                self.last = -INF
                state.pos = int0.pos
                state.last.int = int1
                
                
                state.last.region = self.inner[self.closest]
                
                return True
            
            
            dr = self.find_safe_dr(state)
            

            if dr < eps:
                state.ray_trace = True
                continue
                
            if dr > L:
                state.pos.x += L*state.dire.x
                state.pos.y += L*state.dire.y
                state.pos.z += L*state.dire.z
                
                self.last = -INF
                return False
            
            
            state.pos.x += dr*state.dire.x
            state.pos.y += dr*state.dire.y
            state.pos.z += dr*state.dire.z
            L -= dr

    
    
    cdef void* find_inner(self, double x, double y, double z):
        cdef int i
        for i in range(self.Ninner):
            if (<Volume> self.inner[i]).SDF(x, y, z) < 0:
                return (<Volume> self.inner[i]).find_inner(x, y, z)
            
        return <void*> self
        
       
    cdef bint check_exit(self):
        cdef double r = self.SDF(self.particle.x, self.particle.y, self.particle.z)
        if r < 0: return False
        
        
   
   
    # cdef Volume cross(self, particle):
        
    #     cdef double dL = self.find_save_step()
        
    #     cdef double r, dr
    #     r = 0
        
    #     while r < L:
    #         if self.SDF(particle.x, particle.y, particle.z) > 0:
    #             return self.outer.cross(particle)
            
    #         dr = self.find_safe_step()
    #         particle.move(dr)
    #         r += dr
        
        
        
    #     if self.SDF(x, y, z) > 0:
    #         return self.outer.cross(x, y, z)
            
            
    #     for i in range(self.Ninners):
    #         region = self.inner[i]
    #         if region.SDF(x, y, z) < 0:
    #             return region.search_inner(x, y, z)
        
    #     return self
        
        
        



    



cdef class EmptyVolume(Volume):
    
    cdef bint move(self, double L):
    
        
        self.p.move(L)

        
       
       
        L = self.SDF(self.p.x, self.p.y, self.p.z)

        if L < 0: return False
        if L > 0: self.move_to_surface()
        #self.p.move(-L + 1000*MINdr)
        #else: self.move_to_surface()

        self.p.current_region = <Volume> self.outer.find_inner(self.p.x, self.p.y, self.p.z)
        self.p.current_region.p = self.p
        return True
       
       
       
       
cdef class InfiniteVolume(Volume):
    cdef bint move(self, double L):
        self.p.move(L)
        #self.p.x += L*self.p.ezx
        #self.p.y += L*self.p.ezy
        #self.p.z += L*self.p.ezz
        return False
        
    
cdef class bInfiniteVolume(BoundingVolume):
    

    def __init__(self):
        super(bInfiniteVolume, self).__init__()

    cdef bint move(self, double L):
        
        
        cdef double dr
       
       
       

        while 1:
        
            dr = fmin(self.find_safe_dr(self.p.x, self.p.y, self.p.z), L)
            
            if dr == L:
                self.p.move(L)
                break
                
        
            if dr < MINdr:
                self.p.move(1000*MINdr)
                self.p.current_region = <Volume>(<Volume> self.inner[self.closest]).find_inner(self.p.x, self.p.y, self.p.z)
                self.p.current_region.p = self.p
                return True
            self.p.move(dr)
            L -= dr
            
        return False
        
        
    
    
    
       
cdef class VACCUM(Volume):
    pass
       
       
       
      
        
      
cdef struct NODE:
    cdef double L
    cdef NODE* next
    
cdef struct LIST:
    cdef NODE* start
    cdef int N #number of elements
    
    

    
## probly should be classified as empty union or bounding union
cdef class UNION(Volume):
    """
    SDF:
        max(SDF_a, SDF_b)
        
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
    def __init__(self, Volume A, Volume B):
        self.R = B
        self.L = A
        
    cdef double SDF(self, double3 *pos):
        
        
        
        ## <--- transform point
        
        return fmax(self.R.SDF(pos), self.L.SDF(pos))
    
    
    cdef void multithreaded_intersect(self, Ray ray, cpplist[double3] *intPTR):

        
        cdef int i
        for i in prange(2, ):
            if i == 0:
                self.L.intersect(ray, intPTR)
            
            if i == 1:
                self.R.intersect(ray, intPTR)
    
    
    cdef void intersect(self, Ray ray, cpplist[double3] *intPTR) nogil:
        
        #cdef Ray ray = deref(ray)
        

        
        ## <-- transform point
        
        ## <-- transform direction
        
        self.N = 0
        
        self.L.intersect(pos, dire, intPTR)
        self.R.intersect(pos, dire, intPTR)
        
        cdef cpplist[double3].iterator it = intPTR.begin()
        
        cdef int i
        
        #testing L aginst R
        for i in range(self.L.N):
            if sign(self.R.SDF(it) == 1:
                self.N += 1
                
                ## <-- inverse transform on point
                
                preincrement(it)
                continue
            
            intPTR.erase(it)
            preincrement(it)
            
        #testing R aginst L
        for i in range(self.R.N):
            if sign(self.L.SDF(it)) == 1:
                self.N += 1
                
                ## <-- inverse transform on point 
                
                preincrement(it)
                continue
            
            intPTR.erase(it)
            preincrement(it)
            
            
            
cdef class INTERSECTION(Volume):
    """
    SDF:
        min(SDF_a, SDF_b)
        
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
    def __init__(self, Volume L, Volume R):
        self.L = L
        self.R = R
        
    cdef double SDF(double3 *pos):
        return fmin(self.L.SDF(pos), self.R.SDF(pos))
    
    cdef void intersect(self, double3 pos, double3 dire, cpplist[double3] *intPTR):
        
        self.N = 0
 
        
        self.L.intersect(pos, dire, intPTR)
        self.R.intersect(pos, dire, intPTR)
        
        cdef cpplist[double3].iterator it = intPTR.begin()
        
        cdef int i
        
        #testing L aginst R
        for i in range(self.L.N):
            if sign(self.R.SDF(it) == -1:
                self.N += 1
                preincrement(it)
                continue
            intPTR.erase(it)
            preincrement(it)
            
        #testing R aginst L
        for i in range(self.R.N):
            if sign(self.L.SDF(it)) == -1:
                self.N += 1
                preincrement(it)
                continue
            intPTR.erase(it)
            preincrement(it)
            
        
                   
cdef class DIFFERENCE(Volume):
    """
    SDF:
        min(SDF_a, SDF_b)
        
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
    def __init__(self, Volume L, Volume R):
        self.L = L
        self.R = R
        
        
        
        
        
        
        
        
        
        
        
        

    cdef double SDF(double3 *pos):
        return fmin(self.L.SDF(pos), self.R.SDF(pos))

    
    
    
    

    
    cdef void intersect(self, double3 *pos, double3 *dire, cpplist[intersection] *intPTR):
        
        self.N = 0
 
        
        self.L.intersect(pos, dire, intPTR)
        self.R.intersect(pos, dire, intPTR)
        
        
        cdef cpplist[intersection].iterator it = intPTR.begin()
        
        cdef int i
        
        #testing L aginst R
        for i in range(self.L.N):
            if sign(self.R.SDF(it)) == 1:
                self.N += 1
                preincrement(it)
                continue
            
            intPTR.erase(it)
            preincrement(it)
            
        #testing R aginst L
        for i in range(self.R.N):
            if sign(self.L.SDF(it)) == -1:
                self.N += 1
                preincrement(it)
                continue
            
            intPTR.erase(it)
            preincrement(it)
            

            
            
        
        
        







cdef class eSphere(EmptyVolume):
    cdef double r
    def __init__(self, double radius):
        super(eSphere, self).__init__()

        self.r = radius
        self.r2 = radius*radius
        
        import pyvista as pv
        self.mesh = pv.Sphere(radius = radius)
        
    
    
    cdef double SDF(self, double3 pos):
        
        if self.tr:
            self.transform_point(&pos)

        return sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z) - self.r
    
    
    
    cdef void transform_point(self, double3* P):
        pass
    
    cdef void INVtransform_point(self, double3* P):
        pass
    
    
    cdef void intersect(self, double3 P, double3 D, cpplist[double3] *intPTR):
        
        if self.tr:
            self.transform_point(&P)
            self.transform_direction(&D)
        
        
        cdef double DELTA
        
        cdef double b = P.x*D.x + P.y*D.y + P.z*D.z
        
        DELTA = b*b - P.x*P.x - P.y*P.y - P.z*P.z + self.r2
        
        self.N = 0
        
        if DELTA <= 0:
            return
        
        self.N += 2
        
        DELTA = sqrt(DELTA)
        b *= -1
        cdef intersection I1
        I1.d = b - DELTA
      
        I1.pos.x = P.x + D.x*d
        I1.pos.y = P.y + D.y*d
        I1.pos.z = P.z + D.z*d


        cdef intersection I2
        I2.d = b + DELTA
      
        I2.pos.x = P.x + D.x*d
        I2.pos.y = P.y + D.y*d
        I2.pos.z = P.z + D.z*d

        if self.tr:
            self.INVtransform_point(&I1.pos)
            self.INVtransform_point(&I2.pos)
        
        intPTR.push_back(I1)
        intPTR.push_back(I2)
        
        
    

    
    
    def show(self):
        pass
        
    
    
cdef class bSphere(BoundingVolume):
    
    
    cdef double r
    def __init__(self, double radius):
        super(bSphere, self).__init__()

        import pyvista as pv
        self.mesh = pv.Sphere(radius = radius)
        
        self.r = radius
        
    cdef double SDF(self, double x, double y, double z):
        x *= self.tr.a.x
        y *= self.tr.a.y
        z *= self.tr.a.z
        
        x += self.tr.b.x
        y += self.tr.b.y
        z += self.tr.b.z
        return sqrt(x*x + y*y + z*z) - self.r
    

    
cdef class eBox(EmptyVolume):
    cdef double Rx, Ry, Rz
    
    def __init__(self, double Rx, double Ry, double Rz):
        super(eBox, self).__init__()
        import pyvista as pv
        self.mesh = pv.Box(bounds = (-Rx, Rx, -Ry, Ry, -Rz, Rz), quads = False, level = 2)

        self.Rx = Rx
        self.Ry = Ry
        self.Rz = Rz
        
    cdef double SDF(self, double x, double y, double z):
        x *= self.tr.a.x
        y *= self.tr.a.y
        z *= self.tr.a.z
        
        x += self.tr.b.x
        y += self.tr.b.y
        z += self.tr.b.z
        
        x = abs(x) - self.Rx
        y = abs(y) - self.Ry
        z = abs(z) - self.Rz
        
        return sqrt(fmax(x, 0)**2 + fmax(y, 0)**2 + fmax(z, 0)**2) + fmin(fmax(fmax(x, y), z), 0)
    
    
    def show(self):
        pass
        
    
    
cdef class bBox(BoundingVolume):
    cdef double Rx, Ry, Rz
    
    def __init__(self, double Rx, double Ry, double Rz):
        super(bBox, self).__init__()
        
        
        import pyvista as pv
        self.mesh = pv.Box(bounds = (-Rx, Rx, -Ry, Ry, -Rz, Rz), quads = False, level = 2)
        
        self.Rx = Rx
        self.Ry = Ry
        self.Rz = Rz
        
    cdef double SDF(self, double x, double y, double z):
        x *= self.tr.a.x
        y *= self.tr.a.y
        z *= self.tr.a.z
        
        x += self.tr.b.x
        y += self.tr.b.y
        z += self.tr.b.z
        
        
        x = abs(x) - self.Rx
        y = abs(y) - self.Ry
        z = abs(z) - self.Rz
        
        return sqrt(fmax(x, 0)**2 + fmax(y, 0)**2 + fmax(z, 0)**2) + fmin(fmax(fmax(x, y), z), 0)
