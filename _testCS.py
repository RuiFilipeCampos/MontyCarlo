
from tools.vectors cimport Vector
from geometry.CSG import Sphere
from particles.electrons cimport Electron



from time import perf_counter


print("SETTING UP")

space = Sphere(10_000)



cdef Electron el

cdef Vector pos, ex, ey, ez

pos = Vector(0, 0, 0)
ex = Vector(1, 0, 0)
ey = Vector(0, 1, 0)
ez = Vector(0, 0, 1)

cdef double theta, phi
theta = 0
phi = 0
cdef double t0, tf
cdef int i



cpdef void perf(mat):
    space.fill(mat)
    
    print("DOING TEN RUNS FOR JIT")
    for i in range(10):
        el = Electron._new(space,
                           space,
                           50e3, 
                           pos, 
                           theta, 
                           phi,
                           ex,
                           ey,
                           ez,
                           False,
                           100.)
        el._run()
           
    
    print("START")
    t0 = perf_counter()
    for i in range(10_000):
        el = Electron._new(space,
                           space,
                           50e3, 
                           pos, 
                           theta, 
                           phi,
                           ex,
                           ey,
                           ez,
                           False,
                           100.)
        el._run()
        
    tf = perf_counter()  
    
    print(tf - t0)
    