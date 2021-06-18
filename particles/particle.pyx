# distutils: language = c++
print(">>>>>   IMPORTING PARTICLE")

DEF FULL_RECORD = True

#from numpy.math cimport INFINITY
from ..settings import DEBUG
from libc.math cimport sqrt


#External Imports

import numpy as np
from libc.math cimport sin, cos, log, pow



cimport cython
cimport numpy as cnp

import numpy.random

rand = numpy.random.rand
#from pyquaternion import Quaternion
#import bisect

#Internal Imports
from ..tools.performance import timer
from ..geometry.main cimport Volume




#cdef object choose(list cumul, list items):
#
#    cdef double r = rand()*cumul[-1]
#    cdef int i = 0
#
#
#    for i in range(len(cumul)):
#        if cumul[i] > r:
#            return items[i-1]


#class StopSimulation(Exception):
 #   """
#    Custom Type Error to stop the simulation.
#    """
 #   
 #   pass




# cdef struct double3:
#     cdef double x
#     cdef double y
#     cdef double z
    

# cdef struct PARTICLE_CONTAINER:
#     cdef deque[void*] photons
#     cdef deque[void*] electrons
#     cdef deque[void*] positrons



# cdef struct STATE:
#     # geometric state
#     cdef double3 pos
#     cdef double3 ez
#     cdef double3 ey
#     cdef void* cr   #current region
    
#     # physics state
#     cdef double E
    
#     # multi processing state
#     cdef int PROCESS_INDEX
    
#cdef void double move(double x, double y, double z)
    
from ..types cimport double3

cdef struct STATE:
    mixmax_engine* genPTR
    void *current_region
    double3 pos
    double3 dire
    double3 axis
    double E
    double L 
    double last_displacement

from collections import deque

from cython.operator import dereference as deref
from cython.operator import preincrement as inc


@cython.boundscheck(False)
@cython.initializedcheck(False)
@cython.cdivision(True)
@cython.freelist(100_000_000)
cdef class Particle:

    cdef inline double ENERGY(self):
        raise RuntimeError("ENERGY METHOD CALLED FROM Particles")


    cpdef get_record_pos(self):
      cdef vector[STATE].iterator it = self.state_record.begin()
      
      traj = deque()
      cdef int i
      for i in range(self.state_record.size()):
        traj.append(deref(it).pos)
        inc(it)
      return traj

    cdef void record(self):
        """
        Add current position and energy to the track record.
        """
        
        IF FULL_RECORD:
            self.state_record.push_back(self.state)
        #ELSE:
        #    self.pos_record.push_back(self.state.pos)
        #    self.E_record.push_back(self.state.E)

        
    def getEnergy(self):
        return self.state.E
        
    cdef void move(self, double L):
        # self.state.pos.x += self.state.EZ.x*L
        # self.state.pos.y += self.state.EZ.y*L
        # self.state.pos.z += self.state.EZ.z*L

        self.state.pos.x += self.state.dire.x*L
        self.state.pos.y += self.state.dire.y*L
        self.state.pos.z += self.state.dire.z*L
        
    cdef inline void deposit(self):
        pass
        #self.current_region.voxels[self.index] += self.ENERGY()
    
  
    cdef double rsqrt(self, double x):
        x = 1 - x
        cdef double y = 1 + x*.5  # 0th and 1st order
        x *= x              
        y += x * 0.375            # 2nd order
        x *= x
        y += x * 0.3125
        x *= x
        y += x * 0.2734375
        return y
    
    cdef void normalize(self):
        # cdef double x = self.eyx**2 + self.eyy**2 + self.eyz**2 
        # cdef double norm
  
            
        # x = 1 - x
        # norm = 1 + x*.5  # 0th and 1st order
        # x *= x              
        # norm += x * 0.375            # 2nd order
        # x *= x
        # norm += x * 0.3125
        # x *= x
        # norm += x * 0.2734375
        
        # self.eyx *= norm
        # self.eyy *= norm
        # self.eyz *= norm
        cdef double x
        cdef double norm
        #x = (self.ezx**2 + self.ezy**2 + self.ezz**2 )
        
        x = 1 - (self.state.dire.x*self.state.dire.x + self.state.dire.y*self.state.dire.y + self.state.dire.z*self.state.dire.z )
        norm = 1 + x*.5  # 0th and 1st order
      #  x *= x              
      #  norm += x * 0.375            # 2nd order
     #   x *= x
       # norm += x * 0.3125
      #  x *= x
       # norm += x * 0.2734375
        
        self.state.dire.x *= norm
        self.state.dire.y *= norm
        self.state.dire.z *= norm
        
    cdef void invert_axis(self):
        self.state.axis.x *= -1
        self.state.axis.y *= -1
        self.state.axis.z *= -1


    cdef void invert_dire(self):
        self.state.dire.x *= -1
        self.state.dire.y *= -1
        self.state.dire.z *= -1


    cdef void throwAZIMUTH(self):
        cdef double x, y, a, w1, w2, w3

        while True:
            while True:
                x = 2*self.state.genPTR.get_next_float() - 1
                y = 2*self.state.genPTR.get_next_float() - 1
                a = x**2 + y**2
                if a < 1:
                    break


            w1 = 1 - 2*a 
            a = 2 * sqrt(1 - a)
            w2 = x*a
            w3 = y*a
            a =  w1*self.state.dire.x +  w2*self.state.dire.y + w3*self.state.dire.z
            #w and ez are unit vectors, this should be enough to guarentee that they are not in the same direc
            if a < 1: 
                
                x = 1/sqrt(1 - a**2)
                self.state.axis.x = (w1 - a*self.state.dire.x)*x
                self.state.axis.y = (w2 - a*self.state.dire.y)*x
                self.state.axis.z = (w3 - a*self.state.dire.z)*x
                return
                
                # x = 1/(self.eyx**2 + self.eyy**2 + self.eyz**2)**.5
                
                # self.eyx *= x
                # self.eyy *= x
                # self.eyz *= x
                

                 
        # else: 
        #     import time
            
        #     print("22222111ASdsasadjlkdsja lksadalkdsfjlkdasf lkjdsahf kjsad hfkjda hsflkj hdsaflkj hadlkjf halskjd")
            
        #     time.sleep(1000000)
        #     raise RuntimeError("Exceeded max")
    
    
    
    cdef void rotateTHETA(self, double cos):
        # spent an entire afternoon looking for a serious bug
        # as it happens, a 1.0000000000000002 value of cos will result in NaN values on the rotated vector
        
        
        if cos > 1 :
            cos =  1
        elif cos < -1: 
            cos = -1
        
        

        
        cdef double cos2sq = .5*(1 + cos)
        cdef double cos2 = sqrt(cos2sq)
        
        cdef double sin2sq = 1 - cos2sq
        cdef double sin2 = sqrt(sin2sq)
        cdef double sin2cos2 = 2*sin2*cos2
        
        cdef double eyx2 = self.state.axis.x**2
        cdef double eyy2 = self.state.axis.y**2
        cdef double eyz2 = self.state.axis.z**2
        
        cdef double ezx = self.state.dire.x
        cdef double ezy = self.state.dire.y
        cdef double ezz = self.state.dire.z
        
        #### calculated formulas using sympy's quaternions
        #### copied it here and organized best I could for performance
        self.state.dire.x = \
        cos2sq* ezx + \
        sin2sq* (eyx2 * ezx \
                + 2*self.state.axis.x  * (self.state.axis.y*ezy + self.state.axis.z*ezz) \
                - ezx * (eyy2  + eyz2))         \
        + sin2cos2*(self.state.axis.y*ezz - self.state.axis.z*ezy)
        
        

        self.state.dire.y = \
        cos2sq*ezy + \
        sin2sq*( ezy*( -eyx2 \
                + eyy2 \
                - eyz2)
                + 2*self.state.axis.y*(self.state.axis.x*ezx \
                            + self.state.axis.z*ezz)) \
        + sin2cos2*(-self.state.axis.x*ezz + self.state.axis.z*ezx)

        self.state.dire.z = cos2sq*ezz + \
        sin2sq*(- ezz*   (eyx2 + eyy2 - eyz2) \
                + 2*self.state.axis.z*(self.state.axis.x*ezx  + self.state.axis.y*ezy)) \
        + sin2cos2*(self.state.axis.x*ezy - self.state.axis.y*ezx)


       # print("")
      #  print("ey", self.eyx, self.eyy, self.eyz)
       # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
      #  print("ez", self.ezx, self.ezy, self.ezz)
      #  print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )




# ('ey', 0.304192241115606, 0.7845300950417162, 0.5403513768085548)
# ('norm:', 1.0)
# ('ez', 0.09865239061376316, 0.5382383859203108, -0.8369988923219004)
# ('norm:', 1.0)


        #cdef double x
 
        #x = (self.ezx**2 + self.ezy**2 + self.ezz**2 )
        
        
        #nevermind the names, just using already allocated variables
        cos2sq = 1 - (self.state.dire.x**2 + self.state.dire.y**2 + self.state.dire.z**2 )
        ezz = 1 + cos2sq*.5  # 0th and 1st order
        cos2sq *= cos2sq              
        ezz += cos2sq * 0.375            # 2nd order
        cos2sq *= cos2sq
        ezz += cos2sq * 0.3125
        cos2sq *= cos2sq
        ezz += cos2sq * 0.2734375
        
        self.state.dire.x *= ezz
        self.state.dire.y *= ezz
        self.state.dire.z *= ezz
        
        
        # print("")
        # print("ey", self.eyx, self.eyy, self.eyz)
        # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
        # print("ez", self.ezx, self.ezy, self.ezz)
        # print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )
        
        
        # print("rotateTHETA out")
        





############################################### prev versionof rotate


  ####  cdef void rotateTHETA(self, double cos):
  ####      # spent an entire afternoon looking for a serious bug
  ####      # as it happens, a 1.0000000000000002 value of cos will result in NaN values on the rotated vector
  ####      
####
  ####      if cos > 1 :
  ####          cos =  1
  ####      elif cos < -1: 
  ####          cos = -1
  ####      
  ####      
####
  ####      
  ####      cdef double cos2sq = .5*(1 + cos)
  ####      cdef double cos2 = sqrt(cos2sq)
  ####      
  ####      cdef double sin2sq = 1 - cos2sq
  ####      cdef double sin2 = sqrt(sin2sq)
  ####      cdef double sin2cos2 = 2*sin2*cos2
  ####      
  ####      cdef double eyx2 = self.eyx**2
  ####      cdef double eyy2 = self.eyy**2
  ####      cdef double eyz2 = self.eyz**2
  ####      
  ####      cdef double ezx = self.ezx
  ####      cdef double ezy = self.ezy
  ####      cdef double ezz = self.ezz
  ####      
  ####      #### calculated formulas using sympy's quaternions
  ####      #### copied it here and organized best I could for performance
  ####      self.ezx = \
  ####      cos2sq* ezx + \
  ####      sin2sq* (eyx2 * ezx \
  ####              + 2*self.eyx  * (self.eyy*ezy + self.eyz*ezz) \
  ####              - ezx * (eyy2  + eyz2))         \
  ####      + sin2cos2*(self.eyy*ezz - self.eyz*ezy)
  ####      
  ####      
####
  ####      self.ezy = \
  ####      cos2sq*ezy + \
  ####      sin2sq*( ezy*( -eyx2 \
  ####              + eyy2 \
  ####              - eyz2)
  ####              + 2*self.eyy*(self.eyx*ezx \
  ####                          + self.eyz*ezz)) \
  ####      + sin2cos2*(-self.eyx*ezz + self.eyz*ezx)
####
  ####      self.ezz = cos2sq*ezz + \
  ####      sin2sq*(- ezz*   (eyx2 + eyy2 - eyz2) \
  ####              + 2*self.eyz*(self.eyx*ezx  + self.eyy*ezy)) \
  ####      + sin2cos2*(self.eyx*ezy - self.eyy*ezx)
####
####
  ####     # print("")
  ####    #  print("ey", self.eyx, self.eyy, self.eyz)
  ####     # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
  ####    #  print("ez", self.ezx, self.ezy, self.ezz)
  ####    #  print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )
####
####
####
####
# ####('ey', 0.304192241115606, 0.7845300950417162, 0.5403513768085548)
# ####('norm:', 1.0)
# ####('ez', 0.09865239061376316, 0.5382383859203108, -0.8369988923219004)
# ####('norm:', 1.0)
####
####
  ####      #cdef double x
 ####
  ####      #x = (self.ezx**2 + self.ezy**2 + self.ezz**2 )
  ####      
  ####      
  ####      #nevermind the names, just using already allocated variables
  ####      cos2sq = 1 - (self.ezx**2 + self.ezy**2 + self.ezz**2 )
  ####      ezz = 1 + cos2sq*.5  # 0th and 1st order
  ####      cos2sq *= cos2sq              
  ####      ezz += cos2sq * 0.375            # 2nd order
  ####      cos2sq *= cos2sq
  ####      ezz += cos2sq * 0.3125
  ####      cos2sq *= cos2sq
  ####      ezz += cos2sq * 0.2734375
  ####      
  ####      self.ezx *= ezz
  ####      self.ezy *= ezz
  ####      self.ezz *= ezz
  ####      
  ####      
  ####      # print("")
  ####      # print("ey", self.eyx, self.eyy, self.eyz)
  ####      # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
  ####      # print("ez", self.ezx, self.ezy, self.ezz)
  ####      # print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )
  ####      
  ####      
  ####      # print("rotateTHETA out")
  ####      
#################################################

















    cdef void rotateTHETAvers2(self, double cos):
        # spent an entire afternoon looking for a serious bug
        # as it happens, a 1.0000000000000002 value of cos will result in NaN values on the rotated vector
        if cos > 1 :  
            cos =  1
        elif cos < -1: 
            cos = -1
        
        
        
        cdef double cos2sq = .5*(1 + cos)
        cdef double cos2 = sqrt(cos2sq)
        
        cdef double sin2sq = 1 - cos2sq
        cdef double sin2 = sqrt(sin2sq)
        cdef double sin2cos2 = 2*sin2*cos2
        
        cdef double eyx2 = self.state.axis.x**2 # eyx**2
        cdef double eyy2 = self.state.axis.y**2 # eyy**2
        cdef double eyz2 = self.state.axis.z**2 # eyz**2
        
        cdef double ezx = self.state.dire.x # ezx
        cdef double ezy = self.state.dire.y  #ezy
        cdef double ezz = self.state.dire.z # ezz
        
        #### calculated formulas using sympy's quaternions
        #### copied it here and organized best I could for performance
        self.ezx = \
        cos2sq* ezx + \
        sin2sq* (eyx2 * ezx \
                + 2*self.eyx  * (self.eyy*ezy + self.eyz*ezz) \
                - ezx * (eyy2  + eyz2))         \
        + sin2cos2*(self.eyy*ezz - self.eyz*ezy)
        
        

        self.ezy = \
        cos2sq*ezy + \
        sin2sq*( ezy*( -eyx2 \
                + eyy2 \
                - eyz2)
                + 2*self.eyy*(self.eyx*ezx \
                            + self.eyz*ezz)) \
        + sin2cos2*(-self.eyx*ezz + self.eyz*ezx)

        self.ezz = cos2sq*ezz + \
        sin2sq*(- ezz*   (eyx2 + eyy2 - eyz2) \
                + 2*self.eyz*(self.eyx*ezx  + self.eyy*ezy)) \
        + sin2cos2*(self.eyx*ezy - self.eyy*ezx)


       # print("")
      #  print("ey", self.eyx, self.eyy, self.eyz)
       # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
      #  print("ez", self.ezx, self.ezy, self.ezz)
      #  print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )




# ('ey', 0.304192241115606, 0.7845300950417162, 0.5403513768085548)
# ('norm:', 1.0)
# ('ez', 0.09865239061376316, 0.5382383859203108, -0.8369988923219004)
# ('norm:', 1.0)


        #cdef double x
 
        #x = (self.ezx**2 + self.ezy**2 + self.ezz**2 )
        
        
        #nevermind the names, just using already allocated variables
        cos2sq = 1 - (self.ezx**2 + self.ezy**2 + self.ezz**2 )
        ezz = 1 + cos2sq*.5  # 0th and 1st order
        cos2sq *= cos2sq              
        ezz += cos2sq * 0.375            # 2nd order
        cos2sq *= cos2sq
        ezz += cos2sq * 0.3125
        cos2sq *= cos2sq
        ezz += cos2sq * 0.2734375
        
        self.ezx *= ezz
        self.ezy *= ezz
        self.ezz *= ezz
        
        
        # print("")
        # print("ey", self.eyx, self.eyy, self.eyz)
        # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
        # print("ez", self.ezx, self.ezy, self.ezz)
        # print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )
        
        
        # print("rotateTHETA out")
        










    cdef void rotateAZIMUTH(self, double cos):
        # spent an entire afternoon looking for a serious bug
        # as it happens, a 1.0000000000000002 value of cos will result in NaN values on the rotated vector
        if cos > 1 :  
            
            cos =  1
        elif cos < -1: 
            
            cos = -1
        
        
        #cos = (<double> 1e-15 * (<int> (1e15* cos)))
        cdef double tempEYx = self.eyx
        cdef double tempEYy = self.eyy
        cdef double tempEYz = self.eyz
        
        cdef double tempEZx = self.ezx
        cdef double tempEZy = self.ezy
        cdef double tempEZz = self.ezz

        
        self.eyx = tempEZx
        self.eyy = tempEZy
        self.eyz = tempEZz
        
        self.ezx = tempEYx
        self.ezy = tempEYy
        self.ezz = tempEYz





        #print(f"rotateTHETA in -- cos = {cos}")
        
        #cdef double eyx = self.eyx
        #cdef double eyy = self.eyy
        #cdef double eyz = self.eyz
        #print("ey", self.eyx, self.eyy, self.eyz)
       # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
       # print("ez", self.ezx, self.ezy, self.ezz)
      #  print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )
        
        cdef double cos2sq = .5*(1 + cos)
        cdef double cos2 = sqrt(cos2sq)
        
        cdef double sin2sq = 1 - cos2sq
        cdef double sin2 = sqrt(sin2sq)
        cdef double sin2cos2 = 2*sin2*cos2
        
        cdef double eyx2 = self.eyx**2
        cdef double eyy2 = self.eyy**2
        cdef double eyz2 = self.eyz**2
        
        cdef double ezx = self.ezx
        cdef double ezy = self.ezy
        cdef double ezz = self.ezz
        
        #### calculated formulas using sympy's quaternions
        #### copied it here and organized best I could for performance
        self.ezx = \
        cos2sq* ezx + \
        sin2sq* (eyx2 * ezx \
                + 2*self.eyx  * (self.eyy*ezy + self.eyz*ezz) \
                - ezx * (eyy2  + eyz2))         \
        + sin2cos2*(self.eyy*ezz - self.eyz*ezy)
        
        

        self.ezy = \
        cos2sq*ezy + \
        sin2sq*( ezy*( -eyx2 \
                + eyy2 \
                - eyz2)
                + 2*self.eyy*(self.eyx*ezx \
                            + self.eyz*ezz)) \
        + sin2cos2*(-self.eyx*ezz + self.eyz*ezx)

        self.ezz = cos2sq*ezz + \
        sin2sq*(- ezz*   (eyx2 + eyy2 - eyz2) \
                + 2*self.eyz*(self.eyx*ezx  + self.eyy*ezy)) \
        + sin2cos2*(self.eyx*ezy - self.eyy*ezx)


       # print("")
      #  print("ey", self.eyx, self.eyy, self.eyz)
       # print("norm:", self.eyx**2 + self.eyy**2 + self.eyz**2 )
      #  print("ez", self.ezx, self.ezy, self.ezz)
      #  print("norm:", self.ezx**2 + self.ezy**2 + self.ezz**2 )




# ('ey', 0.304192241115606, 0.7845300950417162, 0.5403513768085548)
# ('norm:', 1.0)
# ('ez', 0.09865239061376316, 0.5382383859203108, -0.8369988923219004)
# ('norm:', 1.0)


        #cdef double x
 
        #x = (self.ezx**2 + self.ezy**2 + self.ezz**2 )
        
        
        #nevermind the names, just using already allocated variables
        cos2sq = 1 - (self.ezx**2 + self.ezy**2 + self.ezz**2 )
        ezz = 1 + cos2sq*.5  # 0th and 1st order
        cos2sq *= cos2sq              
        ezz += cos2sq * 0.375            # 2nd order
        cos2sq *= cos2sq
        ezz += cos2sq * 0.3125
        cos2sq *= cos2sq
        ezz += cos2sq * 0.2734375
        
        self.ezx *= ezz
        self.ezy *= ezz
        self.ezz *= ezz



        self.eyx = self.ezx
        self.eyy = self.ezy
        self.eyz = self.ezz
        
        self.ezx = tempEZx
        self.ezy = tempEZy
        self.ezz = tempEZz






    cpdef add_to_cell(self, object old_cell, object old_points,  int Npoints):
        
        
        cdef int N = self.X.size()
        cdef int i
        
        old_cell.append(N)
        
        cdef int x
        for x in range(0, N):
            old_cell.append(x + Npoints)
       # old_cell.extend(x + Npoints for x in range(0, N))
        
        cdef cnp.ndarray points = np.zeros((N, 3))
        
        
        
        for i in range(N):
            points[i, 0] = self.X[i]
            points[i, 1] = self.Y[i]
            points[i, 2] = self.Z[i]
        
        old_points.extend(points)
        #old_cell.extend(the_cell)
            
        return N
        
        
        
        
        
        
    def get_track(self):
        
        cdef int N = self.X.size()
        cdef int i
        
        cdef cnp.ndarray points = np.zeros((N, 3))
        for i in range(N):
            points[i, 0] = self.X[i]
            
            points[i, 1] = self.Y[i]
            
            points[i, 2] = self.Z[i]
        return points
            
        
        

    def track(self):
        return list(self.X), list(self.Y), list(self.Z)





    cdef void update_references(self):
        raise RuntimeError("update_references from particle.pyx was "          \
                         + "not overriden by update_references from photon.pyx")


    cdef void _run(self, mixmax_engine *genPTR):
        raise RuntimeError(">>>> _run called from particle")













# @cython.boundscheck(False)
# @cython.wraparound(False) 
# @cython.initializedcheck(False)
# @cython.cdivision(True)
# @cython.freelist(10_000_000)
# cdef class Particle:
#     """
#     Particle Class:
#         Provides methods for all geometric operations that you can do on
#         a particle.
        
#     Private Methods:
#         move(self, L):
#             moves particle in the direction of self.ez by a displacement L
        
#         change_direction(self, cos, phi):
#             changes direction of particle in local polar coordinates
#             cos = cos(theta) - theta is polar angle
#             phi              - phi is azimuthal displacement
        
#         cross(self, L):
#             takes a displacement, checks if valid intersection has ocurred
#             if so, self.current_region is updated and self.update_references 
#             is called then returns True if not, returns False
            
#         Private staticmethods:
#             INIT(*args):
#                 receives all info and initializes particle
#                 written in cython for speed, not acessible by python
#     """
    
    
    
#     cdef void INIT(Particle self,
#                   double    E, 
#                   Vector    pos,
#                   Vector    ey,
#                   Vector    ez,
#                   Volume    current_region):
#         """
#         Initialize particle without calling the python interpreter.
#         """

        
        
#         self.pos   = pos
#         self.current_region = current_region
#         self.ey = ey
#         self.ez = ez




        
#         self.E = E
#         #if DEBUG == True:
#          #   self.FILE = open("log.txt", "a")
            
    


#     cdef void _run(self):
#         pass
    
#     cdef void move(self, double L):
#         """ 
#         Move particle by displacement L in the current direction of movement.
#         """
#         #print(f"diplacement = {L}")
        
        
        
#         self.ez.fastNORMALIZE()
        
#         self.pos = Vector._new(self.pos.x + self.ez.x*L,
#                                self.pos.y + self.ez.y*L,
#                                self.pos.z + self.ez.z*L)
        
        
#         #self.pos = self.pos.ADD(self.ez.MUL(L))
#         #print(self.pos.x, self.pos.y, self.pos.z)
#         #self.record()


#     # cdef bint cross(self, double L):
#     #     """
#     #     Check for valid intersections. Change region and return true if one is
#     #     found. Else return False.
#     #     """
#     #     #print(">crossing...?")
        
#     #     #self.ex = self.ex.normalize()
#     #     #self.ey = self.ey.normalize()
#     #     #self.ez = self.ez.normalize()
        
        
        
        
#     #     self.P = self.current_region._getIntersection(self.pos, self.ez)
        
        
#     #     if 0 <= self.P.t <= L:
#     #        # print(">> yes")
#     #         self.current_region = self.P.REGION
#     #         self.move(self.P.t + SURFACE_THICKNESS)
#     #         return True
#     #    # print(">> no")
#     #     return False

#     cdef bint cross(self, double L):
        
#         if self.pos.x**2 + self.pos.y**2 + self.pos.z**2 > 50**2:
#             return True
        
        
#         return False
#         cdef Vector proposal = self.pos.ADD(self.ez.MUL(L))
        
#         cdef V = self.current_region.cross(proposal)
        
#         if self.current_region.NEW_VOL is False:
#             return False
        
        
#         return False
#         self.current_region = V
#         self.move(L + V.SDF(proposal))
        
#         return True
                
                
                
        
        


#     cdef void record(self):
#         """
#         Add current position to the track record.
#         """
        

#         self.X.push_back(self.pos.x)
#         self.Y.push_back(self.pos.y)
#         self.Z.push_back(self.pos.z)
#         self.energy.push_back(self.E)
        
#         #to_log = f"""
# #        
# #{self.pos.x} |  {self.pos.y} | {self.pos.z}
# #current_region = {self.current_region}
# #current_mat = {self.current_material}

# #"""
        
        
#         #self.FILE.write(to_log)


#     def track(self):
#         return list(self.X), list(self.Y), list(self.Z)


#     cdef double rsqrt(self, double x):
#         x = 1 - x
#         cdef double y = 1 + x*.5  # 0th and 1st order
#         x *= x              
#         y += x * 0.375            # 2nd order
#         x *= x
#         y += x * 0.3125
#         x *= x
#         y += x * 0.2734375
#         return y
        
        

#     cdef int change_direction(self, double cos, double phi):
#         """
#         Rotate the local frame of reference.
#         """
        
        
#         #norm =  self.rsqrt(self.ey.x**2 + self.ey.y**2 + self.ey.z**2)
#         norm = 1
#         #norm = 1 / (self.ey.x**2 + self.ey.y**2 + self.ey.z**2 )**.5
#         # self.ey.x = self.ey.x*norm
#         # self.ey.y = self.ey.y*norm
#         # self.ey.z = self.ey.z*norm
        
#         # #norm = 1 / (self.ez.x**2 + self.ez.y**2 + self.ez.z**2 )**.5
#         # #norm = self.rsqrt(self.ez.x**2 + self.ez.y**2 + self.ez.z**2)
#         # self.ez.x = self.ez.x*norm
#         # self.ez.y = self.ez.y*norm
#         # self.ez.z = self.ez.z*norm      
        
        
#         #self.ex = self.ex.normalize()
        
#         self.ey.fastNORMALIZE()
#         self.ez.fastNORMALIZE()
        
        
#         #if cos > 1: cos = 0.99999999
        
#         #print(f"change_direction cos = {cos}  phi = {phi}")
        
#         axis = self.ez
#         self.ey = self.ey.rotateAngle(axis, phi)
        
#         axis = self.ey
#         self.ez = self.ez.rotateCos(axis, cos)



#     cdef void update_references(self):
#         raise RuntimeError("update_references from particle.pyx was "          \
#                          + "not overriden by update_references from photon.pyx")
    


























    # cdef int propagate(self) except -1:
            
    #     cdef long double r = rand()
    #     cdef long double L = -log(r)/self.imfp_T
        
    #     if L > 1e10:
    #         raise StopSimulation
        
    #     if self.cross():
    #         return self.propagate()
        
    #     return 0
        






    # def __init__(Particle self, 
    #              Volume   space, 
    #              double   E, 
    #              Vector   pos, 
    #              double   theta, 
    #              double   phi,
    #              Vector   ex,
    #              Vector   ey,
    #              Vector   ez,
    #              bint     simulate_secondary,
    #              Volume   current_region):


    #     self.simulate_secondary = simulate_secondary
        
    #     self.space = space
    #     self.pos   = pos
    #     self.current_region = current_region

    #     self.ex, self.ey, self.ez = ex, ey, ez

    #     self.theta = theta
    #     self.phi   = phi

    #     self.change_direction(cos(self.theta), 
    #                               self.phi)
        
    #     self.E = E
    #     #self.update_ifmp() #probably should be deffered to derived init


    #     #self.children = None










































































































        # #print("entering geo")
        # cdef Point P = self.current_region._getIntersection(self.pos, self.ez)
        # #print(r, L)
       

        # if not (0 <= P.t <= L) or P.t == 1e100:
        #     self.pos = self.pos + self.ez*L
        #     return 0
        
        # # if P.t == 0:
        # #     self.current_region = P.REGION
        # #     self.update_references()
        # #     self.pos = self.pos + self.ez*1e-4
        # #     self.propagate()
        

        # self.current_region = P.REGION
        
        # #print(P.t, "reg")
        # if self.current_region.imp is 1:
        #     raise StopSimulation
        

        # self.pos = self.pos + self.ez*(P.t + 1e-4)
        # self.update_references()
        

        
        #return self.propagate()















# def rotateAngle(vec, axis, theta):
#     """ Rotates vec by an angle theta along the provided axis """
#     return Quaternion(axis = axis, angle = theta).rotate(vec)

#rotate = lambda vec, axis, theta: Quaternion(axis=axis, angle=theta).rotate(vec)


# def rotateCos(vector, axis, c):
#     """Rotate vector along axis with angle defined by c = cos(theta)"""
#     vector = Quaternion(0, *vector)
    
#     c2 = sqrt(.5*(1 + c))
#     s2 = sqrt(.5*(1 - c))
    
#     q = Quaternion(c2, *(s2*axis))
#     return (q * vector * q.conjugate).imaginary


# def choose(*args):
#     """
#     ***NOTE TO SELF: Avoid using this plz.

#     Return elements with provided probabilities.
#     """
#     n = len(args)

#     if n % 2 != 0:
#         raise Exception("Nummber of elements not equal to number of coefs.")

#     elements = [args[2*i]   for i in range(n//2)]
#     coefs    = [args[2*i+1] for i in range(n//2)]

    
#     #elements = [args[i] for i in range(0, n//2)]
#     #coefs    = [args[i] for i in range(n//2, n)]

#     coefs = array(coefs)
#     probs = coefs/sum(coefs)
#     cumul = [sum(probs[0:i]) for i in range(n//2)] + [1]
#     #print(probs)
#     r = rand()
    
#     i = searchsorted(cumul, r, side="left")
#     return elements[i - 1]    


# def dist(object P, object Q):
#     """Find distance between two points."""
#     return sum((P-Q)**2)**.5


    #maybe use this to store theta and phi
##    def change_direction(self, theta, phi):
##        '''Rotates the local frame.'''
##        axis = self.ez
##        self.ey = rotate(self.ey, axis, phi)
##        #self.ex = rotate(self.ex, axis, phi)
##
##        axis = self.ey
##        self.ez = rotate(self.ez, axis, theta)
##        #self.ex = rotate(self.ex, axis, theta







# class Particle(CParticle):
#     """Python wrapper.
#     Cython does not allow keywords.
#     Might need to change this behaviour later though, I need the inits
#     to be as quick as possible.


#     Nevermind, I'm already in the process of changing this stuff.
#     Keyword args should be reserved to sources.


#     Maybe keep this for a user front end, idk"""
                
    
#     def __init__(self,
#                  space  = [],
#                  E      = 6.,
#                  pos    = [0., 0., 0.],
#                  theta  = 0.,
#                  phi    = 0.,
#                  ex = [1., 0., 0.], ey = [0., 1., 0.], ez = [0., 0., 1.],
#                  current_region = None,
#                  simulate_secondary = False):
        

#         pos = Vector(pos[0], pos[1], pos[2])
#         ex = Vector(ex[0], ex[1], ex[2])
#         ey = Vector(ey[0], ey[1], ey[2])
#         ez = Vector(ez[0], ez[1], ez[2])

        

#         super().__init__(space, E, pos, 
#                          theta, phi, 
#                          ex, ey, ez, 
#                          simulate_secondary, current_region)



#     def summary(self):
#         return {'pos': self.pos,
#                 'E'  : self.E,
#                 'ez' : self.ez}


    
#     def plot(self, tube_radius = 1, line_width = 1):
#         """Plot the path of the particle. (uses mayavi)"""
#         history = array(self.pos_)
#         X, Y, Z = history[:,0], history[:,1], history[:,2]
#         plot3d(X, Y, Z, tube_radius = tube_radius, line_width = line_width)
#         show()



#     def _add_plot(self, fig):
#         """
#         Don't use this. Just an attempt of using other modules for visualizing.
#         Keeping it here to remember later that there are other options.
#         (there is also paraview, but it's very time consuming to learn that)s
#         """
#         from vispy import scene, visuals
#         Plot3D = scene.visuals.create_visual_node(visuals.LinePlotVisual)


#         history = array(self.pos_)
#         X,Y,Z = history[:,0], history[:,1], history[:,2]
#         pos = c_[X, Y, Z]
        
#         Plot3D(pos, width = 2.0, color='red', edge_color='w', symbol='o',
#               face_color = (0.2, 0.2, 1, 0.8),
#               parent     = fig)
        
#     def add_plot(self, fig, plot_secondary = False):
#         """
#         Add trajectory to a mayavi figure. Used by a source instance to plot the
#         results of the simulation.
#         """
        
#         name = self.__class__.__name__
#         if name == 'Photon': color = (0, 0, 1)
#         else:                color = (1, 0, 0)
        
#         history = array(self.pos_)
#         X,Y,Z   = history[:,0], history[:,1], history[:,2]
#         plot3d(X, Y, Z, figure=fig, tube_radius=.8, color=color)
        
#         if plot_secondary is True:
#             for p in self.children_[1:]:
#                 p.add_plot(fig, plot_secondary = True)