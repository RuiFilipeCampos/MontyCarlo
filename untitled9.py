# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 10:34:26 2021

@author: Rui Campos
"""


from __future__ import absolute_import, division, print_function
from mayavi import mlab
import numpy as np
import math






from MontyCarlo.geometry.CSG import Sphere
from MontyCarlo.materials.materials import Material




# gold = Material({79:1}, 19.3,
#                 name = "GOLD" , 
#                 Wcc = 1e3, 
#                 Wcr = 1e6, 
#                 C1 = 0.01, 
#                 C2 = 0.01)


args = ({1:2, 8:1}, 1)
params = dict(Wcc = 2000, Wcr = 1e3, C1 = 0.01, C2 = 0.01)
others = dict(name = "WATER")
water = Material(*args, **params, **others)


from MontyCarlo.particles.electrons import Electron


r1 = 5000
#r2 = 50
#r3 = 60
S1 = Sphere(r1)
S1.fill(water)



ee = Electron.new(S1, S1, E = 50e3)
ee.run()
xs, ys, zs = ee.getTrack()


mlab.figure(1, size=(400, 400), bgcolor=(0, 0, 0))
#mlab.points3d(0,0,0,  scale_factor=.025)
plt = mlab.plot3d(xs, ys, zs,  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))
plt2 = mlab.plot3d(xs, ys, zs,  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))
plt3 = mlab.plot3d(xs, ys, zs,  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))
plt4 = mlab.plot3d(xs, ys, zs,  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))
plt5 = mlab.plot3d(xs, ys, zs,  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))


ms = plt.mlab_source
ms2 = plt2.mlab_source
ms3 = plt3.mlab_source
ms4 = plt4.mlab_source
ms5 = plt5.mlab_source

@mlab.animate(delay=10)
def anim():
    #f = mlab.gcf()
    while True:
        ee = Electron.new(S1, S1, E = 50e3)
        ee.run()
        xs, ys, zs = ee.getTrack()
        ms.reset(x=xs, y=ys, z=zs,  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))
        
        
        ee = Electron.new(S1, S1, E = 50e3)
        ee.run()
        xs, ys, zs = ee.getTrack()
        ms2.reset(x=xs, y=ys, z=zs,  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))
        
        
        ee = Electron.new(S1, S1, E = 50e3)
        ee.run()
        xs, ys, zs = ee.getTrack()
        ms3.reset(x=xs, y=ys, z=zs,  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))
        
        ee = Electron.new(S1, S1, E = 50e3)
        ee.run()
        xs, ys, zs = ee.getTrack()
        ms4.reset(x=xs, y=ys, z=zs,  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))        
        
        ee = Electron.new(S1, S1, E = 50e3)
        ee.run()
        xs, ys, zs = ee.getTrack()
        ms5.reset(x=xs, y=ys, z=zs,  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))
        
        
        yield
        
        
        
        # #for (x, y, z) in zip(xs, ys, zs):
        # for i in range(len(xs)-3):
        #       #print('Updating scene...')
        #       ms.reset(x=xs, y=ys, z=zs,  tube_radius = None)
        #       #plt.mlab_source.set(x=xs[:i+2], y=ys[:i+2], z=zs[:i+2])
        #       #f.scene.reset_zoom()
        #       yield


# anim()
# mlab.show()




def f():
    
    ee = Electron.new(S1, S1, E = 50e3)
    ee.run()
    xs, ys, zs = ee.getTrack()   
    xs = np.array(xs)
    ys = np.array(ys)
    zs= np.array(zs)
    mlab.figure(1, size=(400, 400), bgcolor=(0, 0, 0))
    #mlab.points3d(0,0,0,  scale_factor=.025)
    plt = mlab.plot3d(xs[:2], ys[:2], zs[:2],  tube_radius = None, line_width = 1,opacity = 0.4, color=(1, 0, 0))
    ms = plt.mlab_source
    return xs, ys, zs, ms

import time

@mlab.animate(delay=30)
def anim():
    
    while True:
        #f = mlab.gcf()
        xs, ys, zs, ms = f()
        xs1, ys1, zs1, ms1 = f()
        xs2, ys2, zs2, ms2 = f()
        xs3, ys3, zs3, ms3 = f()
        xs4, ys4, zs4, ms4 = f()
        
        print(len(xs))
    
        for i in range(3, len(xs), 5):
            print(i)
    
            ms.reset(x=xs[:i], 
                     y=ys[:i], 
                     z=zs[:i],  tube_radius = None, line_width = 1,opacity = 0.4,  color=(1, 0, 0))
            ms1.reset(x=xs1[:i], 
                      y=ys1[:i], 
                      z=zs1[:i],  tube_radius = None, line_width = 1,opacity = 0.4,  color=(1, 0, 0))
            ms2.reset(x=xs2[:i], 
                      y=ys2[:i], 
                      z=zs2[:i],  tube_radius = None, line_width = 1,opacity = 0.4,  color=(1, 0, 0))
            
            ms3.reset(x=xs3[:i], y=ys3[:i], z=zs3[:i],  tube_radius = None, line_width = 1,opacity = 0.4,  color=(1, 0, 0))
            ms4.reset(x=xs4[:i], y=ys4[:i], z=zs4[:i],  tube_radius = None, line_width = 1,opacity = 0.4,  color=(1, 0, 0))
            
            yield
        
        
        
        # #for (x, y, z) in zip(xs, ys, zs):
        # for i in range(len(xs)-3):
        #       #print('Updating scene...')
        #       ms.reset(x=xs, y=ys, z=zs,  tube_radius = None)
        #       #plt.mlab_source.set(x=xs[:i+2], y=ys[:i+2], z=zs[:i+2])
        #       #f.scene.reset_zoom()
        #       yield


anim()
mlab.show()








# %%
def plot(tracks):
    """
    Parameters.
    ----------
    Ntracks : int
        NUMBER OF TRACKS TO BE PLOTTED.
    view : tuple
        SET THE CAMERA VIEW.
    E : int or float, optional
        INITIAL ENERGY OF PARTICLE. The default is 6.

    Returns
    -------
    None.
    """
    
    X, Y, Z = tracks[0]
    N = len(X)
    for i in range(N)
        mlab.clf()
    
        # _test_mesh(r1)
        # _test_mesh(r2)
        # _test_mesh(r3)
    
        index = 0
        x = list()
        y = list()
        z = list()
        connections = list()
        s = list()
    
        for track in tracks:
            X, Y, Z = track
    
            # removing last 20 points
            # geometry bug makes some trajectories go to infinity
            X = X[i:-20]
            Y = Y[i:-20]
            Z = Z[i:-20]
            # S = ee.getEnergy()
    
            N = len(X)
            S = [10]*N
            S.append(10)
    
            x.append(np.array(X))
            y.append(np.array(Y))
            z.append(np.array(Z))
            s.append(np.array(S))
    
            connections.append(np.vstack(
                                          [np.arange(index,     index + N - 1.5),
                                           np.arange(index + 1, index + N - .5)]
                                    ).T)
    
            index += N
    
        # print("FINISHED SIMULATION")
    
        # Now collapse all positions, scalars and connections in big arrays
        x = np.hstack(x)
        y = np.hstack(y)
        z = np.hstack(z)
        s = np.hstack(s)
        connections = np.vstack(connections)
    
        # Create the points
        src = mlab.pipeline.scalar_scatter(x, y, z)
    
        # Connect them
        src.mlab_source.dataset.lines = connections
        src.update()
        mlab.pipeline.surface(src, color=(1, 0, 0), line_width=1, 
                                                    opacity=opacity, **kwargs)
        
        mlab.view(90, 90, distance=.1e-3, focalpoint=(0, 0, 25e-6))
        yield
    
        #mlab.view(*view, distance=0.5e-3, focalpoint=(0, 0, 0))
    
        # print("FINISHED PLOT")
    
        # mlab.show()
    














    


