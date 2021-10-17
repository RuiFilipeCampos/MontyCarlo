


# Imports will be cleaner eventually...
from MontyCarlo import *
from MontyCarlo.sources import *
from MontyCarlo.geometry.CSG import *


# Define your materials. They are compiled once per project, so don't worry if it takes a while
# the first time you run the script, the second time `Mat` will just load from cache.

water = Mat({1:2, 8:1}, 1, name = "Water")

air = Mat({6:1.50187000E-04,
           7:7.84430000E-01,
           8:2.10748000E-01,
           18:4.67111000E-03}, 
           1.20479e-03,
           C1 = .2, 
           C2 = .2,
           name = "AirDryNearSeaLevel")

gold = Mat({79:1}, 1.93200000E+01, name = "Gold")


# Define geometry
# The indentation tells you exactly how the BVH is being constructed.
# The `.configure` method will be replaced with a smarter system.
with InfiniteVolume() as outer:
    outer.configure("OUTER", render = False)
    outer.fill(gold) 
    
    with Sphere(100) as outer_sphere:
        outer_sphere in outer
        outer_sphere.configure("outer_sphere", render = True)
        outer_sphere.fill(air)

        with Sphere(50) as inner_sphere:
            inner_sphere in outer_sphere
            inner_sphere.configure("inner_sphere", render = True)
            inner_sphere.fill(water)

            with Sphere(25) as inner_sphere2:
                inner_sphere2 in inner_sphere
                inner_sphere2.configure("inner_sphere2", render = True)
                inner_sphere2.fill(air)



photon_beam = Beam(
                   "photon",       # kind of particle 
                   inner_sphere2,   # initial volume
                   E = 10e6,       # initial eneryg in eV
                   N = 1_000,     # number of particles in the source, careful with this number, might break your run and fill your ram
                   pos = (0, 0, 0) # initial position 
                  ) 



# let Plotter handle the run
plotter = Plotter(photon_beam)


# then ask it for a fig
fig = plotter.new_plot()

# use this method to draw the geometry onto the figure (this will be better)
plotter.add_geometry(fig, outer)


fig.show()
