
from .. import particle as pa

class RandomWalker(pa.Particle):
    def __init__(self, space, current_region
                 E     = 6,
                 pos   = Vector(0., 0., 0.),
                 theta = 0,
                 phi   = 0,
                 ex    = Vector(1., 0., 0.), 
                 ey    = Vector(0., 1., 0.), 
                 ez    = Vector(0., 0., 1.),
                 simulate_secondary = False):

        super().__init__(space, E, pos, 
                         theta, phi, 
                         ex, ey, ez, 
                         simulate_secondary,
                         current_region)