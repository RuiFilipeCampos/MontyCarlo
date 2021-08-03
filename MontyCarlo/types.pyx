__doc__ = """
types.pyx
    A place to keep the frequently used data structures. (I will eventually move everything here).
"""

__author__ = "Rui Campos"






print("Importing `.types`")

import numpy as np
nan = np.nan

cdef struct double3:
    double x, y, z

class MAP(dict):
    """A dictionary with the following syntax sugar:
    
    > map = MAP()
    > map.x = 1
    
    is equivelent to:
    
    > map = MAP()
    > map['x'] = 1
    
    """
    
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError

cdef class py_state:
    """The python analog of `STATE`:

    ```
    cdef struct STATE:
        mixmax_engine* genPTR
        void *current_region
        double3 pos
        double3 dire
        double3 axis
        double E
        double L 
        double last_displacement
    ```

    """

    def __init__(self, 
    pos  = np.array([0, 0, 0], dtype = float),
    dir  = np.array([0, 0, 1], dtype = float),
    axis = np.array([0, 1, 0], dtype = float), 
    E = 1e6, L = 0, last_displacement = 0
    ):

        self.pos = pos
        self.dir = dir
        self.axis = axis
        self.E = E
        self.L = L
        self.last_displacement = last_displacement


    def to_cython(self):

        cdef STATE state

        self.state.pos.x = self.pos[0]
        self.state.pos.y = self.pos[1]
        self.state.pos.z = self.pos[2]

        self.state.dir.x = self.dir[0]
        self.state.dir.y = self.dir[1]
        self.state.dir.z = self.dir[2]

        self.state.axis.x = self.axis[0]
        self.state.axis.y = self.axis[1]
        self.state.axis.z = self.axis[2]

        self.state.E = self.E
        self.state.L = self.L
        self.state.last_displacement = self.last_displacement

        return self.state