# cython: profile=True
from numpy import *
import sys

def plotAxes(rangeX, rangeY, rangeZ, fig, tr = .3):
    """Plot axes on mayavi.mlab figure."""
    plot3d    = sys.modules['mayavi'].mlab.plot3d
    a = zeros(2)
    
    
    plot3d(rangeX, a, a, figure=fig, tube_radius=tr)
    plot3d(a, rangeY, a, figure=fig, tube_radius=tr)
    plot3d(a, a, rangeZ, figure=fig, tube_radius=tr)




cdef class Vol:
    """Base class for volumes."""
    cdef public int imp, Nin, Nout, N0
    cdef public double E, E0
    cdef public object material, condition
    cdef public list tally_out

    
    
    def __init__(self, object condition):
        self.imp  = 1 
        self.Nin = 0
        self.Nout = 0
        self.E0 = 0.
        self.N0 = 0
        self.E = 0.
        self.tally_out = []
        self.material  = None
        self.condition = condition


    def scatter_tally(self, fig, size = 1):
        pos = [particle.pos for particle in self.tally]
        pos = array(pos)
        
        X, Y, Z = pos[:, 0], pos[:, 1], pos[:, 2]
            
        points3d(X, Y, Z, figure = fig, scale_factor = size)
    def fill(self, material):
        """Saves Material object in Vol object."""
        self.material = material
    def __str__(self):
        try:
            self.per_track = self.E/self.Nin
        except ZeroDivisionError:
            self.per_track = None
        __info__ = f"""
        -------------------------------------
        TALLY
        -------------------------------------
        Importance: {self.imp}

        Initial number of particles: {self.N0}
        Total initial energy:        {self.E0}
        
        Energy Absorved: {self.E} MeV
        Particles in:    {self.Nin}
        Particles out:   {self.Nout}

        Energy deposited per particle that entered the volume:
        {self.per_track} MeV/particle
        -------------------------------------
        """
        return __info__


    def getTally(self):
        info = {'N0':self.N0,
                'E0':self.E0,
                'Eabs':self.E,
                'Nin':self.Nin,
                'Nout':self.Nout,
                'Eabs per track': self.E/self.Nin}
        return info


    def __getitem__(self, subs):
        
        x0, xf = subs[0].start, subs[0].stop
        y0, yf = subs[1].start, subs[1].stop
        z0, zf = subs[2].start, subs[2].stop

        if x0 is not None and xf is not None:
            def X(x,y,z):
                f = self.condition
                g = lambda x,y,z: logical_and(x>x0, x<xf)
                return logical_and(f(x,y,z), g(x,y,z))
        else: X = self.condition 

        if y0 is not None and yf is not None:
            def Y(x,y,z):
                f = X
                g = lambda x,y,z: logical_and(y>y0, y<yf)
                return logical_and(f(x,y,z), g(x,y,z))
        else: Y = X
            
        if z0 is not None and zf is not None:
            def Z(x,y,z):
                f = Y
                g = lambda x,y,z: logical_and(z>z0, z<zf)
                return logical_and(f(x,y,z), g(x,y,z))
        else: Z = Y

        return Vol(Z)
 

        

        

        
    #LOGIC
    def __or__(self, other):
        """Union of two regions. Example: new_region = region1 | region2"""
        def new_condition(x,y,z):
            return logical_or(self.condition(x, y, z), other.condition(x, y, z))
        return Vol(new_condition)
    
    def __and__(self, other):
        """Intersection of two regions. Example: new_region = region1 & region2"""
        def new_condition(x, y, z):
            return logical_and(self.condition(x,y,z), other.condition(x,y,z))
        return Vol(new_condition)
    
    def __invert__(self):
        """Negation of a region. Example: inverted_region = ~region"""

        def new_condition(double x, double y, double z):
            return logical_not(self.condition(x, y, z))

        return Vol(new_condition)
        
    def __sub__(self, other):
        """
        Subtraction of two regions:
            new_region = region2 - region1
        region1 has been striped of all the points that it had in common with region1
        """
        return (~other) & self
    
    def __contains__(self, object pos):
        """Checks if position is in this region. Example:
        > (x, y, z) in region
        returns True or False.
        """
        cdef float x = pos[0]
        cdef float y = pos[1]
        cdef float z = pos[2]

        return self.condition(x, y, z)



  
    #MOVE VOLUME
    def moveX(self, x0):
        """
        Move region along x-axis. Example:
        > region.moveX(10)
        Returns new volume region.
        """
        def new_condition(x,y,z):
            return self.condition(x-x0, y, z)
        return Vol(new_condition)
    
    def moveY(self, y0):
        """
        Move region along x-axis. Example:
        > region.moveX(10)
        Returns new volume region.
        """
        def new_condition(x,y,z):
            return self.condition(x, y-y0, z)
        return Vol(new_condition)
    
    def moveZ(self, z0):
        """
        Move region along x-axis. Example:
        > region.moveX(10)
        Returns new volume region.
        """
        def new_condition(x,y,z):
            return self.condition(x, y, z-z0)
        return Vol(new_condition)
    
    def __iadd__(self, displacement):
        """
        Move every point in volume by displacement vector.
        EX:
        sphere = Sphere(10) #sphere radius 10 at (0,0,0)
        sphere += [1,2,1]   #center of sphere is now at (1,2,1)
        """
        
        x0, y0, z0 = displacement
        def new_condition(x, y, z):
            return self.condition(x-x0,y-y0,z-z0)

        return Vol(new_condition)
    
    
        
    #SHAPE VOLUME
    def limit(self, axis, k0, kf):
        """
        Limits region along the following axis:
            axis = "x"
            axis = "y"
            axis = "z"

        Example:
        > region.limit("x", 1, 10)
        region lost all its points outside of 1 < x < 10.
        """
        if axis == "x":
            def new_condition(x,y,z):
                f = self.condition
                g = lambda x,y,z: logical_and(x>k0, x<kf)
                return logical_and(f(x,y,z), g(x,y,z))
        if axis == "y":
            def new_condition(x,y,z):
                f = self.condition
                g = lambda x,y,z: logical_and(y>k0, y<kf)
                return logical_and(f(x,y,z), g(x,y,z))
        if axis == "z":
            def new_condition(x,y,z):
                f = self.condition
                g = lambda x,y,z: logical_and(z>k0, z<kf)
                return logical_and(f(x,y,z), g(x,y,z))
        return Vol(new_condition)
        
        
    
    def plotVol(self, *args,
                rangeT  = [-500, 500],
                rangeX  = None,
                rangeY  = None,
                rangeZ  = None,
                N       = 100j,
                Nx      = None,
                Ny      = None,
                Nz      = None,
                opacity = 0.2,
                color   = (0, 0, 1),
                fig     = None,
                Id      = None,): #to do

        """
        A mayavi figure must be provided.


        ** rangeX, rangeY, rangeZ **
        range of plot in each direction
        if no value is provided,
        it is substituted by rangeT  d
        efault value of rangeT is [-50, 50]

        ** Nx, Ny, Nz **
        resolution of plot in each direction
        if no value is provided, it is substituted by N
        default value of N is 100j

        fig   -> if a mayavi.mlab scene is provided, it plots the volume there
        color -> RGB tuple ex: (1, 0, 0) for Red
        """
        
        contour3d = sys.modules['mayavi'].mlab.contour3d
        
        if fig is None:                           raise ValueError("You must provide a figure. ('fig = figura')")
        if fig.__class__.__name__ is not 'Scene': raise ValueError('The fig keyword argument must be a mayavi figure!')
        
        if rangeX is None: rangeX = rangeT
        if rangeY is None: rangeY = rangeT
        if rangeZ is None: rangeZ = rangeT

        if args:
            x1, x2 = args[0], args[1]
            y1, y2 = args[2], args[3]
            z1, z2 = args[4], args[5]
        else:
            x1, x2 = rangeX
            y1, y2 = rangeY
            z1, z2 = rangeZ

        if Nx is None: Nx = N
        if Ny is None: Ny = N
        if Nz is None: Nz = N
        
        grid = mgrid[x1:x2:Nx,
                    y1:y2:Ny,
                    z1:z2:Nz]

        X, Y, Z = grid
        
        _, Nx, Ny, Nz = grid.shape

        scalar = empty((Nx, Ny, Nz))

        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    x = X[i][j][k]
                    y = Y[i][j][k]
                    z = Z[i][j][k]

                    scalar[i][j][k] = int(Vector(x, y, z) in region)


        

        plotAxes(rangeX, rangeY, rangeZ, fig)
        
        contour3d(X, Y, Z,
                  scalar,
                  figure  = fig,
                  opacity = opacity,
                  color   = color)
        return fig


class Cube(Vol):
    def __init__(self, x0, y0, z0):
        X  = lambda x,y,z: logical_and(0<x, x<x0)
        Y  = lambda x,y,z: logical_and(0<y, y<y0)
        Z  = lambda x,y,z: logical_and(0<z, z<z0)
        XY = lambda x,y,z: logical_and(X(x,y,z), Y(x,y,z))
        super().__init__(lambda x, y, z: logical_and(XY(x, y, z), Z(x, y, z)))
        

class Sphere(Vol):
    """Create sphere of radius 'radius' at position (0, 0, 0)."""

    def __init__(self, double radius):


        def condition(double x, double y, double z):
            return x**2 + y**2 + z**2 < radius**2

        super().__init__(condition)

        
        

class CylinderZ(Vol):
    """Create infinite cylinder whose axis is the z-axis."""
    def __init__(self, radius):
        self.condition = lambda x, y, z: x**2+y**2< radius**2
        
class CylinderX(Vol):
    """Create infinite cylinder whose axis is the x-axis."""
    def __init__(self, radius):
        self.condition = lambda x,y,z: z**2+y**2 < radius**2
        
class CylinderY(Vol):
    """Create infinite cylinder whose axis is the y-axis."""
    def __init__(self, radius):
        self.condition = lambda x,y,z: x**2+z**2 < radius**2

class Ellipsoid(Vol):
    def __init__(self, ax, bx, cx):
        """
        Creates an ellipsoid whose center point is (0, 0, 0). Equation of volume:
        
        (x/a)**2 + (y/b)**2 + (z/c)**2<1

        """
        
        self.condition = lambda x,y,z: (x/ax)**2 + (y/bx)**2 + (z/cx)**2 < 1

class CylElipZ(Vol):
    """
    Create an infinite cylinder whose basis is an ellipse.
    """
    def __init__(self, ax, by):
        self.condition = lambda x,y,z: (x/ax)**2+(y/by)**2 < 1
        
class CylElipX(Vol):
    def __init__(self, az, by):
        self.condition = lambda x,y,z: (z/az)**2+(y/by)**2 < 1
        
class CylElipY(Vol):
    def __init__(self, ax, bz):
        self.condition = lambda x,y,z: (x/ax)**2+(z/bz)**2 < 1

class pX(Vol):
    def __init__(self, X, side):
        """ Create a volume with surface defined by x=X.
        side = "+"
        or
        side = "-"
        """
        if side == "+":
            def condition(x, y, z):
                return x>X
            self.condition = condition
        else:
            def condition(x, y, z):
                return x<X
            self.condition = condition

class pY(Vol):
    def __init__(self, Y, side):
        """ Create a volume with surface defined by y=Y.
        side = "+"
        or
        side = "-"
        """
        if side == "+":
            def condition(x, y, z):
                return y>Y
            self.condition = condition
        else:
            def condition(x, y, z):
                return y<Y
            self.condition = condition

class pZ(Vol):
    def __init__(self, Z, side):
        """ Create a volume with surface defined by z=Z.
        side = "+"
        or
        side = "-"
        """
        if side == "+":
            def condition(x, y, z):
                return z>Z
            self.condition = condition
        else:
            def condition(x, y, z):
                return z<Z
            self.condition = condition




print("> Imported geometry!")
