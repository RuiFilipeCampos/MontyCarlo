#Internal Imports
from .primitives import *


#External Imports

from numpy import *
from mayavi.mlab import *
import mayavi.mlab as mlab
import sys
from ..tools.vectors import *


def plotAxes(rangeX, rangeY, rangeZ, fig, tr = .3):
	"""Plot axes on mayavi.mlab figure."""
	plot3d    = sys.modules['mayavi'].mlab.plot3d
	a = zeros(2)
	
	
	plot3d(rangeX, a, a, figure=fig, tube_radius=tr)
	plot3d(a, rangeY, a, figure=fig, tube_radius=tr)
	plot3d(a, a, rangeZ, figure=fig, tube_radius=tr)

def plotVol(region, *args,
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

	contour3d = sys.modules['mayavi'].mlab.contour3d

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

	if not fig: fig = mlab.figure()
	plotAxes(rangeX, rangeY, rangeZ, fig)
		

	contour3d(X, Y, Z,
			scalar,
			figure  = fig,
			opacity = opacity,
			color   = color)
	return fig






