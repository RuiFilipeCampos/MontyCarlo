__doc__ = """
plotter.py
    Where plots are made. xd
"""

print("Importing `.plotter`")								 
										 
		

from collections import deque


import numpy as np
import pyvista as pv

class Plotter:


	def generate(self, par_data, mfp = False):
		electron_cell = deque()
		electron_points = deque()
		electron_Npoints = 0


		
		for record in par_data:
			N = len(record)

			points = np.zeros((N, 3))

			for i in range(N):
				x, y, z = record[i]['x'], record[i]['y'], record[i]['z']


				

				if x**2 + y**2 + z**2 > 1000**2:
					N = i + 1
					break

				points[i, 0] = x
				points[i, 1] = y
				points[i, 2] = z

			points = points[:N, :]
			electron_cell.append(N)
			electron_cell.extend(x + electron_Npoints for x in range(0, N))
			electron_points.extend(points)
			electron_Npoints += N



		return electron_cell, electron_points, electron_Npoints



	def __init__(self, source):
		source.run(record_pos = True)
		self.photon_mesh = False
		self.electron_mesh = False
		self.positron_mesh = False

		
		data = source.pos_record
		if len(data[0]) != 0:
			photons = self.generate(data[0], mfp = True)
			self.photon_mesh = pv.PolyData(np.array(photons[1]),  lines = np.array(photons[0]), n_lines = photons[2])
		
		if len(data[1]) != 0:
			el = self.generate(data[1])
			self.electron_mesh = pv.PolyData(np.array(el[1]),  lines = np.array(el[0]), n_lines = el[2])
		
		if len(data[2]) != 0:
			po = self.generate(data[2])
			self.positron_mesh = pv.PolyData(np.array(po[1]),  lines = np.array(po[0]), n_lines = po[2])
		
	@staticmethod
	def add_geometry(plt, upper):

		for vol in upper:
			plt.add_mesh(vol.get_mesh(), opacity = 0.1, color = [1, 1, 1])
			Plotter.add_geometry(plt, vol)

		#return plotter

	def new_plot(self, ph_opacity = .1,el_opacity = 1, po_opacity = 1):    
		import pyvista as pv
		pv.set_plot_theme("night")
		plotter = pv.Plotter()
		
		
		if self.photon_mesh:
			plotter.add_mesh(self.photon_mesh, line_width = 0.001,    color='white', opacity=ph_opacity)
		
		if self.electron_mesh:
			plotter.add_mesh(self.electron_mesh, line_width = 0.001,  color='blue', opacity=el_opacity)
		if self.positron_mesh:
			plotter.add_mesh(self.positron_mesh,line_width = 0.001,   color='red', opacity=po_opacity)
				
		return plotter




		
