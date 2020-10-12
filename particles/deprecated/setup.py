import os

old_file_name = "_particle.pyx"
new_file_name = "particle.pyx"
os.rename(old_file_name, new_file_name)

old_file_name = "_photons.pyx"
new_file_name = "photons.pyx"
os.rename(old_file_name, new_file_name)


try:
	import setuptools  # important
	from distutils.core import setup
	from Cython.Build import cythonize

	import numpy
	setup( ext_modules  = cythonize(["*.pyx"]    ),
		   include_dirs = [numpy.get_include()]  )
except: pass


old_file_name = "_particle.pyx"
new_file_name = "particle.pyx"
os.rename(new_file_name, old_file_name)

old_file_name = "_photons.pyx"
new_file_name = "photons.pyx"
os.rename(new_file_name, old_file_name)

