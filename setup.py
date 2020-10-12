import os

old_file_name = "_geometry.pyx"
new_file_name = "geometry.pyx"
os.rename(old_file_name, new_file_name)


try:
	import setuptools  # important
	from distutils.core import setup
	from Cython.Build import cythonize

	setup(ext_modules=cythonize("geometry.pyx"))
except: pass

os.rename(new_file_name, old_file_name)

