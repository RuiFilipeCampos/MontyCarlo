from distutils.core import setup, Extension
from Cython.Build import cythonize


ext = Extension("pyRelax", sources = ["pyRelax.pyx"])


# you specify the c source file in the sources list
#ext = Extension('_cppRelaxAPI', sources = ['_cppRelaxAPI.pyx'])
#setup(name="C spam", ext_modules = cythonize([ext]))
print(ext)

#from setuptools import setup
#from Cython.Build import cythonize


a = cythonize([ext])

print(a)
setup(ext_modules=a)