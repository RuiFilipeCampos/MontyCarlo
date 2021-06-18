# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 12:41:36 2020

@author: Rui Campos
"""

# from setuptools import setup
from Cython.Build import cythonize

# setup(
#     ext_modules = cythonize("tools/*.pyx")
# )



try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Distutils import build_ext
import numpy as np


ext_modules = [Extension("tools.*", ["MontyCarlo\\tools\\*.pyx"]),
               Extension("particles.*", ["MontyCarlo\\particles\\*.pyx"]),
               Extension("particles.*", ["MontyCarlo\\particles\\*.pyx"])]
#ext_modules = [Extension("my_code_cython",["my_code_cython.pyx"]),
#               Extension("another_code_cython",["another_code_cython.pyx"])]

setup(
    name= 'Generic model class',
    cmdclass = {'build_ext': build_ext},
    include_dirs = [np.get_include()],
    ext_modules = cythonize(ext_modules))