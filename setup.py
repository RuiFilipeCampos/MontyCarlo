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
    from setuptools import setup, find_packages
    from setuptools import Extension
except ImportError:
    from distutils.core import setup, find_packages
    from distutils.extension import Extension

from Cython.Distutils import build_ext
import numpy as np






args = ["-O2", "-fp:fast"]
ext_modules = [
Extension("tools.*",              ["MontyCarlo\\tools\\*.pyx"],               extra_compile_args = args),
Extension("particles.*",          ["MontyCarlo\\particles\\*.pyx"],           extra_compile_args = args),
Extension("*",                    ["MontyCarlo\\*.pyx"],                      extra_compile_args = args), 
Extension("geometry.*",           ["MontyCarlo\\geometry\\*.pyx"],            extra_compile_args = args),
Extension("materials.electron.*", ["MontyCarlo\\materials\\electron\\*.pyx"], extra_compile_args = args),
Extension("materials.positron.*", ["MontyCarlo\\materials\\positron\\*.pyx"], extra_compile_args = args),
Extension("materials.*",          ["MontyCarlo\\materials\\*.pyx"],           extra_compile_args = args),
Extension("materials.photon.*",   ["MontyCarlo\\materials\\photon\\*.pyx"],   extra_compile_args = args),
Extension("_random.*",             ["MontyCarlo\\_random\\*.pyx"],              extra_compile_args = args)
]
#-ffast-math

#,Extension("geometry.CGS*", ["geometry\\CGS\\*.pyx"])]



#ext_modules = [Extension("my_code_cython",["my_code_cython.pyx"]),
#               Extension("another_code_cython",["another_code_cython.pyx"])]

import os


# COMPILER OPTIONS
import Cython.Compiler.Options


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


#os.environ['CFLAGS'] = '-O3 -Wall -std=c++11 -I"some/custom/paths"'
setup(
    name = "MontyCarlo", version = 0.0.34, author = "Rui Filipe de Sousa Campos",
    description = "A fast general purpose monte carlo particle simulator (photons, electrons and positrons). Written in Cython, Python and C++.",
    long_description = long_description, long_description_content_type="text/markdown",     url="https://github.com/RuiFilipeCampos/MontyCarlo",
    setup_requires=['setuptools_scm'],
    install_requires=['requests', 
                      'gdown', 
                      'numpy', 
                      'scipy', 
                      'matplotlib', 
                      'pyvista', 
                      'numba', 
                      'bs4', 
                      'pandas', 
                      'plotly', 
                      'scikit-image',
                      'Jinja2', 
                      'pyunpack',
                      'patool'],
    include_package_data=True,
    packages=find_packages(),
    cmdclass = {'build_ext': build_ext},
    include_dirs = [".", np.get_include(), "_random"], #, "MontyCarlo\\materials"
    ext_modules = cythonize(ext_modules, 
    						annotate=False,
                            compiler_directives={'profile': False, 'language_level' : "3"})) 






#	print("Downloading positron elastic data...")
#	#from .materials.positron import download_database






# https://journals.sagepub.com/doi/suppl/10.1177/ANIB_39_2
#  https://iopscience.iop.org/article/10.1088/1742-6596/1662/1/012021/pdf
# https://iopscience.iop.org/article/10.1088/0031-9155/51/14/017

 
