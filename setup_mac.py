# -*- coding: utf-8 -*-

import os
from Cython.Build import cythonize

try:
    from setuptools import setup, find_packages
    from setuptools import Extension
except ImportError:
    from distutils.core import setup, find_packages
    from distutils.extension import Extension

from Cython.Distutils import build_ext
import Cython.Compiler.Options # COMPILER OPTIONS
import numpy as np # need to compile it with the extension modules



# MSVC ARGUMENTS
args = [
        "-std=c++11"
       ]

ext_modules = [ 
                Extension("tools.*",               ["MontyCarlo/tools/*.pyx"] ),                extra_compile_args = args),
                Extension("particles.*",           ["MontyCarlo/particles/*.pyx"]),           extra_compile_args = args),
                Extension("*",                     ["MontyCarlo/*.pyx"]),                       extra_compile_args = args), 
                Extension("geometry.*",            ["MontyCarlo/geometry/*.pyx"]),          extra_compile_args = args),
                Extension("materials.electron.*",  ["MontyCarlo/materials/electron/*.pyx"]),extra_compile_args = args),
                Extension("materials.positron.*",  ["MontyCarlo/materials/positron/*.pyx"]), extra_compile_args = args),
                Extension("materials.*",           ["MontyCarlo/materials/*.pyx"]),           extra_compile_args = args),
                Extension("materials.photon.*",    ["MontyCarlo/materials/photon/*.pyx"]),   extra_compile_args = args),
                Extension("_random.*",             ["MontyCarlo/_random/*.pyx"],             extra_compile_args = args)
              ]
 


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

if __name__ == "__main__":
    setup(
        name = "MontyCarlo",
        version = "0.0.35",
        author = "Rui Filipe de Sousa Campos",
        description = "A fast general purpose monte carlo particle simulator (photons, electrons and positrons). Written in Cython, Python and C++.",
        long_description = long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/RuiFilipeCampos/MontyCarlo",


        setup_requires   = ['setuptools_scm'],
        install_requires = ['requests',     # for downloading databases
                            'gdown',        # for downloading databases
                            'numpy',        # for data processing and C-API
                            'scipy',        # for data processing
                            'matplotlib',   # for data visualization
                            'pyvista',      # for geometry debugging and data visualization
                            'numba',        # for JIT compilation, will be deprecated in the future
                            'bs4',          # for some other dependency
                            'pandas',       # for data processing, data reading and for htmlcreator
                            'plotly',       # for htmlcreator
                            'scikit-image',
                            'Jinja2', 
                            'pyunpack',
                            'patool'],

        include_package_data = True,
        packages             = find_packages(),
        cmdclass             = {'build_ext': build_ext},
        include_dirs         = [".", np.get_include(), "_random"],

        ext_modules = cythonize(
                                ext_modules, 
                                            annotate = False, # this is getting overriden locally ._.
                                compiler_directives = {
                                                       'profile'        : False, # this is also getting overriden locally ._.
                                                       'language_level' : "3"
                                                      }
                               )
    )  




# https://journals.sagepub.com/doi/suppl/10.1177/ANIB_39_2
#  https://iopscience.iop.org/article/10.1088/1742-6596/1662/1/012021/pdf
# https://iopscience.iop.org/article/10.1088/0031-9155/51/14/017

 
