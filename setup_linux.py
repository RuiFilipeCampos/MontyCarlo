"""
Build Monty Carlo. This script is called from the Docker File
"""

import os

from setup_version import version

try:
    from setuptools import setup, find_packages
    from setuptools import Extension
except ImportError:
    from distutils.core import setup, find_packages
    from distutils.extension import Extension


from Cython.Build import cythonize
from Cython.Distutils import build_ext
import Cython.Compiler.Options # COMPILER OPTIONS

# Numpy is compiled with the extension modules.
import numpy as np 

# MACOS gcc args
args = ["-Wno-cpp", "-std=c++11"]

ext_modules = [

    # Math Helpers such as **highly** optimized interpolators.
    Extension(
        "tools.*",
        ["MontyCarlo/tools/*.pyx"],
        extra_compile_args = args,
        language = "c++"
    ),

    # Implementation of the random sampling of each particle
    Extension(
        "particles.*",
        ["MontyCarlo/particles/*.pyx"], 
        extra_compile_args = args,
        language = "c++"
    ),

    # 
    Extension(
        "*",
        ["MontyCarlo/*.pyx"],
        extra_compile_args = args,
        language = "c++"
    ), 

    # Geometry Handlers -> CSG using forward ray tracing techinique mixed with signed distance funciton based ray marcher
    Extension(
        "geometry.*",
        ["MontyCarlo/geometry/*.pyx"],
        extra_compile_args = args,
        language = "c++"
    ),


    # Data generation and processing for eletrons.
    Extension(
        "materials.electron.*",
        ["MontyCarlo/materials/electron/*.pyx"],
        extra_compile_args = args,
        language = "c++"
    ),

    # Data generation and processing for positrons.
    Extension(
        "materials.positron.*",
        ["MontyCarlo/materials/positron/*.pyx"],
        extra_compile_args = args,
        language = "c++"
    ),

    # Data generation and processing. pyRelax is included here.
    Extension(
        "materials.*",
        ["MontyCarlo/materials/*.pyx"],
        extra_compile_args = args,
        language = "c++"
    ),

    # Data generation and processing for photons.
    Extension(
        "materials.photon.*",
        ["MontyCarlo/materials/photon/*.pyx"],
        extra_compile_args = args,
        language = "c++"
    ),

    # External libraries (notable package: MIXMAX, same pRNG used in GEANT4)
    Extension(
        "external.*",
        ["MontyCarlo/external/*.pyx"],
        extra_compile_args = args,
        language = "c++"
    )
]
 




if __name__ == "__main__":

    with open("README.md", "r", encoding="utf-8") as fh:
        long_description = fh.read()


    setup(

        entry_points = {
            'console_scripts': ['myco=MontyCarlo.cl.myco:main'],
        },
        
        
        name = "MontyCarlo",
        version = version,
        author = "Rui Filipe de Sousa Campos",
        description = "A fast general purpose monte carlo particle simulator"
                      " (photons, electrons and positrons). Written in Cython," 
                      " Python and C++.",
        long_description = long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/RuiFilipeCampos/MontyCarlo",


        setup_requires   = ['setuptools_scm'],
        install_requires = [
            'requests',     # for downloading databases
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
            'patool',
            'colorama',
            'termcolor',
        ],
        include_package_data = True,
        packages             = find_packages(),
        cmdclass             = {'build_ext': build_ext},
        include_dirs = [
            ".", 
            np.get_include(), 
        ],
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

 
