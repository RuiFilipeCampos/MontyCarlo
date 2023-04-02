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


def make_extension(import_expression, path, extra_compile_args = ["-Wno-cpp", "-std=c++11", "-Wno-format"], language = "c++"):
    
    extension = Extension(
        import_expression, 
        [path], 
        extra_compile_args = extra_compile_args,
        language = language
    )

    return extension



ext_modules = [
    make_extension("tools.*", "MontyCarlo/tools/*.pyx"),
    make_extension("particles.*", "MontyCarlo/particles/*.pyx"),
    make_extension("*", "MontyCarlo/*.pyx"),
    make_extension("geometry.*", "MontyCarlo/geometry/*.pyx"),
    make_extension("materials.electron.*", "MontyCarlo/materials/electron/*.pyx"),
    make_extension("materials.positron.*", "MontyCarlo/materials/positron/*.pyx"),
    make_extension("materials.*", "MontyCarlo/materials/*.pyx"),
    make_extension("materials.photon.*", "MontyCarlo/materials/photon/*.pyx"),
    make_extension("external.*", "MontyCarlo/external/*.pyx"),
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

 
