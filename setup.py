__doc__ = """Build and/or distribute MontyCarlo.


Usage:

Building and Distributing
=========================

Run this script with the following options:

```
python setup.py --os [OS] --cpu [CPU] --version [VERSION]
```

Where:

- OS in ["mac", "win"]
- CPU in ["intel", "amd"]

BUILDING LOCALLY
================

Run this script with the following options:

```
python setup.py build_ext --inplace
```

An additional flag `-jx`, where `x` is an integer, may be provided
for parallel compilation (a lot of compiling going on, takes a while).  

"""


__author__ = "Rui Campos"


import os
import sys

from setup_version import version
from Cython.Build import cythonize

try:
    from setuptools import setup, find_packages
    from setuptools import Extension
except ImportError:
    from distutils.core import setup, find_packages
    from distutils.extension import Extension


from pathlib import Path
from Cython.Distutils import build_ext
import Cython.Compiler.Options # COMPILER OPTIONS
import numpy as np # need to compile it with the extension modules




EXCLUDED_DIR = ['__pycache__', 'elastic', '.ipynb_checkpoints', 'deprecated']

def iter_dir(parent, directory_list):
    """Recursevly explore the directory structure and accumulate all paths to
    directories that may have .pyx/.pxd files.
    """

    directory_list.append(parent)

    for path in parent.iterdir():
        if path.is_dir():
            if path.name in EXCLUDED_DIR: 
                continue

            iter_dir(path, directory_list)



__PATH__ = Path(".")
src_folder = __PATH__/'MontyCarlo'

directory_list = []
iter_dir(src_folder, directory_list)


args = None

if "--win.amd" in sys.argv:
    lock_commands = True
    sys.argv.remove("--win.amd")

    args = [
             "-O2",             # code optimization
             "-fp:fast",        # math optimization -> changes order of math operations for max efficiency
             "-favor:AMD64"
           ]


if "--win.intel" in sys.argv:
    if args: 
        raise RuntimeError("Incompatible options.")

    sys.argv.remove("--win.intel")

    args =  [
             "-O2",             # code optimization
             "-fp:fast",        # math optimization -> changes order of math operations for max efficiency
             "-favor:INTEL64"
             ]

if "--mac" in sys.argv:
    if args: 
        raise RuntimeError("Incompatible options.")

    sys.argv.remove("--mac")

    args = [
            "-Wno-cpp",
            "-std=c++11",
           ]

if args is None:
    raise RuntimeError("Please specify os/cpu combination.")

"""
OS  = ['win',   'mac']
CPU = ['intel', 'amd']




arg_options = { 
                  'win.amd': [
                            "-O2",             # code optimization
                            "-fp:fast",        # math optimization -> changes order of math operations for max efficiency
                            "-favor:AMD64"
                            ],

                 'win.intel': [
                              "-O2",             # code optimization
                              "-fp:fast",        # math optimization -> changes order of math operations for max efficiency
                              "-favor:INTEL64"
                              ],

                      'mac':  [
                              "-Wno-cpp",
                              "-std=c++11",
                              ],
               }
"""

extra_compile_args = args

def to_python(path):
    """Translates a `Path` object to a string of the form `MontyCarlo.module1`.
    """
    res = path.name

    while path != path.parent:
        path = path.parent
        res = path.name + "." + res

    return res[1:]


# Build Extensions
EXTENSIONS = []


def check_pattern(pattern):
    """Check if there are files with the form "path/*.[ext]"
    """
    return pattern.match(str(pattern))

match = lambda pattern: pattern.match(str(pattern))

for path in directory_list:

    pattern = path/'*.pyx'


    if check_pattern(pattern) is False:
        continue

    print(pattern, True)


    ext = Extension(
                    to_python(path),            
                    [str(pattern)],                
                    extra_compile_args = extra_compile_args,
                    )

    EXTENSIONS.append(ext)



with open("README.md", "r", encoding="utf-8") as file:
    long_description = file.read()

if __name__ == "__main__":
    setup(
            name = "MontyCarlo",
            version = version,
            author = "Rui Filipe de Sousa Campos",
            description = "A fast general purpose monte carlo particle simulator (photons, electrons and positrons). Written in Cython, Python and C++.",
            long_description = long_description,
            long_description_content_type="text/markdown",
            url = "https://github.com/RuiFilipeCampos/MontyCarlo",
            setup_requires   = ['setuptools_scm'],
            install_requires = [
                                 'requests',     # for downloading databases
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
                                ],
            include_package_data = True,
            packages             = find_packages(),
            cmdclass             = {'build_ext': build_ext},
            include_dirs         = [".", np.get_include()],

            ext_modules = cythonize(
                                     EXTENSIONS, 
                                     annotate = False, # this is getting overriden locally ._.
                                     compiler_directives =  {
                                                             'profile'        : False, # this is also getting overriden locally ._.
                                                             'language_level' : "3"
                                                            }
                                    )
         )  





# dia 17 -> 14h30