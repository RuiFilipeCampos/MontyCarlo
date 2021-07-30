# MontyCarlo

Monty Carlo is a Python package that simulates the propagation and effects of ionizing radiation (photons, electrons and positrons with 1keV < E < 1GeV) in matter of homogeneous density, filling CSG models. 

**As of yet, this is an unstable version.** This is a thesis project and, as a student, I am still learning! 

# Installation

It is highly recommended that you install MontyCarlo v0.0.41-pre-alpha on a conda virtual environment containing one of the following python versions, and **nothing else**: 3.7, 3.8 or 3.8. To do so, open an anaconda prompt and run the commands:

```bash
conda create --name py39 python=3.9
conda activate py39
```

The installation steps are simple:

```bash
pip install MontyCarlo
python -c "import MontyCarlo"
```

MyCo will detect that it is the first import and will proceed to download all the necessary databases:

- EADL (\*.txt)
- EPDL (\*.txt)
- EEDL (\*.txt)
- Electron Elastic (\*.npy)
- Positron Elastic (\*.npy)


# A first run !

Once you've installed MontyCarlo, clone the following repository: https://github.com/RuiFilipeCampos/MyCo-EXAMPLE1

Inside this repository folder simply run:

```bash
python main.py
````

Have fun exploring high energy particle tracks in a 3d environment!

- White tracks: Photons
- Blue tracks: electrons
- Red tracks: positrons

![ex01](https://user-images.githubusercontent.com/63464503/124515938-880a8f80-ddd8-11eb-9439-409381b5124a.png)

Be sure to zoom in on every detail! 

![ex02](https://user-images.githubusercontent.com/63464503/124516141-ef284400-ddd8-11eb-9481-099947f7e803.png)

<!---

## What to expect
 
### Speed

Although it is a python module this package is written in a happy mix of Python, [Cython](https://cython.org/), C++. A notable example of a package that also does this is [Numpy](https://github.com/numpy/numpy). Most of the initialization and pretty much all the programming user interface is in Python, so while setting up your simulation or handling the results of it, you'll be dealing with Python. However, from the moment you tell MontyCarlo to start simulating, it leaves the world of Python and starts running optimized C code. Each language is therefore placed strategically so that it can play to its strenghts.


### Fun

Using the power of [vtk](https://vtk.org/) through the wonderful work of [mayavi](https://pypi.org/project/mayavi/) remarkable visualizations are easy in Monty Carlo. 

50keV electrons in water (secondary particles off):

![Electrons in Water ](https://user-images.githubusercontent.com/63464503/110106080-20e4bc00-7da1-11eb-953c-d5904ff196f1.png)


10MeV electrons in water (primary in red, secondary photon in green)

![image](https://user-images.githubusercontent.com/63464503/110102562-d9f4c780-7d9c-11eb-8f70-20f3b26d3503.png)




![SSSS250k](https://user-images.githubusercontent.com/63464503/110109261-14626280-7da5-11eb-8f0b-cd46bf08fca0.png)



## Running a first script

The simplest test is to create a material. For that, create a folder structure like so:

- \my_project
   - \mat
   - \geo
   - main.py
 
In main.py write:

```python 
import MontyCarlo as myco
water = myco.Mat({1:2, 8:1}, 1)
```
This will start compiling all the necessary data to simulate photons, electrons and positrons in water. The first argument is a dictionary of the form

``` 
material = {Z_1:#elements of Z_1
            Z_2:#elements of Z_2
            ...
            Z_n:#elements in Z_n}
```

The second argument is the materials density in g/cm^3.

MyCo will create an output file (\*.html) for debugging purposes (the output file is a work in progress). It will also store the material object in the /mat folder. This way the compilation of a given material is only done once per project. Otherwise, creating an application/simulation for MyCo would be too time consuming. Every time ```myco.Mat({1:2, 8:1}, 1)``` is executed, it will read from the cached file. 


A propper example will be shown here: https://github.com/RuiFilipeCampos/MyCo-EXAMPLE1
-->



# Bugs

This is a very early version of a fairly large code. Bugs are guaranteed! Submitting an [issue](https://github.com/RuiFilipeCampos/MontyCarlo/issues) is a great way to contribute to the project at this stage! 



# Possible Future Work

- Sources
- Tallying
  - Energy Deposition (1d, 2d, 3d, 4d(spatial + temporal) )
  - Flux
  - Others
- Variance reduction 
- Image Detectors
- Extension to E < 1keV (e.g. for laser applications)
- Extension to E > 1GeV
- Implementation of other particles
  - Protons
  - Neutrons
  - etc...
- Dedicated graphics engine (w/sphere tracing)
- An [auto-cad](https://www.autodesk.com/products/autocad/overview) like GUI for CSG modeling
- [Geant4](https://github.com/Geant4/geant4) like API
- GPU accelaration
- CPU multiprocessing/multithreading
- Advanced data vizualization (w/ [ParaView](https://www.paraview.org/))
- Distributed Cloud Computing
- Dedicated python notebook (like Jupyter)
