The first stable realease will be version `0.1`.

# MontyCarlo (v0.1a0.dev3)
![](https://img.shields.io/github/v/release/RuiFilipeCampos/MontyCarlo?include_prereleases) ![license](https://img.shields.io/github/license/RuiFilipeCampos/MontyCarlo) ![pyversion](https://img.shields.io/badge/python-%3E%3D3.7-blue) ![architecture](https://img.shields.io/badge/architecture%20-64--bit-blue) ![os](https://img.shields.io/badge/OS-win%2Fmac-blue)

**MontyCarlo** is a python framework for setting up simulations and/or developing applications whose basis is the simulation of particle transport. It simulates the propagation and effects of ionizing radiation ([photons](https://en.wikipedia.org/wiki/Photon), [electrons](https://en.wikipedia.org/wiki/Electron) and [positrons](https://en.wikipedia.org/wiki/Positron) with energies between 1keV and 1GeV) in matter of homogeneous density, filling [constructive solid geometry](https://en.wikipedia.org/wiki/Constructive_solid_geometry) models.

**As of yet, this is an unstable version.** This is a thesis project and, as a student, I am still learning!

- **Documentation** ([work in progress](https://github.com/RuiFilipeCampos/MontyCarlo/tree/docs)): https://ruifilipecampos.github.io/MontyCarlo/
- **PyPI**: https://pypi.org/project/MontyCarlo/
- *Short Overview*: https://github.com/RuiFilipeCampos/MontyCarlo/wiki

This work has a poster presentation in the [3rd European Congress of Medical Physics](https://www.ecmp2020.org/) and has been presented in a [workshop](https://ruifilipecampos.github.io/MontyCarlo/Poster_workshop_medical_physics.pdf) organized by the Faculty of Sciences of the University of Porto and the Ludwig Maximilian University of Munich.

<img src="https://user-images.githubusercontent.com/63464503/127783224-295ea39e-b935-4cbd-b4d9-a1012bc12729.jpg" width="auto" height="250">


<!--
# Installation

It is highly recommended that you install MontyCarlo 0.1a0.dev1 on a conda virtual environment containing one of the following python versions, and **nothing else**: 3.7, 3.8 or 3.9. To do so, open an anaconda prompt and run the commands:

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

-->
# Try it out!

See [here](https://ruifilipecampos.github.io/MontyCarlo/docs/intro#installation) for instructions on installing the latest release.

Once you've installed MontyCarlo, clone or download the following repository: https://github.com/RuiFilipeCampos/MyCo-EXAMPLE1

Inside this repository folder simply run:

```bash
python main.py
````

Have fun exploring high energy particle tracks in a 3d environment!

- White tracks: Photons
- Blue tracks: electrons
- Red tracks: positrons

The innermost sphere contains water, the outer sphere contains air and the rest of space is filled with gold.

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


# Available Features:

- Construction of any material via a stochiometric formula and density `water = Mat({1:2, 8:1}, 1)`;
- Constructed materials are automatically cached in the folder `your_project\mat`.
- Only spheres are available. This will remain as such until all this has been thoroughly tested:
    - Constructive Solid Geometry (CSG) using the `|` `&` and `-` operators;
    - linear transformations on the volumes (translation and rotation);
    - bounding volume hierarchy (BVH) constructed with the aid of the user;
    - a syntatic indication of the BVH using `with` statements;
    - a new method of particle transport that greatly accelerates the simulation of electrons and positrons;
- The volumes surfaces are rendered and cached in `your_project/geo`;
- Three particles are available:
   - **Photons** (analogue simulation);
     - Compton Scattering;
     - Rayleigh Scattering;
     - Photoelectric Effect;
     - Pair Production;
     - Triplet Production;
   - **Electrons** (class II condensed history);
     - Elastic Scattering (atom is not affected): Angular Deflection + Bremstrahlung Production;
     - Inelastic Scattering (atom is affected): Interaction with an individual atom + with the condensed medium as a whole;
   - **Positrons** (class II condensed history);
     - Elastic Scattering (atom is not affected): Angular Deflection + Bremstrahlung Production;
     - Inelastic Scattering (atom is affected): Interaction with an individual atom + with the condensed medium as a whole;
     - Anihilation (positron meets electron);
- The simulation is coupled (e.g. supports secondary particle creation)
- Supports simulation of post-ionization relaxation effects;
- Two particle sources are available:
   - Isotropic point source: emits particles from a point with randomized directions - `IsotropicPoint`
   - Directional point source: emits particles from a point towards a specified direction - `Beam`
- Automated database download on first import;
- 3d plotting of particle trajectories;
- 3d plotting of the constructed geometry;
- simultaneous plotting of both geometry and trajectories;
- One tally is available:
   - `Z_TALLY` - calculates PDD's
- Automatic generation of \*.html output files (work in progress though)

# Possible Future Work

- Sources
- Tallying
  - Energy Deposition (1d, 2d, 3d, 4d(spatial + temporal) )
  - Flux
  - Others
- Variance reduction 
- Image Detectors
- Extension to E < 1keV (for laser applications)
- Extension to E > 1GeV (for thermonuclear applications)
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
