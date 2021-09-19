---
sidebar_position: 1
---



# Tutorial

Let's discover **MontyCarlo in less than 5 minutes**.


:::danger Disclamer
As of now, MontyCarlo is still in pre-alpha development phase.
There is no guarantees given in terms of stability and numerical results until version `0.1` is released.
:::

:::tip However !
It is possible to try it out!
:::

## What does MontyCarlo do?

- MontyCarlo simulates the propagation of high energy particles through media of uniform density.
- It gives you the option to create CSG models that you can fill with any material.
- And it allows the calculation of the energy deposition of the particles in the various media.


## Installation

MontyCarlo is meant to be installed in a **python virtual environment**. 

You can use whatever environment you want, but in this tutorial we'll be using Annaconda Prompt to create and use virtual environments.

**Hardware Requirements**: 64bit CPU

**Software Requirements**:
- One of these python versions: 3.7, 3.8 or 3.9;
- One of these operating systems: Windows or MacOS;

Open an Annaconda Prompt and create a new environment:

```bash
conda create --name myco_env python=3.9
```

This will create an environment with a clean python installation.

Use the following command to activate the virtual environment:

```bash
conda activate myco_env
```

Now you can install MontyCarlo as you would install any other package:

```bash
pip install MontyCarlo
```

The first time you import MontyCarlo, it will automatically download all the necessary databases. This is a one time step:

```bash
python -c "import MontyCarlo"
```

## Running the Demo

Either download or clone the following repository: [MyCo-EXAMPLE1](https://github.com/RuiFilipeCampos/MyCo-EXAMPLE1).

Navigate to the root folder of the `MyCo-EXAMPLE1` project and run the following command:

```bash
python main.py
```

This demo should not take too long. It will open a plot of the trajectory of the particles that you can explore by zooming in and moving around via `shift + mouse dragg`.


## Creating your own project

To create a new MontyCarlo project called `my_project`, run the following command:

```bash
myco new my_project
```

It will create a folder named `/my_project` with the following structure:

- /my_project
	- /mat
	- /geo
	- /out
	\_\_init\_\_.py
	main.py

The `/mat` and `/geo` folders are for caching purposes. When you instruct MontyCarlo to create a new `Material` instance, it will start compiling all the necessary information for simulating in the material you've defined. This takes quite a while to do! To save you the trouble of compiling over and over again each time you run your application, the `Material` instance is saved into the `/mat` folder. The `/geo` folder follows the same logic, but for rendering the CSG models you construct.

The `/out` folder is where the output files go to when you create a new material. These files hold all the information about the material you've constructed and is meant for debugging purposes.

The `main.py` is where you start writing your application.

## Writing a simple `main.py` script

Let's start by importing:

- the top level MontyCarlo module -> configuration stuff
- the constructive solid geometry module -> for building the geometry
- the sources module -> for defining your sources
- the `MatGen` function -> for creating materials

```python
import MontyCarlo as myco
import MontyCarlo.geometry.CSG as csg
import MontyCarlo.sources as src

from MontyCarlo.materials import MatGen
```

The stochiometric formula of water is $H_2 O$, which is represented as the following dictionary:

```python
water_formula = {
	1:2 #H:2
	8:1 #O:1
}
```

The density of water is approximately $\textrm{1g/cm}^3$. To create a `Material` instance,

```python
water = MatGen(water_formula, 1, name="Water(Liquid)")
```

:::tip Important
`MatGen` will first look into the `/mat` folder to check if there is any material matching your parameters. If there is, it create the `Material` instance from that file, otherwise it will compile a new material.
:::



