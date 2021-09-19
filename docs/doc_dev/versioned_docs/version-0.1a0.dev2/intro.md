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

```Bash
conda create --name myco_env python=3.9
```

This will create an environment with a clean python installation.

Use the following command to activate the virtual environment:

```Bash
conda activate myco_env
```

Now you can install MontyCarlo as you would install any other package:

```Bash
pip install MontyCarlo
```

The first time you import MontyCarlo, it will automatically download all the necessary databases. This is a one time step:

```Bash
python -c "import MontyCarlo"
```

## Running the Demo

Either download or clone the following repository: [MyCo-EXAMPLE1](https://github.com/RuiFilipeCampos/MyCo-EXAMPLE1).

Navigate to the root folder of the `MyCo-EXAMPLE1` project and run the following command:

```Bash
python main.py
```

This demo should not take too long. It will open a plot of the trajectory of the particles that you can explore by zooming in and moving around via `shift + mouse dragg`.

