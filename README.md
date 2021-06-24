# MontyCarlo

## Instalation

It is now possible to install an unstested version of MontyCarlo (0.0.34). The instalation consists in two simple steps.

```
pip install MontyCarlo
```

Which should take a while. The second step consists in just doing a first import:

```python 
import MontyCarlo
```

MyCo will detect that it is the first import and will proceed to download all the necessary databases:

- EADL (\*.txt)
- EPDL (\*.txt)
- EEDL (\*.txt)
- Electron Elastic (\*.npy)
- Positron Elastic (\*.npy)

## Running a first script

The simplest test is to create a material. For that, create a folder structure like so:

- my_project
   - mat
   - geo
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

### Bugs

This is a very early version of a fairly large code. Bugs are guaranteed! I would very much appreaciate if you report them to me. 
