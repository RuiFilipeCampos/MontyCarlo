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
