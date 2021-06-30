# Monty Carlo

This is a thesis project, which is yet in development and not ready for packaging. Nevertheless, I uploaded this placeholder package so that the name is not taken in the meantime.

The first version might be available as early as jully. 



## What to expect

### Customizability 

 
### Speed

Although it is a python module this package is written in a happy mix of Python, [Cython](https://cython.org/), C and C++. A notable example of a package that also does this is [Numpy](https://github.com/numpy/numpy). Most of the initialization and pretty much all the programming user interface is in Python, so while setting up your simulation or handling the results of it, you'll be dealing with Python. However, from the moment you tell MontyCarlo to start simulating, it leaves the world of Python and starts running optimized C code. Each language is therefore placed strategically so that it can play to its strenghts.


### Fun

Using the power of [vtk](https://vtk.org/) through the wonderful work of [mayavi](https://pypi.org/project/mayavi/) remarkable visualizations are easy in Monty Carlo. 

50keV electrons in water (secondary particles off):

![Electrons in Water ](https://user-images.githubusercontent.com/63464503/110106080-20e4bc00-7da1-11eb-953c-d5904ff196f1.png)


10MeV electrons in water (primary in red, secondary photon in green)

![image](https://user-images.githubusercontent.com/63464503/110102562-d9f4c780-7d9c-11eb-8f70-20f3b26d3503.png)




![SSSS250k](https://user-images.githubusercontent.com/63464503/110109261-14626280-7da5-11eb-8f0b-cd46bf08fca0.png)


