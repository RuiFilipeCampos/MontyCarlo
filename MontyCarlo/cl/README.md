# Defining the syntax for the markdown


A self closing tag would be something like this

```
<[TAG_NAME] [VAR_NAME0]=[VAR_VALUE0] [VAR_NAME1]=[VAR_VALUE1] ... /> 
```
These guys are meant to get translated to function calls without any positional arguments. For example:

```HTML
<sphere radius=12.2 material={gold} />
```
gets translated to

```Python
sphere[0](
    radius=12.2,
    material=gold,
)
```

In this case, `sphere` is meant to be a tuple only containing the class

```Python
sphere = (Sphere,)
```

A non self closing tag would be, for example

```html
<sphere radius = 10 material = {gold}>
    <sphere radius = 1 material = {water}/>
<sphere />
```

In this example, it gets translated to

```python
sphere[0](
    sphere[0](
        radius = 1,
        material = water,
    ),
    radius = 10,
    material = gold,
)
```

This structure gives me more than enough info to construct the bvh and the csg stuff. The reason for the tuple in the python code is that I want functions to workout like in react:

```python
def custom_shape(*args, **kwargs):
    return (
    <union>
        <sphere radius = 1/>
        <parallelepiped dx=10 dy=1 dz=1
    </union>
    ,)
```








-----


```HTML
<Sphere radius=1 material=gold>
    <Cube dx=1 dy=1 dz=1 material=water />
    <Cube dx=-1 dy=-1 dz=-1 material=water />
</Sphere>
```


```Python
Sphere(
    Cube(
        dx = 1,
        dy = 1,
        dz = 1,
        material = Water,
    ),
    Cube(
        dx = -1,
        dy = -1,
        dz = -1,
        material = Water,
    ),    
    radius = 1,
    material = gold,
)
```


