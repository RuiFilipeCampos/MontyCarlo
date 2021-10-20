# Defining the syntax for the markdown


A self closing tag would be something like this

```HTML
<[TAG_NAME] [VAR_NAME0]=[VAR_VALUE0] [VAR_NAME1]=[VAR_VALUE1] ... /> 
```


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


