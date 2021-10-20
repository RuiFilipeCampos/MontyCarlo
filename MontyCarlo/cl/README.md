# Defining the syntax for the markdown


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


