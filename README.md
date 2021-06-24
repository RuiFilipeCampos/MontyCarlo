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

