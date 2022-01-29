FROM python:3.9.9

# Not necessarily needed, just in case...
RUN /usr/local/bin/python -m pip install --upgrade pip

# MontyCarlo *build* depenencies...
RUN pip install cython
RUN pip install wheel
RUN pip install setuptools
RUN pip install numpy

# Cloning the repository...
RUN git clone --depth 1 -b pre-alpha/master https://github.com/RuiFilipeCampos/MontyCarlo.git

# Building Monty Carlo...
RUN cd MontyCarlo && python setup_linux.py build_ext --inplace >> output.txt
