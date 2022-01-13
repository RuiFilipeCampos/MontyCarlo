FROM python:3.9.9

# Not necessarily needed, just in case...
RUN /usr/local/bin/python -m pip install --upgrade pip

# MontyCarlo *build* depenencies...
RUN pip install cython
RUN pip install numpy
RUN pip install wheel
RUN pip install setuptools

# Cloning the repository...
RUN git clone https://github.com/RuiFilipeCampos/MontyCarlo.git

# Building Monty Carlo...
RUN cd MontyCarlo && python setup_linux.py >> output.txt
