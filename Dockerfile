FROM python:3.9.9

RUN mkdir app

# Not necessarily needed, just in case...
RUN /usr/local/bin/python -m pip install --upgrade pip

# MontyCarlo *build* depenencies...
RUN pip install cython
RUN pip install wheel
RUN pip install numpy
RUN pip install setuptools

COPY MontyCarlo app/MontyCarlo
COPY setup_linux.py app
COPY requirements.txt app
COPY setup_version.py app
COPY setup.cfg app
COPY README.md app





# Building Monty Carlo...

RUN ls && python app/setup_linux.py build_ext -j6 --inplace
RUN pip install -r requirements.txt

WORKDIR /app


