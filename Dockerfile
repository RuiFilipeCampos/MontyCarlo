FROM python:3.9.9

# Not necessarily needed, just in case...
RUN /usr/local/bin/python -m pip install --upgrade pip

# MontyCarlo *build* depenencies...
RUN pip install cython
RUN pip install wheel
RUN pip install numpy
RUN pip install setuptools

COPY MontyCarlo /MontyCarlo
COPY setup_linux.py .
COPY requirements.txt .
COPY README.md .
COPY setup_version.py .
COPY setup.cfg .


# Building Monty Carlo...

RUN mkdir app


RUN ls && python setup_linux.py build_ext -j6 -b ./app/
RUN cd app && pip install -r requirements.txt

WORKDIR /app


