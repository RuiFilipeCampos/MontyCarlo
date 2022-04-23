FROM python:3.9.9

RUN mkdir app

# Not necessarily needed, just in case...
RUN /usr/local/bin/python -m pip install --upgrade pip



COPY MontyCarlo app/MontyCarlo
COPY setup_linux.py app
COPY requirements.txt app
COPY setup_version.py app
COPY setup.cfg app/setup.cfg
COPY README.md app


RUN apt-get update && apt-get install -y python3-opencv
RUN pip install opencv-python

# MontyCarlo *build* depenencies...
RUN pip install cython
RUN pip install wheel
RUN pip install setuptools

RUN cd app && pip install -r requirements.txt


# Building Monty Carlo...
RUN cd app && ls && python setup_linux.py build_ext -j6 --inplace
RUN cd app && python -c "import MontyCarlo"

WORKDIR /app
