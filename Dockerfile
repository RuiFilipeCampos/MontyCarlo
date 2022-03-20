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
COPY server /server


# Building Monty Carlo...


RUN ls && python setup_linux.py build_ext -j6 -b ./server/engine/
RUN cd server && pip install -r requirements.txt

WORKDIR /server
EXPOSE 1000
CMD ["python", "manage.py", "runserver", "0.0.0.0:1000"]

# EXPOSE 1000
# STOPSIGNAL SIGTERM

# RUN cd server && python manage.py runserver 1000
