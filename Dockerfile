FROM python:3.9.9

RUN mkdir app

RUN /usr/local/bin/python -m pip install --upgrade pip

RUN apt-get update && apt-get install -y python3-opencv
RUN pip install opencv-python

RUN pip install cython
RUN pip install wheel
RUN pip install setuptools

COPY requirements.txt app
RUN cd app && pip install -r requirements.txt

COPY MontyCarlo app/MontyCarlo
COPY setup_linux.py app
COPY setup_version.py app
COPY setup.cfg app/setup.cfg
COPY README.md app

WORKDIR /app


# Building Monty Carlo...
RUN python setup_linux.py build_ext -j6 --inplace
ENV PYTHONPATH "${PYTHONPATH}:/app"
RUN python -c "import MontyCarlo"

# Setting up jupyter
# Install Jupyter
RUN pip install jupyter

# Expose the Jupyter port
EXPOSE 8888

RUN mkdir /notebooks
RUN mkdir /notebooks/geo
RUN mkdir /notebooks/mat
RUN mkdir /notebooks/out

WORKDIR /notebooks

RUN apt install -y libgl1-mesa-glx xvfb
RUN pip install panel

COPY jupyterStyle.css /app/jupyterStyle.css
RUN jupyter notebook --generate-config
RUN echo "c.NotebookApp.extra_static_paths = ['/app/jupyterStyle.css']" | tee -a ~/.jupyter/jupyter_notebook_config.py


# Start Jupyter Notebook
# CMD ["jupyter", "notebook", "--ip='0.0.0.0'", "--port=8888", "--no-browser", "--allow-root"]

# jupyter notebook --ip='0.0.0.0' --port=8888 --no-browser --allow-root