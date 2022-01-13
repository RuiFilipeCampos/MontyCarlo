FROM python:3.9.9
RUN /usr/local/bin/python -m pip install --upgrade pip

RUN git clone https://github.com/RuiFilipeCampos/MontyCarlo.git
RUN cd MontyCarlo && python setup_linux.py
