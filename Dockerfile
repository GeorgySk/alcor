ARG PYTHON3_VERSION=3

FROM python:${PYTHON3_VERSION}

RUN apt-get update && \
    apt-get install -y gfortran \
                       unzip \
                       # scipy's dependencies
                       libblas-dev \
                       liblapack-dev \
                       libatlas-base-dev

WORKDIR /alcor

ARG FORTRAN_COMPILER_OPTIONS
ENV FORTRAN_COMPILER_OPTIONS=${FORTRAN_COMPILER_OPTIONS}

COPY ./test_project test_project
RUN cd test_project && \
    # unzipping seed files
    unzip -o forts.zip && \
    gfortran main.f -o main.e ${FORTRAN_COMPILER_OPTIONS}

COPY ./requirements.txt requirements.txt
RUN python3 -m pip install -r requirements.txt

COPY ./alcor alcor
COPY ./tests tests
COPY ./README.rst README.rst
COPY ./setup.py setup.py
COPY ./setup.cfg setup.cfg
RUN python3 -m pip install .

COPY ./manage.py manage.py

RUN cd tests/tables && \
    unzip -o fort_files.zip

COPY ./docker-entrypoint.sh docker-entrypoint.sh
ENTRYPOINT ["/alcor/docker-entrypoint.sh"]
