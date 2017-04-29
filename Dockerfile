FROM python:3.5

WORKDIR /alcor
COPY . /alcor/
RUN python3 -m pip install .
