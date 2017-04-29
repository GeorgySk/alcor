FROM python:3

WORKDIR /alcor
COPY . /alcor/
RUN python3 -m pip install .
