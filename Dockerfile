ARG PYTHON3_VERSION="3"

FROM python:${PYTHON3_VERSION}

RUN apt-get update && \
    apt-get install -y gfortran \
                       unzip

WORKDIR /alcor

COPY ./requirements.txt /alcor/requirements.txt
RUN python3 -m pip install -r requirements.txt

COPY . /alcor/
RUN python3 -m pip install .

ENTRYPOINT ["python3", "manage.py"]
