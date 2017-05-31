FROM python:3

RUN apt-get update && \
    apt-get install -y gfortran \
                       unzip && \
#   bokeh requirement:
    apt-get install -y nodejs

WORKDIR /alcor

COPY ./requirements.txt /alcor/requirements.txt
RUN python3 -m pip install -r requirements.txt

COPY . /alcor/
RUN python3 -m pip install .

ENTRYPOINT ["python3", "manage.py"]
