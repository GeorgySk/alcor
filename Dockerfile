FROM python:3

WORKDIR /alcor

COPY ./requirements.txt /alcor/requirements.txt
RUN python3 -m pip install -r requirements.txt

COPY . /alcor/
RUN python3 -m pip install .

ENTRYPOINT ["python3", "manage.py"]
