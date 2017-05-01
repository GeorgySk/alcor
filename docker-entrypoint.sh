#!/usr/bin/env bash
set -e

cd test_project
# unzipping seed files
unzip -o forts.zip
# compiling Fortran files
gfortran main.f -o main.e
cd -

python3.6 manage.py "$@"
