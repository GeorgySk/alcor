#!/usr/bin/env bash

set -ex

cd test_project
# unzipping seed files
unzip -o forts.zip

gfortran main.f -o main.e ${FORTRAN_COMPILER_OPTIONS}
cd -

python3 manage.py "$@"
