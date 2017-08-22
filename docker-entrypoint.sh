#!/usr/bin/env bash
set -e

cd test_project
# unzipping seed files
unzip -o forts.zip

gfortran main.f -o main.e ${FORTRAN_COMPILER_OPTIONS}
cd -

python3.6 manage.py "$@"
