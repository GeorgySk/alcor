#!/usr/bin/env bash
set -e

cd test_project
# compiling Fortran files
gfortran main.f -o main.e
cd -

python3.6 manage.py "$@"
