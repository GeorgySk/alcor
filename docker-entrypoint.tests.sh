#!/usr/bin/env bash

f2py -c ./fortran_tests/tables/da_cooling.f90 -m da_cooling
f2py -c ./fortran_tests/tables/db_cooling.f90 -m db_cooling
f2py -c ./fortran_tests/tables/da_color.f90 -m da_color
cp ./fortran_tests/tables/fort_files/* .

python3 setup.py test

rm *.so
rm fort.*
