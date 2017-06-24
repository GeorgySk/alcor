#!/usr/bin/env bash
set -e

cd test_project
# unzipping seed files
unzip -o forts.zip

ffpe_trap_list=$(echo "zero"\
                      "overflow"\
                      "underflow"\
                      "invalid"\
                 | tr [:space:] ,)
fortran_compiler_options=$(echo "-g"\
                                "-fbacktrace"\
                                "-ffpe-trap=${ffpe_trap_list::-1}"\
                                "-Wall"\
                                "-Wextra"\
                                "-Wconversion"\
                                "-fbounds-check"\
                                "-fcheck=all")
# compiling Fortran files
gfortran main.f -o main.e ${fortran_compiler_options}
cd -

python3.6 manage.py "$@"
