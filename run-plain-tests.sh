#!/usr/bin/env bash

set -ex

cd tests/tables
# unzipping Fortran link files
unzip -o fort_files.zip
cd -

python3 setup.py test
