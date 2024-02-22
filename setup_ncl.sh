#!/bin/bash

conda create -n ncl_stable -c conda-forge ncl
source activate ncl_stable
ncl -V

# ncl xxx.ncl
