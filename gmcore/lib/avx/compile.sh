#!/bin/bash

source ../../../env.sh

mpiicx -mavx512f -c avx_interface.cpp -o avx.o 
