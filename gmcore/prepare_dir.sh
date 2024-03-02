#!/bin/bash
# lxy: I tried to do this in select_run.sh, but it didn't work while run on sbatch. So I'm doing it here. 
if [ ! -d "output" ]; then
    mkdir "output"
fi

if [ ! -d "run" ]; then
    mkdir "run"
fi

pushd run
if [ ! -d "GMCORE-TESTBED" ]; then
    mkdir "GMCORE-TESTBED"
    cp -r /data/gomars_libs/GMCORE-TESTBED/* GMCORE-TESTBED
fi
popd