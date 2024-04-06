#!/bin/bash

set -e
cd "$(dirname $0)" || exit 1

source ./env_xyw.sh

cd gmcore

if [ ! -d fox ]; then
  # git clone https://github.com/andreww/fox.git
  cp -r /data/gomars_data/fox ./
fi

cd fox

export FC=mpiifx
export FCFLAGS="-O3 -g"

./configure --enable-fast --enable-sax --enable-wcml --enable-wxml --enable-wkml --enable-dom --prefix=$(pwd)

make clean

make -j64
make check
make install