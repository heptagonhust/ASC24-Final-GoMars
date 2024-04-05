#!/bin/bash

cd "$(dirname $0)" || exit 1
set -e

source ./env.sh

gptl_dir="$(pwd)/gmcore/gptl"
mkdir -p "$gptl_dir"

cd "$gptl_dir"
pwd
if [ x"$(ls -A .)" = x"" ]; then
  # wget https://github.com/jmrosinski/GPTL/releases/download/v8.1.1/gptl-8.1.1.tar.gz
  # tar -zxvf gptl-8.1.1.tar.gz
  cp -r /data/gomars_data/gptl-8.1.1 .
fi

cd gptl-8.1.1
# wget https://gist.githubusercontent.com/bonfus/21dec6b966859f5f509b935f8b055a7f/raw/macros.make
# ./configure --enable-none --disable-openmp --prefix=$gptl_dir --enable-shared=no --enable-static=yes --enable-pmpi
# delete '--enable-none' because configure: WARNING: unrecognized options: --enable-none
./configure --disable-openmp --prefix=$gptl_dir --enable-shared=no --enable-static=yes --enable-pmpi

# make check
make install
# wordaround: remove all shared libs
# rm ${gptl_dir}/lib/*so*
