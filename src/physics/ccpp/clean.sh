#!/bin/bash

./framework/scripts/ccpp_prebuild.py \
  --config host/ccpp_prebuild_config.py \
  --clean

rm host/CCPP_STATIC_API.cmake
