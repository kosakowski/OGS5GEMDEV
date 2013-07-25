#!/bin/bash

cd sources/
mkdir Build
cd Build/
cmake .. -DOGS_FEM=ON
make

