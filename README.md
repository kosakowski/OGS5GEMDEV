# OGS5GEMDEV

This is the OpenGeoSys (version 5) coupling to the GEMS3K thermodynamic solver.

You need boost and boost development libraries to be installed (especially boost-threads).
For the new GEMS3K kernel version one needs nlohmann json, eigen3, pybind11, spdlog, thermofun  library installed (best via conda/mamba).
In an ubuntu 22.04.3 system, the packages `nlohmann-json3-dev` and `libspdlog-dev` have to be installed.

## Basic procedure for compilation:

```bash
git clone https://github.com/kosakowski/OGS5GEMDEV.git
```

```bash
cd OGSGEMDEV/sources
mkdir build 
cd build
cmake .. -DOGS_FEM_GEMS=ON
make
sudo make install
```

optional CMAKE arguments: `-DBOOST_ROOT=/path/to/boost/`

For (sometimes outdated) instructions on usage see the ogs5gem-docs and manual directories. Best is to start with the OGS5GEM specific benchmarks in the benchmarks/C directory. 
