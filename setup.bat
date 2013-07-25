::: The Visual Studio project file will be created in
::: /sources/Build/OGS-FEM-5.sln
set USE_VC_EXPRESS=FALSE

::: Create build directory :::
rd /S /Q sources\Build
mkdir sources\Build

cd sources\Build
cmake .. -DOGS_FEM=ON %CMAKE_GENERATOR%
cmake ..