::: This script will build all the needed libraries in /Libs
::: which will take some hours.
::: The Visual Studio project file will be created in
::: /sources/Build/OGS-5.sln
set USE_VC_EXPRESS=FALSE

::: Create build directory :::
rd /S /Q sources\Build
mkdir sources\Build

call setup_symlinks.bat

cd Libs
call build_all_but_qt_and_vtk_win.bat
cd ..
cd sources\Build
cmake .. -DOGS_FEM=ON %CMAKE_GENERATOR%
cmake ..