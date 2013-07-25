:: Setup Visual Studio environment
call setup_vs.bat %1

:: Goto sources directory
cd ..

:: Cleanup
rd /S /Q Release build_fem

:: Build FEM
mkdir build_fem
cd build_fem
cmake -G %generator% -DOGS_FEM=ON -DOGS_DONT_USE_QT=ON -DOGS_PACKAGING=ON -DOGS_PACKAGING_ZIP=ON ..
cmake ..
devenv OGS-FEM-5.sln /Build Release
cmake ..
devenv OGS-FEM-5.sln /Build Release /Project PACKAGE
cd ..
cd scripts