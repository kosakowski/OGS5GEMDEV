:: Setup Visual Studio environment
call setup_vs.bat %1

:: Goto sources directory
cd ..

:: Cleanup
rd /S /Q Release build_fem build_gems build_pqc build_brns build_lis

:: Executables will copied to Release directory
mkdir Release

:: Build FEM
mkdir build_fem
cd build_fem
cmake -G %generator% -DOGS_FEM=ON -DOGS_DONT_USE_QT=ON ..
cmake ..
devenv OGS-FEM-5.sln /Build Release
cd bin\Release
copy /Y ogs.exe ..\..\..\Release\
cd ..\..\..\


:: Build FEM_GEMS
mkdir build_gems
cd build_gems
cmake -G %generator% -DOGS_FEM_GEMS=ON -DOGS_DONT_USE_QT=ON ..
cmake ..
devenv OGS-FEM-5-GEMS.sln /Build Release
cd bin\Release
ren ogs.exe ogs_gems.exe
copy /Y ogs_gems.exe ..\..\..\Release\
cd ..\..\..\

:: Build FEM_PQC
mkdir build_pqc
cd build_pqc
cmake -G %generator% -DOGS_FEM_PQC=ON -DOGS_DONT_USE_QT=ON ..
cmake ..
devenv OGS-FEM-5-PQC.sln /Build Release
cd bin\Release
ren ogs.exe ogs_pqc.exe
copy /Y ogs_pqc.exe ..\..\..\Release\
cd ..\..\..\

:: Build FEM_BRNS
mkdir build_brns
cd build_brns
cmake -G %generator% -DOGS_FEM_BRNS=ON -DOGS_DONT_USE_QT=ON ..
cmake ..
devenv OGS-FEM-5-BRNS.sln /Build Release
cd bin\Release
ren ogs.exe ogs_brns.exe
copy /Y ogs_brns.exe ..\..\..\Release\
cd ..\..\..\

:: Build FEM_LIS
mkdir build_lis
cd build_lis
cmake -G %generator% -DOGS_FEM_LIS=ON -DOGS_DONT_USE_QT=ON ..
cmake ..
devenv OGS-FEM-5-LIS.sln /Build Release
cd bin\Release
ren ogs.exe ogs_lis.exe
copy /Y ogs_lis.exe ..\..\..\Release\
cd ..\..\..\

:: Build FEM_MPI
mkdir build_mpi
cd build_mpi
cmake -G %generator% -DOGS_FEM_MPI=ON -DOGS_DONT_USE_QT=ON -DMPI_LIBRARY:FILEPATH="C:/Program Files (x86)/MPICH2/lib/mpi.lib" -DMPI_INCLUDE_PATH:PATH="C:/Program Files (x86)/MPICH2/include" -DMPIEXEC:FILEPATH="C:/Program Files (x86)/MPICH2/bin/mpiexec.exe" ..
cmake ..
devenv OGS-FEM-5-MPI.sln /Build Release
cd bin\Release
ren ogs.exe ogs_mpi.exe
copy /Y ogs_mpi.exe ..\..\..\Release\
cd ..\..\..\

cd scripts