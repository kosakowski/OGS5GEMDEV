#! /bin/bash

cd h_tri
echo "Running PETSc benchmark: h_tri"
mpirun -np 3 $1 h_tri &> log.txt 
cd ..

cd McWhorter
echo "Running PETSc benchmark: McWhorter"
mpirun -np 4 $1 mcwt &> log.txt 
cd ..

cd Richards
echo "Running PETSc benchmark: Richards"
mpirun -np 4 $1 h_us_quad &> log.txt 
cd ..

cd T_tri
echo "Running PETSc benchmark: T_tri"
mpirun -np 4 $1 t_tri &> log.txt 
cd ..

cd KueperProblem-PS
echo "Running PETSc benchmark: KueperProblem-PS"
mpirun -np 3 $1 kueper &> log.txt 
cd ..
