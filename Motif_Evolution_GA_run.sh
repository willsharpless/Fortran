#! /bin/bash
#
gfortran -g -c -fbacktrace -fcheck=all -w Motif_Evolution_GA.f90 
# -lmwlapack -lmwblas
#-framework Accelerate
# Might need to do Accelerate instead of vecLib
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran Motif_Evolution_GA.o M_scramble.o rkf45_WASedit.o -framework Accelerate
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm Motif_Evolution_GA.o
#
mv a.out Motif_Evolution_GA
# ./Motif_Evolution_GA

for i in 1 2 3
do
   ./Motif_Evolution_GA > Evolution_runs/Motif_Evolution_300its_run${i}_010721.txt
done

# ./rkf45_test > rkf45_test.txt
# if [ $? -ne 0 ]; then
#   echo "Run error."
#   exit
# fi
# # rm rkf45_test
# #
# echo "Normal end of execution."

# LAPACK and BLAS, Linearalgebra, BNNS, BigNum, Sparse,
#/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework