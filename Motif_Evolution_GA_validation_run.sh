#! /bin/bash
#
gfortran -g -c -fbacktrace -fcheck=all -w Motif_Evolution_GA_validation.f90 
# -lmwlapack -lmwblas
#-framework Accelerate
# Might need to do Accelerate instead of vecLib
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran Motif_Evolution_GA_validation.o M_scramble.o rkf45_WASedit.o -framework Accelerate
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm Motif_Evolution_GA_validation.o
#
mv a.out Motif_Evolution_GA_validation
./Motif_Evolution_GA_validation > Motif_Evolution_GA_validation.txt

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