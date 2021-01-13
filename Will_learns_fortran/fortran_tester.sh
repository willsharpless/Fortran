#! /bin/bash
#
gfortran -c -fcheck=all -Wall ft_tester.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran ft_tester.o rkf45_WASedit.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm ft_tester.o
#
mv a.out ft_tester
./ft_tester
# ./rkf45_test > rkf45_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm ft_tester
#
echo "Normal end of execution."