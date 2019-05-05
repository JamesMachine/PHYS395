#!/bin/bash
echo complile problem1... 
gfortran -O3 -fdefault-real-8 prob1_even.f90 -o a && ./a
gfortran -O3 -fdefault-real-8 prob1_odd.f90 -o a && ./a
python plot1.py

echo complile problem2... 
echo list of different energy states
gfortran -O3 -fdefault-real-8 prob2.f90 -o a && ./a
python plot2.py

echo complile problem3... 
echo list of different energy states
gfortran -O3 -fdefault-real-8 prob3.f90 -o a && ./a
python plot3.py

echo complile problem4... 
echo list of different energy states
gfortran -O3 -fdefault-real-8 prob4.f90 -llapack -o a && ./a
python plot4.py

echo complile problem5... 
echo list of different energy states
gfortran -O3 -fdefault-real-8 prob5.f90 -llapack -o a && ./a
python plot5.py

