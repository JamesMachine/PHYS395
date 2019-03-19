#!/bin/bash

gfortran -o3 -fdefault-real-8 newtmet.f90 -o newtmet
./newtmet > data_newtmet

gfortran -o3 -fdefault-real-8 goldenmin.f90 -o goldenmin
./goldenmin > data_goldenmin

gfortran -o3 -fdefault-real-8 levmar.f90 -o levmar -llapack
./levmar 3 < a2.dat > data_levmar

gnuplot plot.gpl

