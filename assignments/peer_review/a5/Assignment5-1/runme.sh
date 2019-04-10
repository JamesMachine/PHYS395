#!/bin/bash

gfortran -O3 -fdefault-real-8 Q1odd.f90 -o Q1odd
./Q1odd > Q1odddata

gfortran -O3 -fdefault-real-8 Q1odd2.f90 -o Q1odd2
./Q1odd2 > Q1odd2data

gfortran -O3 -fdefault-real-8 Q1even.f90 -o Q1even
./Q1even > Q1evendata

gfortran -O3 -fdefault-real-8 Q1even2.f90 -o Q1even2
./Q1even2 > Q1even2data

gfortran -O3 -fdefault-real-8 Q2.f90 -o Q2
./Q2

gfortran -O3 -fdefault-real-8 Q2evenplot.f90 -o Q2evenplot
./Q2evenplot > Q2evenplotdata

gfortran -O3 -fdefault-real-8 Q2oddplot.f90 -o Q2oddplot
./Q2oddplot > Q2oddplotdata

gfortran -O3 -fdefault-real-8 Q2evenplotsq.f90 -o Q2evenplotsq
./Q2evenplotsq > Q2evenplotsqdata

gfortran -O3 -fdefault-real-8 Q2oddplotsq.f90 -o Q2oddplotsq
./Q2oddplotsq > Q2oddplotsqdata

gfortran -O3 -fdefault-real-8 Q3.f90 -o Q3
./Q3

gfortran -O3 -fdefault-real-8 Q3plot.f90 -o Q3plot
./Q3plot > Q3plotdata

gfortran -O3 -fdefault-real-8 text.f90 -o text
./text

gfortran -O3 -fdefault-real-8 Q4v.f90 -llapack -o Q4v
./Q4v

gfortran -O3 -fdefault-real-8 Q4plot.f90 -llapack -o Q4plot
./Q4plot > Q4plotdata

gfortran -O3 -fdefault-real-8 Q4plotsq.f90 -llapack -o Q4plotsq
./Q4plotsq > Q4plotsqdata

gfortran -O3 -fdefault-real-8 Q5v.f90 -llapack -o Q5v
./Q5v

gfortran -O3 -fdefault-real-8 Q5plot.f90 -llapack -o Q5plot
./Q5plot > Q5plotdata

gfortran -O3 -fdefault-real-8 Q5plotsq.f90 -llapack -o Q5plotsq
./Q5plotsq > Q5plotsqdata


chmod u+x plotme.gpl
./plotme.gpl
