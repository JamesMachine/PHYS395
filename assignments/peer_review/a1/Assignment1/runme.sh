#!/bin/bash
#to run, first enter chmod u+x runme.sh
#then enter ./runme.sh

gfortran -O3 -fdefault-real-8 -o A1 Assignment1_first.f90
./A1 > DATA1
echo "10 terms done"
gfortran -O3 -fdefault-real-8 -o A1 Assignment1_second.f90
./A1 > DATA2
echo "100 terms done"
echo "Calulating uniform grid error..."
gfortran -O3 -fdefault-real-8 -o A1 Assignment1_third.f90
./A1

gfortran -O3 -fdefault-real-8 -o A1 Assignment1_fourth.f90
./A1 > DATA4
echo "10 terms with zeros"
gfortran -O3 -fdefault-real-8 -o A1 Assignment1_fifth.f90
./A1 > DATA5
echo "100 terms with zeros"
echo "Calculating zeros grid error..."
gfortran -O3 -fdefault-real-8 -o A1 Assignment1_sixth.f90
./A1

chmod u+x plot1.gpl plot2.gpl plot3.gpl plot4.gpl
./plot1.gpl

./plot2.gpl

./plot3.gpl

./plot4.gpl
echo "Plots have been saved"

echo "All Done"
