# Things may not work right out of the box, as bash might not have
# permissions to run this file. The terminal will complain something
# like "./runme.sh: permission denied". If that is the case, please
# enter:
# chmod u+x runme.sh
# in the terminal while in the folder containing runme.sh

# Run the programs that answer questions #1-4
gfortran -O3 -fdefault-real-8 LeastSquares3.f90 -llapack
./a.out < Assignment\ #2.dat
mv fort.1 n=3.txt

gnuplot -persist <<-EOFMarker
set terminal pdf
set title "Least Squares Fit (without dfelss(), n = 3)
set key
set output "n=3.pdf"
p'Assignment #2.dat' w d t 'Data', 'n=3.txt' w l t 'SVD Fit'
EOFMarker

gfortran -O3 -fdefault-real-8 LeastSquares7.f90 -llapack
./a.out < Assignment\ #2.dat
mv fort.1 n=7.txt

gnuplot -persist <<-EOFMarker
set terminal pdf
set title "Least Squares Fit (without dfelss(), n = 7)
set key
set output "n=7.pdf"
p'Assignment #2.dat' w d t 'Data', 'n=7.txt' w l t 'SVD Fit'
EOFMarker

# Question 5
gfortran -O3 -fdefault-real-8 LapackLeastSquares.f90 -llapack
./a.out < Assignment\ #2.dat
mv fort.1 dgelss.txt

gnuplot -persist <<-EOFMarker
set terminal pdf
set title "Least Squares Fit by dgelss() (n = 3)
set key
set output "dgelss().pdf"
p'Assignment #2.dat' w d t 'Data', 'dgelss.txt' w l t 'dgelss() Fit'
EOFMarker
