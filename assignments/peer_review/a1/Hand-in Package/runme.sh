echo "compiling..."
gfortran -c -fdefault-real-8 gaussjordan.f90
gfortran -c -fdefault-real-8 polyapprox.f90
gfortran -c -fdefault-real-8 tester.f90

echo "running..."
gfortran tester.o gaussjordan.o polyapprox.o -o tester && ./tester > output && cat output

echo "plotting..."
gnuplot prob2.gpl
gnuplot prob3.gpl
gnuplot prob5.gpl
echo "plots complete, saved as PDFs in directory."
