echo "compiling..."
gfortran -c -fdefault-real-8 newtonMethod.f90 -llapack
gfortran -c -fdefault-real-8 bracketSearch.f90 -llapack
gfortran -c -fdefault-real-8 LMFit.f90 -llapack
gfortran -c -fdefault-real-8 tester.f90 -llapack

echo "running..."
gfortran -g -llapack tester.o newtonMethod.o bracketSearch.o LMFit.o -o tester &&./tester > output && cat output

echo "testing code complete"
echo "plotting..."
gnuplot Plot.gpl

echo "all actions complete"
echo "outputs can be found at 'Output' and 'Plot.pdf'"
