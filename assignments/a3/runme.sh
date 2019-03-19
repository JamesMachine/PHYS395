echo 'Program started!...'


# compile codes and run for Q1 - Q3
gfortran -g -O3 -fdefault-real-8 a3_newton.f90 -llapack -o a3_newton && ./a3_newton
gfortran -g -O3 -fdefault-real-8 a3_bracket.f90 -llapack -o a3_bracket && ./a3_bracket


#reading the # of lines of the data --> put 9130 in fortran code for the # of measurement points
array=($(wc data.dat))
lines=${array[0]}

#for Q4 and Q5, compile and run
gfortran -g -O3 -fdefault-real-8 a3_nonlinear_fit_LM.f90 -llapack -o a3_nonlinear_fit_LM && ./a3_nonlinear_fit_LM

# showing results
echo "Q1: roots (x,y) are:"
cat output1_roots

echo "Q2,Q3: Minima (x,y) are:"
cat output2_minima

echo "Q5: Showing chi results and coefficients:"
echo "# of data points in data.dat: " $lines
cat output3_parameterResults


echo "Note that graphs for Q1 to Q5 are all plotted into pdf!"
gnuplot plotme.gpl

echo "Program ends!..."


# Other methods are following (LM method works and the best..., other methods under development)
#cd GausseNewtonAndGradDec_underDev
#gfortran -g -O3 -fdefault-real-8 a3_nonlinear_fit_gradDec.f90 -llapack -o a3_nonlinear_fit_gradDec && ./a3_nonlinear_fit_gradDec
#gfortran -g -O3 -fdefault-real-8 a3_nonlinear_fit_gaussNewton.f90 -llapack -o a3_nonlinear_fit_gaussNewton && ./a3_nonlinear_fit_gaussNewton
