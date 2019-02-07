
echo 'Assignment 2 Program start!...'
echo "Note that in my source code, I set n=4 and n=8 respectively since my n starts from 1 but not 0, but here I indicate n=3 and n=7 as starting from 0"


echo 'For n=3 (4 terms) ----------------------------------------'
gfortran -g -O3 -fdefault-real-8 a2code_n3.f90 -llapack -o a2code_n3 && ./a2code_n3
echo 'This result indicates that our fitting is good.'
echo
echo

echo 'For n=7 (8 terms) ----------------------------------------'
gfortran -g -O3 -fdefault-real-8 a2code_n7.f90 -llapack -o a2code_n7 && ./a2code_n7
echo 'This result indicates that our fitting is better than the case of n=3.'
echo
echo

echo "Next, repeat for n=3 and n=7 using dgelss..."
echo 'For n=3 (4 terms) using dgelss ----------------------------------------'
gfortran -g -O3 -fdefault-real-8 a2code_n3_lss.f90 -llapack -o a2code_n3_lss && ./a2code_n3_lss
echo 'This result indicates that our fitting is good.'
echo
echo


echo 'For n=7 (4 terms) using dgelss ----------------------------------------'
gfortran -g -O3 -fdefault-real-8 a2code_n7_lss.f90 -llapack -o a2code_n7_lss && ./a2code_n7_lss
echo 'This result indicates that our fitting is better than the case of n=3.'
echo
echo "svd method agrees very closely with the method using dgelss! Note error % is (slightly) low for svd method for n=3 but low for dgelss method for n=7. In"
echo "In other words, accuracy is high uing svd method for lower order n=3 but high using dgelss for higher order n=7!"


gnuplot plotme.gpl
echo 'Program end!...'
