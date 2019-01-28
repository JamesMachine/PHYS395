
echo "program start!..."
gfortran -fdefault-real-8 q1.f90 -o q1 && ./q1 > output1 && cat output1
gfortran -fdefault-real-8 q2.f90 -o q2 && ./q2
gfortran -fdefault-real-8 q2_2.f90 -o q2_2 && ./q2_2


echo "plotting start!..."
gnuplot q2.gpl
echo "plotting end! Saved into pdf files..."

echo "program end!..."
