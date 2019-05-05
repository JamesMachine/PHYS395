
echo "program start!..."

mkdir -p output
#mkdir -p ./output/output_q1
mkdir -p ./output/output_q2
mkdir -p ./output/output_q3
mkdir -p ./output/output_q4
mkdir -p ./output/output_q5


gfortran -fdefault-real-8 q1.f90 -llapack -o q1 && ./q1
gfortran -fdefault-real-8 q2.f90 -llapack -o q2 && ./q2
gfortran -fdefault-real-8 q2_2.f90 -llapack -o q2_2 && ./q2_2
gfortran -fdefault-real-8 q3.f90 -llapack -o q3 && ./q3
gfortran -fdefault-real-8 q3_2.f90 -llapack -o q3_2 && ./q3_2
gfortran -fdefault-real-8 q4.f90 -llapack -o q4 && ./q4
gfortran -fdefault-real-8 q4_2.f90 -llapack -o q4_2 && ./q4_2
gfortran -fdefault-real-8 q5.f90 -llapack -o q5 && ./q5
gfortran -fdefault-real-8 q5_2.f90 -llapack -o q5_2 && ./q5_2



echo "plotting start!..."
gnuplot q2.gpl
gnuplot q3.gpl
gnuplot q4.gpl
gnuplot q5.gpl
echo "plotting end! Saved into pdf files..."

echo "program end!..."
