
echo "Program started..."

#delete output folder by rm -r output
#create output folder if they don't exist
mkdir -p output
mkdir -p ./output/output_q1
mkdir -p ./output/output_q2
mkdir -p ./output/output_q3
mkdir -p ./output/output_q4
mkdir -p ./output/output_q5


gfortran -g -O3 -fdefault-real-8 q1.f90 -llapack -o q1 && ./q1
gfortran -g -O3 -fdefault-real-8 q2.f90 -llapack -o q2 && ./q2
gfortran -g -O3 -fdefault-real-8 q3.f90 -llapack -o q3 && ./q3
gfortran -g -O3 -fdefault-real-8 q4.f90 -llapack -o q4 && ./q4 #--> you can run but this fails
#gfortran -g -O3 -fdefault-real-8 q5.f90 -llapack -o q5 && ./q5



gnuplot plotme.gpl

echo "Program ends..."
