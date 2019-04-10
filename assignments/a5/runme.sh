#Run everything with command:
#time bash runme.sh


echo "Program starts..."
echo ""

#delete output folder by rm -r output
#create output folder if they don't exist
mkdir -p output
mkdir -p ./output/output_q1
mkdir -p ./output/output_q2
mkdir -p ./output/output_q3
mkdir -p ./output/output_q4
mkdir -p ./output/output_q5

gfortran -O3 -fdefault-real-8 shrodinger_shooting_q1.f90 -llapack -o shrodinger_shooting_q1 && ./shrodinger_shooting_q1
gfortran -O3 -fdefault-real-8 shrodinger_shooting_q2.f90 -llapack -o shrodinger_shooting_q2 && ./shrodinger_shooting_q2
gfortran -O3 -fdefault-real-8 shrodinger_shooting_q3.f90 -llapack -o shrodinger_shooting_q3 && ./shrodinger_shooting_q3
gfortran -O3 -fdefault-real-8 shrodinger_dgeev_q4.f90 -llapack -o shrodinger_dgeev_q4 && ./shrodinger_dgeev_q4
gfortran -O3 -fdefault-real-8 shrodinger_dgeev_q5.f90 -llapack -o shrodinger_dgeev_q5 && ./shrodinger_dgeev_q5

echo "Plotting starts..."
#python plotme.py #python takes 10x slower but easy string manipulation
gnuplot plotme.gpl #gnuplot much faster but harder to code than python
echo "Plotting ends..."

echo ""
echo""

echo "Comment:"
echo "The psi (VR) for Q4 and Q5 are already normalized by dgeev. For all psi from \
Q1 to Q5 are all normalized properly."
echo "Probability amplitude (psi amplitude) for Q2,Q3 and Q4,Q5 are different because \
they have different x range of numerical computation (Note that you cannot compute till infinity). \
These different range of computation will lead to the different results. However, they are \
all normalized to 1, so basically same answer."
echo "I can also use integrator for Q4 and Q5 by found eigenvalues using dgeev, but why not just using dgeev!"
echo "Finally, my performnace is really fast huh? : D"
echo ""

echo "Program ends..."
