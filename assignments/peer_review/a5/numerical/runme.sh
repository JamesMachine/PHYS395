# Things may not work right out of the box, as bash might not have
# permissions to run this file. The terminal will complain something
# like "./runme.sh: permission denied". If that is the case, please
# enter:
# chmod u+x runme.sh
# in the terminal while in the folder containing runme.sh

########################################################################
# QUESTION #1
########################################################################

gfortran -O3 -fdefault-real-8 wrong_harmonic_shooter.f90 -o q1.out
./q1.out

# even and odd solution for non eigenvalue E
gnuplot -persist <<-EOFMarker
set terminal pdf
set key on
set xrange [0:3]
set title "Non-Eigenvalue Solution (E=1.0)"
set output "Non-Eigenvalue Solution (E=1.0).pdf"
plot 'fort.1' index 0 u 1:2 w l title "Even",\
	'fort.1' index 1 u 1:2 w l title "Odd"
EOFMarker

########################################################################
# QUESTION #2
########################################################################

gfortran -O3 -fdefault-real-8 harmonic_shooter.f90 -o q2.out

./q2.out 0.4 1.6
cat fort.4 > Eigenvalues_Shooting.txt
cat fort.1 > Harmonic_Shooting.txt

./q2.out 2.4 3.6
cat fort.4 >> Eigenvalues_Shooting.txt
cat fort.1 >> Harmonic_Shooting.txt

./q2.out 4.4 5.6
cat fort.4 >> Eigenvalues_Shooting.txt
cat fort.1 >> Harmonic_Shooting.txt

./q2.out 6.4 7.6
cat fort.4 >> Eigenvalues_Shooting.txt
cat fort.1 >> Harmonic_Shooting.txt

./q2.out 8.4 9.6
cat fort.4 >> Eigenvalues_Shooting.txt
cat fort.1 >> Harmonic_Shooting.txt

# plot wave functions for eigenvalues
gnuplot -persist <<-EOFMarker
set terminal pdf
set key on
set xrange [0:7]
set title "Harmonic Wave Function (shooting)"
set output "Harmonic Wave Functions (shooting).pdf"
plot 'Harmonic_Shooting.txt' index 0 u 1:2 w l title "0.5",\
	'Harmonic_Shooting.txt' index 1 u 1:2 w l title "1.5",\
	'Harmonic_Shooting.txt' index 2 u 1:2 w l title "2.5",\
	'Harmonic_Shooting.txt' index 3 u 1:2 w l title "3.5",\
	'Harmonic_Shooting.txt' index 4 u 1:2 w l title "4.5",\
	'Harmonic_Shooting.txt' index 5 u 1:2 w l title "5.5",\
	'Harmonic_Shooting.txt' index 6 u 1:2 w l title "6.5",\
	'Harmonic_Shooting.txt' index 7 u 1:2 w l title "7.5",\
	'Harmonic_Shooting.txt' index 8 u 1:2 w l title "8.5",\
	'Harmonic_Shooting.txt' index 9 u 1:2 w l title "9.5"
EOFMarker

# plot normalized probability densities for eigenvalues
gnuplot -persist <<-EOFMarker
set terminal pdf
set key on
set xrange [0:7]
set title "Harmonic Probability Densities (shooting)"
set output "Harmonic Probability Densities (shooting).pdf"
plot 'Harmonic_Shooting.txt' index 0 u 1:3 w l title "0.5",\
	'Harmonic_Shooting.txt' index 1 u 1:3 w l title "1.5",\
	'Harmonic_Shooting.txt' index 2 u 1:3 w l title "2.5",\
	'Harmonic_Shooting.txt' index 3 u 1:3 w l title "3.5",\
	'Harmonic_Shooting.txt' index 4 u 1:3 w l title "4.5",\
	'Harmonic_Shooting.txt' index 5 u 1:3 w l title "5.5",\
	'Harmonic_Shooting.txt' index 6 u 1:3 w l title "6.5",\
	'Harmonic_Shooting.txt' index 7 u 1:3 w l title "7.5",\
	'Harmonic_Shooting.txt' index 8 u 1:3 w l title "8.5",\
	'Harmonic_Shooting.txt' index 9 u 1:3 w l title "9.5"
EOFMarker

########################################################################
# question #3
########################################################################

gfortran -O3 -fdefault-real-8 anharmonic_shooter.f90 -o q3.out

./q3.out 0.3 1.6
cat fort.4 > Eigenvalues_Anharmonic_Shooting.txt
cat fort.1 > Anharmonic_Shooting.txt

./q3.out 2.9 4.7
cat fort.4 >> Eigenvalues_Anharmonic_Shooting.txt
cat fort.1 >> Anharmonic_Shooting.txt

./q3.out 6.4 8.5
cat fort.4 >> Eigenvalues_Anharmonic_Shooting.txt
cat fort.1 >> Anharmonic_Shooting.txt

./q3.out 10.4 12.8
cat fort.4 >> Eigenvalues_Anharmonic_Shooting.txt
cat fort.1 >> Anharmonic_Shooting.txt

./q3.out 14.9 17.6
cat fort.4 >> Eigenvalues_Anharmonic_Shooting.txt
cat fort.1 >> Anharmonic_Shooting.txt

# plot wave functions for eigenvalues
gnuplot -persist <<-EOFMarker
set terminal pdf
set key on
set xrange [0:5]
set title "Anharmonic Wave Functions (shooting)"
set output "Anharmonic Wave Functions (shooting).pdf"
plot 'Anharmonic_Shooting.txt' index 0 u 1:2 w l title "0.42",\
	'Anharmonic_Shooting.txt' index 1 u 1:2 w l title "1.51",\
	'Anharmonic_Shooting.txt' index 2 u 1:2 w l title "2.96",\
	'Anharmonic_Shooting.txt' index 3 u 1:2 w l title "4.62",\
	'Anharmonic_Shooting.txt' index 4 u 1:2 w l title "6.45",\
	'Anharmonic_Shooting.txt' index 5 u 1:2 w l title "8.43",\
	'Anharmonic_Shooting.txt' index 6 u 1:2 w l title "10.53",\
	'Anharmonic_Shooting.txt' index 7 u 1:2 w l title "12.74",\
	'Anharmonic_Shooting.txt' index 8 u 1:2 w l title "15.05",\
	'Anharmonic_Shooting.txt' index 9 u 1:2 w l title "17.45"
EOFMarker

# plot normalized probability densities for eigenvalues
gnuplot -persist <<-EOFMarker
set terminal pdf
set key on
set xrange [0:5]
set title "Anharmonic Probability Densities (shooting)"
set output "Anharmonic Probability Densities (shooting).pdf"
plot 'Anharmonic_Shooting.txt' index 0 u 1:3 w l title "0.42",\
	'Anharmonic_Shooting.txt' index 1 u 1:3 w l title "1.51",\
	'Anharmonic_Shooting.txt' index 2 u 1:3 w l title "2.96",\
	'Anharmonic_Shooting.txt' index 3 u 1:3 w l title "4.62",\
	'Anharmonic_Shooting.txt' index 4 u 1:3 w l title "6.45",\
	'Anharmonic_Shooting.txt' index 5 u 1:3 w l title "8.43",\
	'Anharmonic_Shooting.txt' index 6 u 1:3 w l title "10.53",\
	'Anharmonic_Shooting.txt' index 7 u 1:3 w l title "12.74",\
	'Anharmonic_Shooting.txt' index 8 u 1:3 w l title "15.05",\
	'Anharmonic_Shooting.txt' index 9 u 1:3 w l title "17.45"
EOFMarker

########################################################################
# QUESTION #4
########################################################################
gfortran -O3 -fdefault-real-8 harmonic_rayleigh.f90 -llapack -o q4.out
./q4.out
cat fort.4 > Eigenvalues_dgeev.txt
cat fort.1 > Harmonic_dgeev.txt

# plot wave functions for dgeev solution to harmonic oscillator
gnuplot -persist <<-EOFMarker
set terminal pdf
set key on
set xrange [0:7]
set title "Harmonic Wave Functions (dgeev)"
set output "Harmonic Wave Functions (dgeev).pdf"
plot 'Harmonic_dgeev.txt' index 0 u 1:2 w l title "0.5",\
	'Harmonic_dgeev.txt' index 1 u 1:2 w l title "1.5",\
	'Harmonic_dgeev.txt' index 2 u 1:2 w l title "2.5",\
	'Harmonic_dgeev.txt' index 3 u 1:2 w l title "3.5",\
	'Harmonic_dgeev.txt' index 4 u 1:2 w l title "4.5",\
	'Harmonic_dgeev.txt' index 5 u 1:2 w l title "5.5",\
	'Harmonic_dgeev.txt' index 6 u 1:2 w l title "6.5",\
	'Harmonic_dgeev.txt' index 7 u 1:2 w l title "7.5",\
	'Harmonic_dgeev.txt' index 8 u 1:2 w l title "8.5",\
	'Harmonic_dgeev.txt' index 9 u 1:2 w l title "9.5"
EOFMarker

# plot normalized probability densities for dgeev solution to harmonic oscillator
gnuplot -persist <<-EOFMarker
set terminal pdf
set key on
set xrange [0:7]
set title "Harmonic Probability Densities (dgeev)"
set output "Harmonic Probability Densities (dgeev).pdf"
plot 'Harmonic_dgeev.txt' index 0 u 1:3 w l title "0.5",\
	'Harmonic_dgeev.txt' index 1 u 1:3 w l title "1.5",\
	'Harmonic_dgeev.txt' index 2 u 1:3 w l title "2.5",\
	'Harmonic_dgeev.txt' index 3 u 1:3 w l title "3.5",\
	'Harmonic_dgeev.txt' index 4 u 1:3 w l title "4.5",\
	'Harmonic_dgeev.txt' index 5 u 1:3 w l title "5.5",\
	'Harmonic_dgeev.txt' index 6 u 1:3 w l title "6.5",\
	'Harmonic_dgeev.txt' index 7 u 1:3 w l title "7.5",\
	'Harmonic_dgeev.txt' index 8 u 1:3 w l title "8.5",\
	'Harmonic_dgeev.txt' index 9 u 1:3 w l title "9.5"
EOFMarker

########################################################################
# QUESTION #5
########################################################################

gfortran -O3 -fdefault-real-8 anharmonic_rayleigh.f90 -llapack -o q5.out
./q5.out

cat fort.4 > Eigenvalues_Anharmonic_dgeev.txt
cat fort.1 > Anharmonic_dgeev.txt

# plot wave functions for dgeev solution to harmonic oscillator
gnuplot -persist <<-EOFMarker
set terminal pdf
set key on
set xrange [0:5]
set title "Anharmonic Wave Functions (dgeev)"
set output "Anharmonic Wave Functions (dgeev).pdf"
plot 'Anharmonic_dgeev.txt' index 0 u 1:2 w l title "0.42",\
	'Anharmonic_dgeev.txt' index 1 u 1:2 w l title "1.51",\
	'Anharmonic_dgeev.txt' index 2 u 1:2 w l title "2.96",\
	'Anharmonic_dgeev.txt' index 3 u 1:2 w l title "4.62",\
	'Anharmonic_dgeev.txt' index 4 u 1:2 w l title "6.45",\
	'Anharmonic_dgeev.txt' index 5 u 1:2 w l title "8.43",\
	'Anharmonic_dgeev.txt' index 6 u 1:2 w l title "10.53",\
	'Anharmonic_dgeev.txt' index 7 u 1:2 w l title "12.74",\
	'Anharmonic_dgeev.txt' index 8 u 1:2 w l title "15.05",\
	'Anharmonic_dgeev.txt' index 9 u 1:2 w l title "17.45"
EOFMarker

# plot normalized probability densities for dgeev solution to harmonic oscillator
gnuplot -persist <<-EOFMarker
set terminal pdf
set key on
set xrange [0:5]
set title "Anharmonic Probability Densities (dgeev)"
set output "Anharmonic Probability Densities (dgeev).pdf"
plot 'Anharmonic_dgeev.txt' index 0 u 1:3 w l title "0.42",\
	'Anharmonic_dgeev.txt' index 1 u 1:3 w l title "1.51",\
	'Anharmonic_dgeev.txt' index 2 u 1:3 w l title "2.96",\
	'Anharmonic_dgeev.txt' index 3 u 1:3 w l title "4.62",\
	'Anharmonic_dgeev.txt' index 4 u 1:3 w l title "6.45",\
	'Anharmonic_dgeev.txt' index 5 u 1:3 w l title "8.43",\
	'Anharmonic_dgeev.txt' index 6 u 1:3 w l title "10.53",\
	'Anharmonic_dgeev.txt' index 7 u 1:3 w l title "12.74",\
	'Anharmonic_dgeev.txt' index 8 u 1:3 w l title "15.05",\
	'Anharmonic_dgeev.txt' index 9 u 1:3 w l title "17.45"
EOFMarker

# put everything away in a nice folder
mkdir data
mv *.pdf ./data
mv *.txt ./data

echo NOTE:
echo Some of the dgeev wave functions are upside down, but it\'s a valid solution so I don\'t care
echo Also the dgeev wave functions are normalized, which I don\'t bother to do in the shooting method
