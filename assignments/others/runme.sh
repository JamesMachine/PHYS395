echo "Program starts..."



#gfortran -O3 -fdefault-real-8 test.f90 -o test && ./test

#gfortran -O3 -fdefault-real-8 my_diffuse.f90 -o my_diffuse && ./my_diffuse
#gfortran -O3 -fdefault-real-8 my_wave.f90 -o my_wave && ./my_wave
#gfortran -O3 -fdefault-real-8 my_fft.f90 -lfftw3 -o my_fft &&./my_fft > ./output/my_fft.dat
#gfortran -O3 -fdefault-real-8 fft_bvp1d.f90 -lfftw3 -o fft_bvp1d &&./fft_bvp1d > ./output/fft_bvp1d.dat



#gnuplot animate_diffuse.gpl
#gnuplot animate_wave.gpl
gnuplot plotme.gpl

echo "Program ends..."
