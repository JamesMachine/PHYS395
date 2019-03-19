
# make this false if you don't animate to run faster
doanimate=true

echo "Program starts..."

#note that without using -fdefault-real-8, it will be faster
gfortran -O3 -fdefault-real-8 -fopenmp dbpend.f90 -llapack -lcfitsio -o dbpend && ./dbpend
gfortran -O3 -fopenmp dbpend_phsp.f90 -llapack -lcfitsio -o dbpend_phsp && ./dbpend_phsp

echo "plotting starts..."
gnuplot plotme.gpl
python plotfit.py
echo "plotting ends..."


if $doanimate ; then
  echo "animation starts..."
  gnuplot animate.gpl
  echo "animation ends..."
fi


echo "Program ends..."


echo "Step size dt used: 0.038 --> within error 10e-12"
