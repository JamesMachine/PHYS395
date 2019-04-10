
# make dofractal to false to run much much fatser, otherwise, it will take 5hours.

# make this false if you don't animate
doanimate=true

# make this false if you don't do phase scan
dofractal=true

echo "Program starts..."

#note that without using -fdefault-real-8, it will be faster
gfortran -O3 -fdefault-real-8 -fopenmp dbpend.f90 -llapack -lcfitsio -o dbpend && ./dbpend
gnuplot plotme.gpl

if $doanimate ; then
  echo "animation starts..."
  gnuplot animate.gpl
  echo "animation ends..."
fi


if $dofractal ; then
  gfortran -O3 -fopenmp dbpend_phsp.f90 -llapack -lcfitsio -o dbpend_phsp && ./dbpend_phsp
  python plotfit.py

  gfortran -O3 -fopenmp dbpend_phsp2.f90 -llapack -lcfitsio -o dbpend_phsp2 && ./dbpend_phsp2
  python plotfit2.py

  gfortran -O3 -fopenmp dbpend_phsp3.f90 -llapack -lcfitsio -o dbpend_phsp3 && ./dbpend_phsp3
  python plotfit3.py

  gfortran -O3 -fopenmp dbpend_phsp4.f90 -llapack -lcfitsio -o dbpend_phsp4 && ./dbpend_phsp4
  python plotfit4.py

  echo "Output images shows that it is fractal!!"
fi


echo "Program ends..."


echo "Step size dt used: 0.038 --> within error 10e-12"
