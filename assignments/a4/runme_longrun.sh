
echo "Long run program for Q3 starts..."

gfortran -O3 -fopenmp dbpend_longrun.f90 -llapack -lcfitsio -o dbpend_longrun && ./dbpend_longrun

echo "plotting starts..."
python plotfit_longrun.py
echo "plotting ends..."

echo "Long run program program ends..."
