
echo "Long run program for Q3 starts..."

gfortran -O3 -fopenmp dbpend_phsp_longrun.f90 -llapack -lcfitsio -o dbpend_phsp_longrun && ./dbpend_phsp_longrun

echo "plotting starts..."
python plotfit_longrun.py
echo "plotting ends..."

echo "Long run program program ends..."
