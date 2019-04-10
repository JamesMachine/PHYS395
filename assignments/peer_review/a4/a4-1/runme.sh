#!/bin/bash
#to run, first enter the following:
#chmod u+x runme.sh animating.gpl save.py savezoom1.py savezoom2.py savezoom3.py savezoom4.py
#then enter ./runme.sh

echo 'This may take some time, plots are included in archive if you do not want to wait'
echo 'Run animating.gpl to just see the double pendulum animation'

echo "Exploring Phase Space..."
gfortran -O3 -fdefault-real-8 -fopenmp -o A4 A4Q3.f90 -lcfitsio
./A4
echo "Done, saving image..."
./save.py
echo "Zooming in..."
gfortran -O3 -fdefault-real-8 -fopenmp -o A4 A4Q3Zoom1.f90 -lcfitsio
./A4
echo "Done, saving image..."
./savezoom1.py
echo "Zooming in..."
gfortran -O3 -fdefault-real-8 -fopenmp -o A4 A4Q3Zoom2.f90 -lcfitsio
./A4
echo "Done, saving image..."
./savezoom2.py
echo "Zooming in..."
gfortran -O3 -fdefault-real-8 -fopenmp -o A4 A4Q3Zoom3.f90 -lcfitsio
./A4
echo "Done, saving image..."
./savezoom3.py
echo "Zooming in..."
gfortran -O3 -fdefault-real-8 -fopenmp -o A4 A4Q3Zoom4.f90 -lcfitsio
./A4
echo "Done, saving image..."
./savezoom4.py
echo 'Now showing Energy loss in integration algorithm...'
./plot1.gpl
echo "Now animating pendulum..."
gfortran -O3 -fdefault-real-8 -fopenmp -o A4 A4.f90 -lcfitsio
./A4 > 'DATA'

./animating.gpl
