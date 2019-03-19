#!/bin/bash

fbuild () {
    gfortran -llapack -fdefault-real-8 -O3 -fimplicit-none -Wall $@
}

echo "Hello. Beginning production run of Assignment 02."
echo "Building the utilities..."
fbuild -c utilities.f90
echo "Utilities built."
echo

echo "Now compiling problem 1..."
fbuild -c question1.f90
fbuild -o question1 utilities.o question1.o

echo "Compiling problem 1 done. Running problem 1..."
./question1 > question1.dat

echo "Running problem 1 finished. Answer acquired is:"
cat question1.dat

echo

echo "Now compiling problem 2 ..."
fbuild -c question2.f90
fbuild -o question2 utilities.o question2.o

echo "Compiling problem 2 done. Running problem 2..."
./question2 > question2.dat

echo "Running problem 2 finished. Answer acquired is:"
cat question2.dat

echo

echo "Now compiling problem 4..."

fbuild -c question4.f90
fbuild -o question4 utilities.o question4.o

echo "Compiling problem 4 done. Running problem 5..."
./question4 < DATA_FILES/DATA.dat > question4.dat

echo "Running problem 4 finished. Answer acquired is:"
cat question4.dat

echo "Making plots..."
python plot.py

echo

# we are finished with the executables. Delete all of them
echo "Done with executables. Removing them from EXISTENCE..."
rm question? *.o

# we are done with the data files. Delete all of them
echo "Done with the data files. Removing them from EXISTENCE..."
rm *.dat

echo "Finished executing. Have a great day!"