rm -r results
# generate samples
gfortran -c r1279/r1279.f90 r1279/ran2.f code/model.f code/metropolis.f
chmod +x metropolis.o model.o r1279.o ran2.o
gfortran metropolis.o model.o r1279.o ran2.o -o metropolis.out
rm *.o
./metropolis.out
# calculate observables
gfortran -c r1279/r1279.f90 r1279/ran2.f code/model.f code/observables.f
chmod +x observables.o model.o r1279.o ran2.o
gfortran observables.o model.o r1279.o ran2.o -o observables.out
rm *.o
rm *.mod
./observables.out
rm *.out
# plot observables
python3 code/mz.py
open results/*.pdf