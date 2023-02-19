rm -r results/observables
gfortran -c r1279/r1279.f90 r1279/ran2.f code/model.f code/observables.f
chmod +x observables.o model.o r1279.o ran2.o
gfortran observables.o model.o r1279.o ran2.o -o observables.out
rm *.o
./observables.out
rm *.out