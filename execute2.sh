# reconstruction 
gfortran -c r1279/r1279.f90 r1279/ran2.f code/model.f code/pseudolikelihood.f
chmod +x pseudolikelihood.o model.o r1279.o ran2.o
gfortran pseudolikelihood.o model.o r1279.o ran2.o -o pseudolikelihood.out
rm *.o
rm *.mod
./pseudolikelihood.out