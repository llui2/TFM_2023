# SAMPLE
gfortran -c r1279/r1279.f90 r1279/ran2.f quantum/model.f quantum/metropolis.f
chmod +x metropolis.o model.o r1279.o ran2.o
gfortran metropolis.o model.o r1279.o ran2.o -o metropolis.out
# OBSERVABLE
gfortran -c r1279/r1279.f90 r1279/ran2.f quantum/model.f quantum/observables.f
chmod +x observables.o model.o r1279.o ran2.o
gfortran observables.o model.o r1279.o ran2.o -o observables.out
# RECONSTRUCTION 
gfortran -c -g r1279/r1279.f90 r1279/ran2.f quantum/model.f quantum/pseudolikelihood.f
chmod +x pseudolikelihood.o model.o r1279.o ran2.o
gfortran pseudolikelihood.o model.o r1279.o ran2.o -o pseudolikelihood.out
# PERFORMANCE 
gfortran -c r1279/r1279.f90 r1279/ran2.f quantum/model.f quantum/performance.f
chmod +x performance.o model.o r1279.o ran2.o
gfortran performance.o model.o r1279.o ran2.o -o performance.out
rm *.o
rm *.mod
qsub job.sh