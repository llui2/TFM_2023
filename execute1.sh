rm -r results
# SAMPLE
gfortran -c r1279/r1279.f90 r1279/ran2.f classic/model.f classic/metropolis.f
chmod +x metropolis.o model.o r1279.o ran2.o
gfortran metropolis.o model.o r1279.o ran2.o -o metropolis.out
rm *.o
rm *.mod
./metropolis.out
rm *.out
# OBSERVABLE
gfortran -c r1279/r1279.f90 r1279/ran2.f classic/model.f classic/observables.f
chmod +x observables.o model.o r1279.o ran2.o
gfortran observables.o model.o r1279.o ran2.o -o observables.out
rm *.o
rm *.mod
./observables.out
rm *.out
# PLOT OBSERVABLE
python3 plots/mz.py
open results/mz_*.pdf
# RECONSTRUCTION 
gfortran -c -g r1279/r1279.f90 r1279/ran2.f classic/model.f classic/pseudolikelihood.f
chmod +x pseudolikelihood.o model.o r1279.o ran2.o
gfortran pseudolikelihood.o model.o r1279.o ran2.o -o pseudolikelihood.out
rm *.o
rm *.mod
./pseudolikelihood.out
rm *.out
# PERFORMANCE 
gfortran -c r1279/r1279.f90 r1279/ran2.f classic/model.f classic/performance.f
chmod +x performance.o model.o r1279.o ran2.o
gfortran performance.o model.o r1279.o ran2.o -o performance.out
rm *.o
rm *.mod
./performance.out
rm *.out
#PLOT PERFORMANCE
python3 plots/g.py
open results/g_*.pdf