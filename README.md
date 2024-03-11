# Inverse inference of a quantum spin glass

## About

`classic` contains all the scripts for the classical model.

`quantum` contains all the scripts for the quantum model.

`r1279` contains the random number generator scripts.

`plot` contains the scripts to plot the results.

### Execute in console
In order to run the simulations in the console, follow these steps: Set the desired parameters in the `input.txt` file. Execute `execute1.sh` to obtain the classical results or `execute2.sh` to obtain the quantum results. After running the script, a new directory called `results` will be created. Inside the `results` directory, you will find the samples stored as `S_p_SEED.bin` in a subdirectory named `sample/Txxx_Γxxx/`, where xxx corresponds to the specific parameter values used, the relative reconstruction errors stored as `g_p_SEED.bin` in a subdirectory named `accuracy/Txxx_Γxxx/`, and a subdirectory named `data/` containing the longitudinal magnetization, the transverse magnetization, and the averaged relative reconstruction error as a function of (T,Γ,p).