#!/bin/bash

# Job name:
#$ -N metropolis

# Output log file:
#$ -o $JOB_NAME-$JOB_ID.log

# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory /home/username
#$ -cwd
#$ -o output-$JOB_ID.log
#$ -e error-$JOB_ID.err

##$ -m e
##$ -M  lltorreh4@alumnes.ub.edu


# Now comes the commands to be executed
# Copy files to the local disk on the node
cp input.txt $TMPDIR/
cp metropolis.out $TMPDIR/
cp observables.out $TMPDIR/

# Change to the execution directory
cd $TMPDIR/
# And run the exe
./metropolis.out
./observables.out
# Finally, we copy back all important output to the working directory
scp -r results nodo00:$SGE_O_WORKDIR


