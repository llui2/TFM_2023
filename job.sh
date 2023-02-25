#!/bin/bash
# Job name:
#$ -N Metropolis
# Output log file:
#$ -o $JOB_NAME-$JOB_ID.log
# Combining output/error messages into one file
##$ -j y
# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory /home/username
#$ -cwd
# With -V you pass the env variables, it's necessary. And the unset module is needed to remove some errors
#$ -V
# Uncomment the following line if you want to know in which host your job was executed
# echo "Running on " `hostname`

# Now comes the commands to be executed
# Copy files to the local disk on the node
cp input.txt $TMPDIR/
cp execute.sh $TMPDIR/
cp -r code $TMPDIR/
cp -r r1279 $TMPDIR/

# Change to the execution directory
cd $TMPDIR/
# And run the exe
./execute.sh
# Finally, we copy back all important output to the working directory
scp -r results nodo00:$SGE_O_WORKDIR


