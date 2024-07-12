#!/bin/bash --login
#SBATCH --job-name=long_job
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=02:00:00
##SBATCH --gres=gpu:1
#SBATCH --mem=32
#SBATCH --export=NONE

# To compile with the GNU toolchain
module load gcc/9.4.0
module load matlab
# leave in, it lists the environment loaded by the modules
module list

#  Note: SLURM_JOBID is a unique number for every job.
#  These are generic variables
EXECUTABLE=test1.m
SCRATCH=$MYSCRATCH/test1.m/$SLURM_JOBID
RESULTS=$MYGROUP/test1_results/$SLURM_JOBID

###############################################
# Creates a unique directory in the SCRATCH directory for this job to run in.
if [ ! -d $SCRATCH ]; then 
    mkdir -p $SCRATCH 
fi 
echo SCRATCH is $SCRATCH

###############################################
# Creates a unique directory in your GROUP directory for the results of this job
if [ ! -d $RESULTS ]; then 
    mkdir -p $RESULTS 
fi
echo the results directory is $RESULTS

################################################
# declare the name of the output file or log file
mkdir -p output_folder
OUTPUT=$output_folder

#############################################
#   Copy input files to $SCRATCH
#   then change directory to $SCRATCH

cp $EXECUTABLE $SCRATCH
cp KnudsenDiffusivity.m $SCRATCH
cp convectivemasstransfercoefficient.m $SCRATCH
cp binarydiffusivity.m $SCRATCH
cp ReactionRateConstant.m $SCRATCH

cd $SCRATCH

srun -n 1 ./$EXECUTABLE >> ${OUTPUT}

#############################################
#    $OUTPUT file to the unique results dir
# note this can be a copy or move  
mv  $OUTPUT ${RESULTS}

cd $HOME

###########################
# Clean up $SCRATCH 

rm -r $SCRATCH

echo test1_jxu job finished at  `date`


