#!/bin/bash
## stan actual run, intercept, origin, and rear


## Job name
#SBATCH --job-name=yak-ff-cov
#SBATCH --mail-type=ALL

## Partition and Allocation
#SBATCH -p ckpt
#SBATCH -A stf

## Resources
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --ntasks=28
#SBATCH --mem=120G

## Specify the working directory for this job
#SBATCH --chdir=.

## We don't need any modules beyond the default loaded ones

##  Scripts to be executed here
apptainer exec --bind ../../:/gscratch/scrubbed/mmin ../../../../../gscratch/stf/mmin/cmdstanr_latest.sif Rscript 14.7-yak-ff-cov-ckpt.R
