#!/bin/bash
## stan actual run, intercept only


## Job name
#SBATCH --job-name=stan_actual_int

## Partition and Allocation
#SBATCH -p stf
#SBATCH -A stf

## Resources
#SBATCH --nodes=3
#SBATCH --time=24:00:00
#SBATCH --ntasks=3
#SBATCH --mem=100G

## Specify the working directory for this job
#SBATCH --chdir=.

## Import any modules here
module load contrib/gcc/6.2.0 
module load contrib/R-4.2.0/4.2.0

##  Scripts to be executed here
Rscript 01_stan_actual_int_only_mox_500iter.R
