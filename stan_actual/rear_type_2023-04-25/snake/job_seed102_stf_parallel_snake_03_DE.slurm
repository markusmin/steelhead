#!/bin/bash
## stan actual run, intercept only


## Job name
#SBATCH --job-name=102_snake

## Partition and Allocation
#SBATCH -p stf
#SBATCH -A stf

## Resources
#SBATCH --nodes=1
#SBATCH --time=144:00:00
#SBATCH --ntasks=28
#SBATCH --mem=120G

## Specify the working directory for this job
#SBATCH --chdir=.

## Import any modules here
module load contrib/gcc/6.2.0 
module load contrib/R-4.2.0/4.2.0

##  Scripts to be executed here
Rscript seed102_parallel_snake_03_stan_actual_int_origin_mox_deteff.R
