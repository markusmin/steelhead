#!/bin/bash
## JAGS run - full dataset, intercept only

## Job name
#SBATCH --job-name=intonly_full

## Partition and Allocation
#SBATCH -p ckpt
#SBATCH -A stf-ckpt

## Resources
#SBATCH --nodes=1
#SBATCH --time=16:00:00
#SBATCH --ntasks=3
#SBATCH --mem=100G

## Specify the working directory for this job
#SBATCH --chdir=.

## Import any modules here
module load contrib/R-4.2.0/4.2.0

##  Scripts to be executed here
Rscript 12.1_hyak_fullmodel_intercept_only_1.1_2kiter.R
