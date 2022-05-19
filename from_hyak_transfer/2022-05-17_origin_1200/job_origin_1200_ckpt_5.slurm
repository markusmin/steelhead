#!/bin/bash
## JAGS run - 1200 fish, origin only simulation

## Job name
#SBATCH --job-name=1200_origin_hyak_5

## Partition and Allocation
#SBATCH -p ckpt
#SBATCH -A stf-ckpt

## Resources
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --ntasks=3
#SBATCH --mem=100G

## Specify the working directory for this job
#SBATCH --chdir=.

## Import any modules here
module load contrib/R-4.2.0/4.2.0

##  Scripts to be executed here
Rscript origin_1200_hyak_5.R