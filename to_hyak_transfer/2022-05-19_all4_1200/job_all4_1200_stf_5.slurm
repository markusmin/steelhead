#!/bin/bash
## JAGS run - 1200 fish, all4 (temp, flow, origin, rear) simulation

## Job name
#SBATCH --job-name=1200_all4_hyak_5

## Partition and Allocation
#SBATCH -p stf
#SBATCH -A stf

## Resources
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --ntasks=3
#SBATCH --mem=100G

## Specify the working directory for this job
#SBATCH --chdir=.

## Import any modules here
module load contrib/R-4.2.0/4.2.0

##  Scripts to be executed here
Rscript all4_1200_hyak_5.R
