#!/bin/bash
## JAGS run - 1200 fish, temperature only simulation

## Job name
#SBATCH --job-name=complete_det_hist

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
Rscript 03_hyak_complete_detection_histories_v3_part2.R
