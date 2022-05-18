#!/bin/bash
## JAGS run - 1200 fish, origin only simulation

## Job name
#SBATCH --job-name=1200_origin_hyak

## Partition and Allocation
#SBATCH -p ckpt
#SBATCH -A stf-ckpt

## Resources
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --ntasks=15
#SBATCH --mem=100G

## Specify the working directory for this job
#SBATCH --chdir=.

## Import any modules here
module load contrib/R-4.2.0/4.2.0

##  Scripts to be executed here
## Note: We are using 3 tasks (CPUs) for each script, since we are running 3 chains in parallel
srun --ntasks 3 Rscript origin_rear_1200_hyak_1.R
srun --ntasks 3 Rscript origin_rear_1200_hyak_2.R
srun --ntasks 3 Rscript origin_rear_1200_hyak_3.R
srun --ntasks 3 Rscript origin_rear_1200_hyak_4.R
srun --ntasks 3 Rscript origin_rear_1200_hyak_5.R
