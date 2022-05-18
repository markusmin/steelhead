#!/bin/bash
## JAGS run - 1200 fish, with covariates simulation

## Job name
#SBATCH --job-name=1200cov_hyak

## Partition and Allocation
#SBATCH -p stf
#SBATCH -A stf

## Resources
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --ntasks=3
#SBATCH --mem=100G

## Specify the working directory for this job
#SBATCH --chdir=.

## Import any modules here


##  Scripts to be executed here
Rscript 1200cov_hyak.R