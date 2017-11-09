#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1
#SBATCH --ntasks=1
R CMD BATCH --no-save parallel_eval.R eval.out 
