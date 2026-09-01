#!/bin/bash
#SBATCH --job-name=BIMMC_Fit
#SBATCH --chdir=/data/home/whayes2/metamemory-model/Code/Shell
#SBATCH --output=/data/home/whayes2/metamemory-model/Code/Shell/BIMMC_Fit.out
#SBATCH --error=/data/home/whayes2/metamemory-model/Code/Shell/BIMMC_Fit.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=whayes2@binghamton.edu
#SBATCH --mail-type=ALL

module load gnu14/14.2.0
module load R/4.5.0

cd /data/home/whayes2/metamemory-model/Code/Model/

Rscript model_fitting.R
