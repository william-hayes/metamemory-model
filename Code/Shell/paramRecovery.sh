#!/bin/bash
#SBATCH --job-name=BIMMC_Recovery
#SBATCH --output=BIMMC_Recovery.txt
#SBATCH --error=BIMMC_Recovery.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mail-user=whayes2@binghamton.edu
#SBATCH --mail-type=ALL

module load gnu14/14.2.0
module load R/4.4.1

Rscript param_recovery.R