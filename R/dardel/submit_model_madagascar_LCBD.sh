#!/bin/bash -l
#SBATCH -A naiss2025-22-41
#SBATCH -p shared
#SBATCH -c 1
#SBATCH -t 05:00:00
#SBATCH --mem 20000
#SBATCH -e logfiles/error_madagascar_LCBD.e
#SBATCH --mail-type=FAIL
#SBATCH -J fit_madagascar_LCBD

module load PDC/23.12    
module load R/4.4.0
Rscript fit_model_madagascar_LCBD.R > logfiles/out_madagascar_LCBD.out
