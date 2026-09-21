#!/bin/bash
#SBATCH --job-name=deff_eval
#SBATCH --time=00:10:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/deff_eval_%j.out
#SBATCH --error=logs/deff_eval_%j.err

mkdir -p logs

module load r/4.5.0

Rscript src/R/estimate_deff.R
