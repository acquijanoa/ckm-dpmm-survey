#!/bin/bash
#SBATCH --job-name=test_sw
#SBATCH --time=00:05:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/test_sandwich_%j.out
#SBATCH --error=logs/test_sandwich_%j.err

module load r/4.5.0
Rscript src/R/test_sandwich_math.R
