#!/bin/bash
# Generate the finite population for one clustering strength.
#   sbatch --export=ALL,CLUSTER_C=0.3 run_generate_population.sbatch

#SBATCH --job-name=gen_pop
#SBATCH --time=00:15:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/gen_pop_%j.out
#SBATCH --error=logs/gen_pop_%j.err

export CLUSTER_C="${CLUSTER_C:-0.3}"
export POP_SCENARIO="${POP_SCENARIO:-base}"
mkdir -p logs

module load r/4.5.0

Rscript src/R/generate_population.R
