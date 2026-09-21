#!/bin/bash
# CLUSTER_C selects the clustering-strength grid point (default 0.3).
#   sbatch --export=ALL,CLUSTER_C=0.3 run_generate_population.sbatch
#   sbatch --export=ALL,CLUSTER_C=1   run_generate_population.sbatch
#
# Multi-population variance-decomposition study: submit as an array job to
# generate one independent superpopulation replicate per task (POP_ID
# defaults to $SLURM_ARRAY_TASK_ID when not set explicitly), e.g.
#   sbatch --array=1-25 --export=ALL,CLUSTER_C=0.3 run_generate_population.sbatch

#SBATCH --job-name=gen_pop
#SBATCH --time=00:15:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/gen_pop_%j.out
#SBATCH --error=logs/gen_pop_%j.err

CLUSTER_C="${CLUSTER_C:-0.3}"
export POP_ID="${POP_ID:-$SLURM_ARRAY_TASK_ID}"
mkdir -p logs "data/pop_c${CLUSTER_C}" "data/samples_c${CLUSTER_C}"

module load r/4.5.0

Rscript src/R/generate_population.R
