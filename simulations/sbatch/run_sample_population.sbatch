#!/bin/bash
# CLUSTER_C selects the clustering-strength grid point (default 0.3),
# matching whichever population run_generate_population.sbatch produced.
# N_PSU selects the PSUs-per-stratum grid point (default 4, the original
# design); other values sample more/fewer, smaller/larger clusters at the
# same overall sampling fraction f. N_SAMPLES sets how many samples to draw
# (default 1000).
#   sbatch --export=ALL,CLUSTER_C=0.3 run_sample_population.sbatch
#   sbatch --export=ALL,CLUSTER_C=1   run_sample_population.sbatch
#   sbatch --export=ALL,CLUSTER_C=1,N_PSU=20 run_sample_population.sbatch
#
# Multi-population variance-decomposition study: submit as an array job
# matching the population array (POP_ID defaults to $SLURM_ARRAY_TASK_ID),
# with N_SAMPLES=100 per population, e.g.
#   sbatch --array=1-25 --export=ALL,CLUSTER_C=0.3,N_SAMPLES=100 run_sample_population.sbatch

#SBATCH --job-name=sample_pop
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/sample_pop_%j.out
#SBATCH --error=logs/sample_pop_%j.err

CLUSTER_C="${CLUSTER_C:-0.3}"
N_PSU="${N_PSU:-4}"
export POP_ID="${POP_ID:-$SLURM_ARRAY_TASK_ID}"
N_SAMPLES="${N_SAMPLES:-1000}"
mkdir -p logs

module load r/4.5.0

Rscript src/R/sample_population.R 1 "$N_SAMPLES"
