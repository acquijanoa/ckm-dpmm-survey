#!/bin/bash
# Draw N_SAMPLES two-stage samples from one population.
#   sbatch --export=ALL,CLUSTER_C=0.3,N_PSU=4 run_sample_population.sbatch

#SBATCH --job-name=sample_pop
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/sample_pop_%j.out
#SBATCH --error=logs/sample_pop_%j.err

export CLUSTER_C="${CLUSTER_C:-0.3}"
export POP_SCENARIO="${POP_SCENARIO:-base}"
export N_PSU="${N_PSU:-4}"
N_SAMPLES="${N_SAMPLES:-1000}"
mkdir -p logs

module load r/4.5.0

Rscript src/R/sample_population.R 1 "$N_SAMPLES"
