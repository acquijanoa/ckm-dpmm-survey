#!/bin/bash
# Submit one array job per arm (1000 Monte Carlo replicates each):
#   sbatch --array=1-1000 --export=ALL,DPMM_ARM=informed run_dpmm_array.sbatch
#   sbatch --array=1-1000 --export=ALL,DPMM_ARM=naive    run_dpmm_array.sbatch
#
# Each array task runs exactly one Monte Carlo replicate. The R script
# reads $SLURM_ARRAY_TASK_ID and $DPMM_ARM, derives a deterministic seed
# from (DPMM_ARM, task id), runs one replicate, and saves the result (or,
# on failure, the seed needed to reproduce it) under results/.
#
# After both arrays complete, aggregate via a quick batch/srun call:
#   srun --time=00:10:00 --mem=4G Rscript -e '
#     source("src/R/dpmm_survey_simulation.R")
#     print(aggregate_results("informed"))
#     print(aggregate_results("naive"))'

#SBATCH --job-name=dpmm_sim
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/dpmm_%A_%a.out
#SBATCH --error=logs/dpmm_%A_%a.err

mkdir -p logs results .rcpp_cache

module load r/4.5.0

Rscript src/R/dpmm_survey_simulation.R
