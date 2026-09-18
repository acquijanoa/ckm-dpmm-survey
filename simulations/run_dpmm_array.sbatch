#!/bin/bash
# Submit one array job per arm:
#   sbatch --array=1-500 --export=ALL,DPMM_ARM=informed run_dpmm_array.sbatch
#   sbatch --array=1-500 --export=ALL,DPMM_ARM=naive    run_dpmm_array.sbatch
#
# Each array task runs exactly one Monte Carlo replicate. The R script
# reads $SLURM_ARRAY_TASK_ID and $DPMM_ARM, derives a deterministic seed
# from (DPMM_ARM, task id), runs one replicate, and saves the result (or,
# on failure, the seed needed to reproduce it) under results/.
#
# After both arrays complete, aggregate interactively in R:
#   source("dpmm_survey_simulation.R")   # SLURM_ARRAY_TASK_ID unset locally
#                                         # -> only defines functions, no run
#   print(aggregate_results("informed"))
#   print(aggregate_results("naive"))
#
# To debug a specific failed task (e.g. naive, task 340):
#   rep <- reproduce_replicate("naive", 340)

#SBATCH --job-name=dpmm_sim
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/dpmm_%A_%a.out
#SBATCH --error=logs/dpmm_%A_%a.err

mkdir -p logs results

module load R   # adjust to the cluster's module system

Rscript dpmm_survey_simulation.R
