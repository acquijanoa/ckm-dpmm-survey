#!/bin/bash
# Validate the sparse-Dirichlet + unweighted-z fit + H_hat-fixed sandwich
# on many independent samples from a SINGLE FIXED population
# (data/pop/population.rds + data/samples/sample_*.rds), isolating true
# survey-sampling variance from the simulation's population-regeneration
# noise.
#   sbatch --array=1-100 run_dpmm_fixedpop.sbatch

#SBATCH --job-name=dpmm_fixedpop
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/dpmm_%A_%a.out
#SBATCH --error=logs/dpmm_%A_%a.err

mkdir -p logs results_fixedpop .rcpp_cache_sparse

module load r/4.5.0

Rscript src/R/dpmm_fixedpop_validate.R
