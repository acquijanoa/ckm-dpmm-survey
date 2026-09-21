#!/bin/bash
# Raw-mixture (marginal-likelihood, unmerged fitted components) Godambe
# correction variant -- see dpmm_fixedpop_validate_rawmix.R header.
#   sbatch --array=1-N run_dpmm_fixedpop_rawmix.sbatch

#SBATCH --job-name=dpmm_fixedpop_rawmix
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/rawmix_%A_%a.out
#SBATCH --error=logs/rawmix_%A_%a.err

mkdir -p logs results_fixedpop_rawmix .rcpp_cache_sparse

module load r/4.5.0

Rscript src/R/dpmm_fixedpop_validate_rawmix.R
