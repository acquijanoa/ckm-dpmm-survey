#!/bin/bash
# One array task = one replicate of src/R/dpmm_fit.R (see its header for the
# method and options). Production defaults need no --export; arms differ by
# MERGE / THRESHOLD / DESIGN_DF, and each lands under its own
# RESULT_TAG (e.g. c0.3_datadriven_thresh_df).
#
#   sbatch --array=1-1000 run_fit.sbatch                                   # production
#   sbatch --array=1-1000 --export=ALL,MERGE=truth run_fit.sbatch          # truth-merge ceiling
#   sbatch --array=1-1000 --export=ALL,THRESHOLD=0 run_fit.sbatch          # no threshold
#   sbatch --array=1-1000 --export=ALL,DESIGN_DF=0 run_fit.sbatch          # no design-df
#
# Measured per task: ~2 min, <1 GB.

#SBATCH --job-name=dpmm_fit
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/fit_%A_%a.out
#SBATCH --error=logs/fit_%A_%a.err

export CLUSTER_C="${CLUSTER_C:-0.3}"
export N_PSU="${N_PSU:-4}"
export MERGE="${MERGE:-datadriven}"
export THRESHOLD="${THRESHOLD:-1}"
export DESIGN_DF="${DESIGN_DF:-1}"
# Same BLAS kernels on every node type, so fits are reproducible across nodes.
export OPENBLAS_CORETYPE=Haswell
mkdir -p logs results .rcpp_cache

module load r/4.5.0

Rscript src/R/dpmm_fit.R
