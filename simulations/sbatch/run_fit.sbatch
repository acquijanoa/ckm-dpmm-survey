#!/bin/bash
# One array task = one replicate of src/R/dpmm_fit.R (see its header for the
# method and options). Production defaults need no --export; arms differ by
# the options listed in dpmm_fit.R's header, and each lands under its own RESULT_TAG.
#
#   sbatch --array=1-1000 run_fit.sbatch                            # c0.3 (production)
#   sbatch --array=1-1000 --export=ALL,MERGE=truth run_fit.sbatch   # c0.3_truth
#   sbatch --array=1-1000 --export=ALL,THRESHOLD=0 run_fit.sbatch   # c0.3_nothresh
#   sbatch --array=1-1000 --export=ALL,WEIGHT_SUM=neff run_fit.sbatch   # c0.3_neff
#   sbatch --array=1-200 --export=ALL,MERGE=nonempty,ALPHA0=0.01 run_fit.sbatch   # c0.3_nonempty_a0.01
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
export WEIGHT_SUM="${WEIGHT_SUM:-n}"
export ALPHA0="${ALPHA0:-1}"
export POP_SCENARIO="${POP_SCENARIO:-base}"
export L="${L:-10}"
export KAPPA0="${KAPPA0:-0.5}"
export CHAIN="${CHAIN:-1}"
# Same BLAS kernels on every node type, so fits are reproducible across nodes.
export OPENBLAS_CORETYPE=Haswell
mkdir -p logs results .rcpp_cache

module load r/4.5.0

Rscript src/R/dpmm_fit.R
