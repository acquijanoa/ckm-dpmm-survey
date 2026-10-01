#!/bin/bash
# One array task = one replicate of src/R/dpmm_fit.R (see its header for the
# method and options). Production defaults need no --export; arms differ by
# MERGE / THRESHOLD / N_ITER / N_BURNIN, and each lands under its own RESULT_TAG.
#
#   sbatch --array=1-1000 run_fit.sbatch                            # c0.3 (production)
#   sbatch --array=1-1000 --export=ALL,MERGE=truth run_fit.sbatch   # c0.3_truth
#   sbatch --array=1-1000 --export=ALL,THRESHOLD=0 run_fit.sbatch   # c0.3_nothresh
#   sbatch --array=1-200 --mem=2G --time=00:15:00 --export=ALL,N_ITER=4000,N_BURNIN=2000 run_fit.sbatch   # c0.3_it4000b2000
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
export N_ITER="${N_ITER:-1500}"
export N_BURNIN="${N_BURNIN:-500}"
# Same BLAS kernels on every node type, so fits are reproducible across nodes.
export OPENBLAS_CORETYPE=Haswell
mkdir -p logs results .rcpp_cache

module load r/4.5.0

Rscript src/R/dpmm_fit.R
