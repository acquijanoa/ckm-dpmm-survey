#!/bin/bash
# Raw-mixture + total-variance-recombination Godambe correction fit --
# see src/R/dpmm_fit.R header for the derivation. This is the only active
# fitting method; superseded variants live under archive/.
# CLUSTER_C selects the clustering-strength grid point (matches
# generate_population.R / sample_population.R); defaults to 0.3 if unset.
# N_PSU selects the PSUs-per-stratum grid point (default 4); other values
# must match a sample_population.R run submitted with the same N_PSU.
# APPLY_CORR (0/1, default 0) toggles the design-df sandwich scale (c_corr)
# and t-quantile width corrections in sandwich_correct_prevalence.
# ORACLE (0/1, default 0) replaces the fitted-component plug-in a_ik with
# the TRUE superpopulation Gaussian densities, isolating whether the
# residual coverage gap is attributable to unpropagated (mu_k,Sigma_k)
# posterior uncertainty.
# UNKNOWN_K0 (0/1, default 0) replaces the truth-based raw-component merge
# with a data-driven (posterior-similarity-matrix) merge, for testing
# whether the method still recovers K0 and calibration without ground
# truth -- as would be required on real data.
#   sbatch --array=1-N --export=ALL,CLUSTER_C=0.3 \
#     --output=logs/fit_c0.3_%A_%a.out \
#     --error=logs/fit_c0.3_%A_%a.err \
#     run_fit.sbatch
#   sbatch --array=1-N --export=ALL,CLUSTER_C=1 \
#     --output=logs/fit_c1_%A_%a.out \
#     --error=logs/fit_c1_%A_%a.err \
#     run_fit.sbatch
#   sbatch --array=1-N --export=ALL,CLUSTER_C=1,N_PSU=20 \
#     --output=logs/fit_c1_npsu20_%A_%a.out \
#     --error=logs/fit_c1_npsu20_%A_%a.err \
#     run_fit.sbatch

#SBATCH --job-name=dpmm_fit
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/fit_%A_%a.out
#SBATCH --error=logs/fit_%A_%a.err

CLUSTER_C="${CLUSTER_C:-0.3}"
N_PSU="${N_PSU:-4}"
export APPLY_CORR="${APPLY_CORR:-0}"
export ORACLE="${ORACLE:-0}"
export UNKNOWN_K0="${UNKNOWN_K0:-0}"
export MERGE_THRESHOLD="${MERGE_THRESHOLD:-0}"
export MERGE_MIN_WEIGHT="${MERGE_MIN_WEIGHT:-0.02}"
export MERGE_MIN_SEP="${MERGE_MIN_SEP:-1.5}"
mkdir -p logs results .rcpp_cache

module load r/4.5.0

Rscript src/R/dpmm_fit.R
