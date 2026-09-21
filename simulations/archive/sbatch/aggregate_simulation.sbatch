#!/bin/bash
#SBATCH --job-name=dpmm_eval
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/evaluation_1000_final.out
#SBATCH --error=logs/evaluation_1000_final.err

module load r/4.5.0

Rscript -e '
source("src/R/dpmm_survey_simulation.R")

cat("\n========================================================================\n")
cat("      DPMM SURVEY SIMULATION STUDY: 1000-REPLICATE FINAL EVALUATION        \n")
cat("========================================================================\n")

true_prev <- get_true_prevalence()
cat("\n[1] TRUE SUPERPOPULATION PREVALENCES (pi_0):\n")
for (k in seq_along(true_prev)) {
  cat(sprintf("   Phenotype %d: %.4f (%.2f%%)\n", k, true_prev[k], true_prev[k] * 100))
}

evaluate_arm <- function(arm_name) {
  files <- list.files("results", pattern = sprintf("^result_%s_.*\\.rds$", arm_name), full.names = TRUE)
  n_files <- length(files)
  cat(sprintf("\n------------------------------------------------------------------------\n"))
  cat(sprintf("   EVALUATION: %s ARM (n = %d replicates)\n", toupper(arm_name), n_files))
  cat(sprintf("------------------------------------------------------------------------\n"))

  if (n_files == 0) {
    cat("No files found.\n")
    return(NULL)
  }

  covered_all <- matrix(NA, n_files, 4)
  ci_width_all <- matrix(NA, n_files, 4)
  pi_mat <- matrix(NA, n_files, 4)
  K_active_vec <- integer(n_files)
  elapsed_vec <- numeric(n_files)

  for (i in seq_along(files)) {
    r <- readRDS(files[i])
    
    # Adjusted draws for first 3 phenotypes
    adj_draws <- r$adjusted_draws[, 1:3, drop = FALSE]
    
    # Derived draws for 4th phenotype: 1 - sum(pi_1, pi_2, pi_3)
    p4_draws <- 1 - rowSums(adj_draws)
    all_draws <- cbind(adj_draws, p4_draws)
    
    ci_lo <- apply(all_draws, 2, quantile, probs = 0.025)
    ci_hi <- apply(all_draws, 2, quantile, probs = 0.975)
    
    for (k in 1:4) {
      covered_all[i, k] <- (true_prev[k] >= ci_lo[k]) && (true_prev[k] <= ci_hi[k])
      ci_width_all[i, k] <- ci_hi[k] - ci_lo[k]
    }

    if (!is.null(r$pi_hat)) {
      pi_mat[i, ] <- r$pi_hat[1:4]
    }
    K_active_vec[i] <- r$K_active
    if (!is.null(r$elapsed_seconds)) elapsed_vec[i] <- r$elapsed_seconds
  }

  mean_pi <- colMeans(pi_mat, na.rm = TRUE)
  cov_rates <- colMeans(covered_all, na.rm = TRUE)
  mean_ci_width <- colMeans(ci_width_all, na.rm = TRUE)

  cat("\n1. Empirical 95% Coverage & Interval Width (Sandwich-Adjusted):\n")
  for (k in 1:4) {
    cat(sprintf("   Phenotype %d: Coverage = %5.1f%% | Mean 95%% CI Width = %.4f (Nominal: 95.0%%)\n",
                k, cov_rates[k] * 100, mean_ci_width[k]))
  }

  cat("\n2. Point Estimates & Bias vs Truth:\n")
  for (k in 1:4) {
    bias <- mean_pi[k] - true_prev[k]
    rel_bias <- (bias / true_prev[k]) * 100
    cat(sprintf("   Phenotype %d: Est = %.4f | Truth = %.4f | Bias = %+.4f (%+.1f%%)\n",
                k, mean_pi[k], true_prev[k], bias, rel_bias))
  }

  cat("\n3. Clustering & Runtime:\n")
  cat(sprintf("   Mean Active Components (K*): %.2f\n", mean(K_active_vec)))
  valid_elapsed <- elapsed_vec[elapsed_vec > 0]
  if (length(valid_elapsed) > 0) {
    cat(sprintf("   Mean Replicate Runtime:       %.1f s (%.2f min)\n",
                mean(valid_elapsed), mean(valid_elapsed) / 60))
  }
}

evaluate_arm("informed")
evaluate_arm("naive")

cat("\n========================================================================\n\n")
'
