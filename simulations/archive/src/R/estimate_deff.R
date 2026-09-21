# =============================================================================
# Estimate Design Effects (DEFF) from 1,000 Complex Survey Samples
# =============================================================================

suppressPackageStartupMessages({
  library(survey)
})

POP_FILE    <- "data/pop/population.rds"
SAMPLES_DIR <- "data/samples"

cat("========================================================================\n")
cat(" ESTIMATING DESIGN EFFECT (DEFF) OBSERVED ACROSS 1,000 SAMPLES\n")
cat("========================================================================\n")

pop <- readRDS(POP_FILE)
true_pi <- pop$true_prevalence
cat(sprintf("Population Size: N = %d elements across %d PSUs in %d strata\n",
            pop$N_total, pop$total_clusters, pop$H_STRATA))
cat("True Population Prevalences (pi_0):\n")
for (k in seq_along(true_pi)) {
  cat(sprintf("  Phenotype %d: %.4f (%.2f%%)\n", k, true_pi[k], true_pi[k] * 100))
}

sample_files <- list.files(SAMPLES_DIR, pattern = "^sample_.*\\.rds$", full.names = TRUE)
n_samples <- length(sample_files)
cat(sprintf("\nAnalyzing %d completed sample files...\n", n_samples))

# Matrices to store estimates
pi_hat_mat <- matrix(NA_real_, n_samples, 4)
y_mean_mat <- matrix(NA_real_, n_samples, 3)
deff_w_vec <- numeric(n_samples)
deff_svy_mat <- matrix(NA_real_, n_samples, 4)
n_vec <- integer(n_samples)

for (i in seq_along(sample_files)) {
  s <- readRDS(sample_files[i])
  n <- s$n
  n_vec[i] <- n
  w <- s$w

  # 1. Weighted prevalence estimates
  for (k in 1:4) {
    pi_hat_mat[i, k] <- sum(w * (s$z_true == k)) / sum(w)
  }

  # Weighted mean of outcomes Y
  for (p in 1:3) {
    y_mean_mat[i, p] <- sum(w * s$Y[, p]) / sum(w)
  }

  # 2. Kish weighting design effect: 1 + CV(w)^2
  cv_sq <- var(w) / (mean(w)^2)
  deff_w_vec[i] <- 1 + cv_sq

  # 3. Model-based survey design effect via survey package (first 100 samples for speed)
  if (i <= 100) {
    df_s <- data.frame(
      stratum = s$stratum,
      psu = s$psu_id,
      w = s$w,
      z1 = as.numeric(s$z_true == 1),
      z2 = as.numeric(s$z_true == 2),
      z3 = as.numeric(s$z_true == 3),
      z4 = as.numeric(s$z_true == 4)
    )
    des <- svydesign(ids = ~psu, strata = ~stratum, weights = ~w, data = df_s, nest = TRUE)
    svy_means <- svymean(~z1 + z2 + z3 + z4, des, deff = TRUE)
    deff_svy_mat[i, ] <- deff(svy_means)
  }
}

mean_n <- mean(n_vec)
cat(sprintf("\nAverage Sample Size per Replicate: n = %.1f (Sampling Fraction f = %.4f)\n",
            mean_n, mean_n / pop$N_total))

# -----------------------------------------------------------------------
# A. EMPIRICAL MONTE CARLO DESIGN EFFECT
# -----------------------------------------------------------------------
# Var_complex(pi_hat) = empirical variance across the 1,000 independent samples
# Var_srs(pi_hat)     = pi_0 * (1 - pi_0) / n
mc_var_complex <- apply(pi_hat_mat, 2, var)
mc_mean_pi     <- colMeans(pi_hat_mat)
srs_var        <- true_pi * (1 - true_pi) / mean_n
deff_empirical <- mc_var_complex / srs_var

cat("\n------------------------------------------------------------------------\n")
cat(" 1. EMPIRICAL MONTE CARLO DESIGN EFFECT (across 1,000 replicates):\n")
cat("------------------------------------------------------------------------\n")
for (k in 1:4) {
  bias <- mc_mean_pi[k] - true_pi[k]
  rel_bias <- (bias / true_pi[k]) * 100
  cat(sprintf("  Phenotype %d:\n", k))
  cat(sprintf("    True pi_0 = %.4f | Mean Est = %.4f | Bias = %+.4f (%+.2f%%)\n",
              true_pi[k], mc_mean_pi[k], bias, rel_bias))
  cat(sprintf("    Var_complex = %.4e | Var_srs = %.4e\n",
              mc_var_complex[k], srs_var[k]))
  cat(sprintf("    --> EMPIRICAL DEFF = %.3f (Design Factor DEFT = %.3f)\n",
              deff_empirical[k], sqrt(deff_empirical[k])))
}

# -----------------------------------------------------------------------
# B. KISH WEIGHTING DEFF & SURVEY DESIGN DEFF
# -----------------------------------------------------------------------
cat("\n------------------------------------------------------------------------\n")
cat(" 2. DECOMPOSITION OF DESIGN EFFECTS:\n")
cat("------------------------------------------------------------------------\n")
mean_deff_w <- mean(deff_w_vec)
cat(sprintf("  Kish Weighting Effect (DEFF_w = 1 + CV(w)^2): %.3f\n", mean_deff_w))
cat("  (Variation in survey weights due to PPS on x2 inflates variance by ~%.1f%%)\n",
    (mean_deff_w - 1) * 100)

cat("\n  Survey Package Taylor Linearized DEFF (average across samples):\n")
for (k in 1:4) {
  mean_deff_lin <- mean(deff_svy_mat[1:100, k], na.rm = TRUE)
  cat(sprintf("    Phenotype %d: Linearized DEFF = %.3f\n", k, mean_deff_lin))
}

cat("\n------------------------------------------------------------------------\n")
cat(" 3. OUTCOME VARIABLE (Y) EMPIRICAL DEFF:\n")
cat("------------------------------------------------------------------------\n")
# For continuous outcome Y, compare Monte Carlo var to pop_var / n
for (p in 1:3) {
  var_y_complex <- var(y_mean_mat[, p])
  var_y_pop     <- var(pop$Y[, p])
  var_y_srs     <- var_y_pop / mean_n
  deff_y        <- var_y_complex / var_y_srs
  cat(sprintf("  Marker Y%d: Var_complex = %.4e | Var_srs = %.4e | DEFF = %.3f\n",
              p, var_y_complex, var_y_srs, deff_y))
}

cat("\n========================================================================\n")
