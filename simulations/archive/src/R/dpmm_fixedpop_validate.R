# =============================================================================
# Validation on a FIXED population with many independent samples
# (data/pop/population.rds + data/samples/sample_*.rds)
#
# Unlike dpmm_survey_simulation_sparse.R, which regenerates a brand-new
# finite population on every replicate (conflating population-regeneration
# noise with actual survey-sampling variance in the between-replicate
# variance used to diagnose the Godambe sandwich correction), this script
# holds the population FIXED and only varies the sample -- the exact
# quantity the sandwich correction is designed to estimate. It reuses the
# sparse-Dirichlet + unweighted-z fit (dpmm_helpers_sparse.cpp) and the
# H_hat-fixed sandwich correction from dpmm_survey_simulation_sparse.R,
# applied to samples with real ~1,250-cluster PSU structure (200 PSUs
# realized per sample) rather than the simulation's simpler 200-PSU design.
# =============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
  library(label.switching)
  library(coda)
  library(Rcpp)
  library(RcppArmadillo)
})

# Directory layout: this file lives in <PROJECT_ROOT>/src/R/, the Rcpp
# helper in <PROJECT_ROOT>/src/cpp/, and data/results/cache stay at
# <PROJECT_ROOT> (the simulations/ working directory) regardless of CWD.
SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
PROJECT_ROOT <- dirname(dirname(SCRIPT_DIR))
cpp_file <- file.path(dirname(SCRIPT_DIR), "cpp", "dpmm_helpers_sparse.cpp")
if (!file.exists(cpp_file)) cpp_file <- "dpmm_helpers_sparse.cpp"
cache_dir <- file.path(PROJECT_ROOT, ".rcpp_cache_sparse")
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

loaded <- FALSE
for (attempt in 1:6) {
  tryCatch({
    Rcpp::sourceCpp(cpp_file, cacheDir = cache_dir)
    loaded <- TRUE
  }, error = function(e) {
    Sys.sleep(runif(1, 0.2, 2.0))
  })
  if (loaded) break
}
if (!loaded) Rcpp::sourceCpp(cpp_file)

POP_PATH <- file.path(PROJECT_ROOT, "data", "pop", "population.rds")
SAMPLES_DIR <- file.path(PROJECT_ROOT, "data", "samples")

pop <- readRDS(POP_PATH)
K0 <- length(pop$MU0_LIST)
P_DIM <- ncol(pop$Y)
MU0_LIST <- pop$MU0_LIST
TRUE_PREVALENCE <- pop$true_prevalence

L_TRUNC       <- 10
N_MCMC_ITER   <- 1500
N_BURNIN      <- 500
SIGMA_FLOOR   <- 0.05
ALPHA0_SPARSE <- 1
B_BOOTSTRAP   <- 200

floor_eigen <- function(Sigma, floor = SIGMA_FLOOR) {
  eig <- eigen(Sigma, symmetric = TRUE)
  vals <- pmax(eig$values, floor)
  eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
}

# Unweighted-z classification + sparse Dirichlet(alpha0/L) mixture weights,
# identical mechanism to dpmm_survey_simulation_sparse.R's naive arm.
fit_wdpmm <- function(Y, w, L = L_TRUNC, n_iter = N_MCMC_ITER, n_burnin = N_BURNIN,
                       alpha0 = ALPHA0_SPARSE, kappa0 = 0.5, nu0 = P_DIM + 2) {
  n <- nrow(Y); p <- ncol(Y)
  mu0 <- colMeans(Y)
  Psi0 <- diag(apply(Y, 2, var)) * (nu0 - p - 1)

  z <- sample.int(L, n, replace = TRUE)
  mu_k <- lapply(seq_len(L), function(k) mu0 + rnorm(p, 0, 1))
  Sigma_k <- lapply(seq_len(L), function(k) diag(p))

  keep <- (n_burnin + 1):n_iter
  z_draws <- matrix(NA_integer_, length(keep), n)
  min_eig_trace <- numeric(n_iter)
  draw_i <- 0

  for (iter in seq_len(n_iter)) {
    pi0 <- sparse_dirichlet_weights_cpp(z, w, as.integer(L), alpha0)
    pi_mat <- matrix(pi0, n, L, byrow = TRUE)

    logp <- matrix(-Inf, n, L)
    for (k in seq_len(L)) {
      dens <- tryCatch(
        mvtnorm::dmvnorm(Y, mean = mu_k[[k]], sigma = Sigma_k[[k]], log = TRUE),
        error = function(e) rep(-1e10, n)
      )
      logp[, k] <- log(pmax(pi_mat[, k], 1e-12)) + dens
    }
    logp <- logp - apply(logp, 1, max)
    probs <- exp(logp); probs <- probs / rowSums(probs)
    z <- sample_z_cpp(probs)

    for (k in seq_len(L)) {
      idx <- which(z == k)
      if (length(idx) == 0) {
        mu_k[[k]] <- mu0 + MASS::mvrnorm(1, rep(0, p), Psi0 / (nu0 - p - 1))
        Sigma_k[[k]] <- floor_eigen(solve(rWishart(1, nu0, solve(Psi0))[, , 1]))
        next
      }
      wk <- w[idx]; Yk <- Y[idx, , drop = FALSE]
      nk_eff <- sum(wk)
      ybar_k <- colSums(Yk * wk) / nk_eff
      d <- sweep(Yk, 2, ybar_k)
      Sk <- crossprod(d * sqrt(wk))
      kappa_n <- kappa0 + nk_eff
      nu_n    <- nu0 + nk_eff
      mu_n    <- (kappa0 * mu0 + nk_eff * ybar_k) / kappa_n
      d0 <- ybar_k - mu0
      Psi_n <- Psi0 + Sk + (kappa0 * nk_eff / kappa_n) * outer(d0, d0)
      Sigma_draw <- tryCatch(
        solve(rWishart(1, nu_n, solve(Psi_n))[, , 1]),
        error = function(e) diag(p)
      )
      Sigma_draw <- floor_eigen(Sigma_draw)
      Sigma_k[[k]] <- Sigma_draw
      mu_k[[k]] <- MASS::mvrnorm(1, mu_n, Sigma_draw / kappa_n)
    }

    if (iter %in% keep) {
      draw_i <- draw_i + 1
      z_draws[draw_i, ] <- z
    }
  }
  list(z_draws = z_draws, mu_k = mu_k, Sigma_k = Sigma_k)
}

relabel_and_summarize <- function(fit, w, Y, K_active_min, coverage = 0.99) {
  z_draws <- fit$z_draws
  occ_sorted <- sort(table(z_draws), decreasing = TRUE)
  cum_frac <- cumsum(occ_sorted) / sum(occ_sorted)
  K_active <- max(K_active_min, min(which(cum_frac >= coverage)))
  K_active <- min(K_active, length(occ_sorted))
  top_k <- as.integer(names(occ_sorted)[seq_len(K_active)])
  remap <- setNames(seq_along(top_k), top_k)
  z_collapsed <- matrix(sapply(z_draws, function(x) {
    if (as.character(x) %in% names(remap)) remap[as.character(x)] else NA
  }), nrow(z_draws), ncol(z_draws))

  valid_rows <- which(rowSums(is.na(z_collapsed)) == 0)
  if (length(valid_rows) < 10) {
    K_active <- min(K_active + 2, length(occ_sorted))
    top_k <- as.integer(names(occ_sorted)[seq_len(K_active)])
    remap <- setNames(seq_along(top_k), top_k)
    z_collapsed <- matrix(sapply(z_draws, function(x) {
      if (as.character(x) %in% names(remap)) remap[as.character(x)] else NA
    }), nrow(z_draws), ncol(z_draws))
    valid_rows <- which(rowSums(is.na(z_collapsed)) == 0)
  }
  z_valid <- z_collapsed[valid_rows, , drop = FALSE]

  m <- nrow(z_valid); n <- ncol(z_valid)
  p_alloc <- array(0, dim = c(m, n, K_active))
  for (i in seq_len(m)) {
    for (k in seq_len(K_active)) p_alloc[i, , k] <- as.numeric(z_valid[i, ] == k)
  }
  class_mat <- matrix(as.integer(z_valid), m, n)

  stephens_res <- label.switching::stephens(p_alloc)
  perm <- stephens_res$permutations

  prev_draws <- matrix(NA_real_, m, K_active)
  z_relabeled_last <- NULL
  for (i in seq_len(m)) {
    cur_perm <- perm[i, ]
    z_relab <- cur_perm[class_mat[i, ]]
    tab <- numeric(K_active)
    for (k in seq_len(K_active)) tab[k] <- sum(w[z_relab == k])
    prev_draws[i, ] <- tab / sum(w)
    if (i == m) z_relabeled_last <- z_relab
  }

  centroids <- matrix(NA_real_, K_active, ncol(Y))
  for (k in seq_len(K_active)) {
    idx <- which(z_relabeled_last == k)
    if (length(idx) > 0) {
      centroids[k, ] <- colSums(Y[idx, , drop = FALSE] * w[idx]) / sum(w[idx])
    } else {
      centroids[k, ] <- rep(1e6, ncol(Y))
    }
  }

  cost <- matrix(0, K_active, K0)
  for (k in seq_len(K_active)) {
    for (j in seq_len(K0)) {
      cost[k, j] <- sqrt(sum((centroids[k, ] - MU0_LIST[[j]])^2))
    }
  }

  assigned_fitted <- integer(K0)
  avail_fitted <- seq_len(K_active)
  for (j in seq_len(K0)) {
    best_k <- avail_fitted[which.min(cost[avail_fitted, j])]
    assigned_fitted[j] <- best_k
    avail_fitted <- setdiff(avail_fitted, best_k)
  }

  component_to_phenotype <- integer(K_active)
  component_to_phenotype[assigned_fitted] <- seq_len(K0)
  if (length(avail_fitted) > 0) {
    component_to_phenotype[avail_fitted] <-
      apply(cost[avail_fitted, , drop = FALSE], 1, which.min)
  }

  prev_draws_aligned <- matrix(0, m, K0)
  for (k in seq_len(K_active)) {
    j <- component_to_phenotype[k]
    prev_draws_aligned[, j] <- prev_draws_aligned[, j] + prev_draws[, k]
  }
  z_point_aligned <- component_to_phenotype[z_relabeled_last]

  list(prev_draws = prev_draws_aligned, z_point_estimate = as.integer(z_point_aligned),
       K_active = K0, K_raw_active = K_active)
}

# H_hat-FIXED Godambe sandwich correction (sum(w)*I bread matrix).
sandwich_correct_prevalence <- function(prev_draws, w, psu_id, stratum, z_point,
                                         K_active, B = B_BOOTSTRAP) {
  pi_bar <- colMeans(prev_draws)[-K_active]
  Kd <- length(pi_bar)

  H_hat <- diag(sum(w), Kd)

  one_hot <- matrix(0, length(z_point), Kd)
  for (k in seq_len(Kd)) one_hot[, k] <- as.numeric(z_point == k)

  unique_psus <- unique(psu_id)
  n_psus <- length(unique_psus)
  psu_map <- setNames(seq_len(n_psus) - 1L, as.character(unique_psus))

  psu_scores <- matrix(0, n_psus, Kd)
  unit_scores <- w * (one_hot - matrix(pi_bar, length(z_point), Kd, byrow = TRUE))
  for (i in seq_along(psu_id)) {
    p_idx <- psu_map[[as.character(psu_id[i])]] + 1L
    psu_scores[p_idx, ] <- psu_scores[p_idx, ] + unit_scores[i, ]
  }

  strata_ids <- unique(stratum)
  stratum_psu_indices <- lapply(strata_ids, function(h) {
    h_psus <- unique(psu_id[stratum == h])
    as.integer(psu_map[as.character(h_psus)])
  })

  J_hat <- fast_sandwich_bootstrap_cpp(psu_scores, stratum_psu_indices, as.integer(B))

  H_inv <- solve(H_hat + diag(1e-8, Kd))
  Sigma_sandwich <- H_inv %*% J_hat %*% H_inv

  draws_reduced <- prev_draws[, -K_active, drop = FALSE]
  draws_centered <- sweep(draws_reduced, 2, colMeans(draws_reduced))
  Sigma_naive <- cov(draws_reduced)
  R1 <- tryCatch(chol(Sigma_sandwich), error = function(e) chol(Sigma_sandwich + diag(1e-6, Kd)))
  R2 <- tryCatch(chol(Sigma_naive), error = function(e) chol(Sigma_naive + diag(1e-6, Kd)))
  rotated <- draws_centered %*% solve(R2) %*% R1
  adjusted_draws <- sweep(rotated, 2, colMeans(draws_reduced), "+")

  list(pi_hat = colMeans(prev_draws), adjusted_draws = adjusted_draws,
       Sigma_sandwich = Sigma_sandwich, Sigma_naive = Sigma_naive)
}

run_one_sample <- function(sample_id) {
  samp_path <- file.path(SAMPLES_DIR, sprintf("sample_%04d.rds", sample_id))
  samp <- readRDS(samp_path)

  fit <- fit_wdpmm(samp$Y, samp$w)
  K_active_min <- K0 + 2
  rl <- relabel_and_summarize(fit, samp$w, samp$Y, K_active_min)
  sw <- sandwich_correct_prevalence(rl$prev_draws, samp$w, samp$psu_id,
                                     samp$stratum, rl$z_point_estimate, rl$K_active)

  list(pi_hat = sw$pi_hat, adjusted_draws = sw$adjusted_draws,
       K_active = rl$K_active, K_raw_active = rl$K_raw_active)
}

RESULTS_DIR <- file.path(PROJECT_ROOT, "results_fixedpop")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

run_and_save_sample <- function(sample_id) {
  out_path <- file.path(RESULTS_DIR, sprintf("result_%04d.rds", sample_id))
  err_path <- file.path(RESULTS_DIR, sprintf("ERROR_%04d.rds", sample_id))

  t0 <- proc.time()
  result <- tryCatch(
    run_one_sample(sample_id),
    error = function(e) {
      saveRDS(list(sample_id = sample_id, error = conditionMessage(e)), err_path)
      NULL
    }
  )
  elapsed <- (proc.time() - t0)["elapsed"]

  if (is.null(result)) {
    cat(sprintf("[FAILED] sample_id=%d elapsed=%.2fs -- see %s\n", sample_id, elapsed, err_path))
    quit(status = 1, save = "no")
  }
  result$sample_id <- sample_id
  result$elapsed_seconds <- elapsed
  saveRDS(result, out_path)
  cat(sprintf("[OK] sample_id=%d elapsed=%.2fs -> %s\n", sample_id, elapsed, out_path))
}

aggregate_fixedpop_results <- function() {
  files <- list.files(RESULTS_DIR, pattern = "^result_.*\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No completed result files found.")

  n_files <- length(files)
  pi_mat <- matrix(NA, n_files, K0)
  covered <- matrix(NA, n_files, K0)
  ci_widths <- matrix(NA, n_files, K0)
  K_raw_vec <- integer(n_files)

  for (i in seq_along(files)) {
    r <- readRDS(files[i])
    pi_mat[i, ] <- r$pi_hat[1:K0]
    adj <- r$adjusted_draws[, seq_len(K0 - 1), drop = FALSE]
    last_col <- 1 - rowSums(adj)
    all_draws <- cbind(adj, last_col)
    ci_lo <- apply(all_draws, 2, quantile, probs = 0.025)
    ci_hi <- apply(all_draws, 2, quantile, probs = 0.975)
    ci_widths[i, ] <- ci_hi - ci_lo
    covered[i, ] <- (TRUE_PREVALENCE >= ci_lo) & (TRUE_PREVALENCE <= ci_hi)
    K_raw_vec[i] <- r$K_raw_active
  }

  mean_pi <- colMeans(pi_mat)
  emp_se <- apply(pi_mat, 2, sd)
  mean_width <- colMeans(ci_widths)
  implied_se <- mean_width / (2 * 1.96)

  list(n_completed = n_files,
       true_prevalence = TRUE_PREVALENCE,
       mean_pi_hat = mean_pi,
       bias = mean_pi - TRUE_PREVALENCE,
       rel_bias_pct = 100 * (mean_pi - TRUE_PREVALENCE) / TRUE_PREVALENCE,
       empirical_se = emp_se,
       implied_se = implied_se,
       coverage_pct = colMeans(covered) * 100,
       mean_K_raw_active = mean(K_raw_vec))
}

# -----------------------------------------------------------------------
# ENTRY POINT
# -----------------------------------------------------------------------
slurm_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
if (!is.na(slurm_task_id)) {
  run_and_save_sample(as.integer(slurm_task_id))
} else {
  cat("Sourced without SLURM_ARRAY_TASK_ID -- defined functions only.\n")
}
