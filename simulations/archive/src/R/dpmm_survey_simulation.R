# =============================================================================
# Survey-weighted DPMM phenotyping: simulation study
#
# Extends Williams & Savitsky (2021, JRSS-B) PPS1/PPS3 informative-sampling
# designs from a single continuous outcome to a K0-component Gaussian
# mixture outcome. Validates the Godambe sandwich correction (their
# Algorithm 1) applied to the relabeled, permutation-invariant phenotype
# prevalence functional of a survey-weighted truncated Dirichlet Process
# Gaussian Mixture Model.
#
# Scope note: the sandwich correction below is implemented on the
# *relabeled hard-assignment prevalence vector*, treated as a weighted
# multinomial estimating-equation problem (score = w_i*(e_{z_i} - pi),
# Fisher info = diag(pi) - pi pi'). This is a deliberate, documented
# simplification of Algorithm 1 -- applying it to the full raw mixture
# parameter vector (means, covariances, weights) would require a much
# larger Hessian and is not needed here, since the target estimand is
# the prevalence functional, not the raw (unidentified, label-switching-
# prone) component parameters themselves.
# =============================================================================

suppressPackageStartupMessages({
  library(MASS)             # mvrnorm
  library(mvtnorm)          # dmvnorm
  library(label.switching)  # Stephens (2000) relabeling algorithm
  library(coda)             # Rhat / effective sample size diagnostics
  library(Rcpp)
  library(RcppArmadillo)
})

# Load compiled C++ routines (with persistent cache and retry for NFS concurrency)
# Directory layout: this file lives in <PROJECT_ROOT>/src/R/, the Rcpp
# helper in <PROJECT_ROOT>/src/cpp/, and data/results/cache stay at
# <PROJECT_ROOT> (the simulations/ working directory) regardless of CWD.
SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
PROJECT_ROOT <- dirname(dirname(SCRIPT_DIR))
cpp_file <- file.path(dirname(SCRIPT_DIR), "cpp", "dpmm_helpers.cpp")
if (!file.exists(cpp_file)) cpp_file <- "dpmm_helpers.cpp"
cache_dir <- file.path(PROJECT_ROOT, ".rcpp_cache")
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

set.seed(20260918)

# -----------------------------------------------------------------------
# 0. GLOBAL SETTINGS
# -----------------------------------------------------------------------

N_POP          <- 2000000 # population size (before sampling) -- kept >> TARGET_N
N_PSU          <- 200     # number of PSUs in the population (PPS3's z2 layer)
STRATA_H       <- 10      # fixed number of strata
P_DIM          <- 3       # dimension of the CKM-marker-like outcome
K0             <- 4       # true number of phenotypes
TARGET_N       <- 2000    # target sample size after PPS selection
L_TRUNC        <- 10      # DPMM truncation level
N_MCMC_ITER    <- 1500
N_BURNIN       <- 500
SIGMA_FLOOR    <- 0.05    # minimum eigenvalue floor for Sigma_k
R_REPLICATES   <- 1000    # Monte Carlo replicates for production evaluation
B_BOOTSTRAP    <- 200     # PSU-bootstrap replicates for J_hat

# True phenotype parameters (fixed "superpopulation" truth, theta_0)
MU0_LIST <- list(
  c(-2, -2,  0),
  c( 2, -2,  1),
  c(-2,  2, -1),
  c( 2,  2,  0)
)
SIGMA0_LIST <- list(
  diag(c(0.6, 0.6, 0.5)),
  diag(c(0.5, 0.7, 0.5)),
  diag(c(0.6, 0.5, 0.6)),
  diag(c(0.5, 0.5, 0.7))
)

# -----------------------------------------------------------------------
# 1. GENERATIVE MODEL: PPS1 (informative selection) + PPS3 (PSU effect)
# -----------------------------------------------------------------------

generate_population <- function(N = N_POP, n_psu = N_PSU,
                                 H = STRATA_H, beta_x2 = 1.2, z2_mean = 0.5,
                                 include_z2 = TRUE) {
  psu_id <- sample(seq_len(n_psu), N, replace = TRUE)
  psu_stratum <- rep(seq_len(H), length.out = n_psu)[sample(n_psu)]
  stratum <- psu_stratum[psu_id]

  z2_psu <- if (include_z2) rexp(n_psu, rate = 1 / z2_mean) else rep(0, n_psu)
  z2     <- z2_psu[psu_id] - z2_mean * include_z2

  x1 <- rnorm(N)
  x2 <- rgamma(N, shape = 2, rate = 1)

  beta_k  <- c(-0.5, 0.5, -0.25, 0.25)
  eta     <- outer(x2 - 2, beta_k)
  probs   <- exp(eta) / rowSums(exp(eta))
  z_true  <- apply(probs, 1, function(p) sample.int(K0, 1, prob = p))

  Y <- matrix(NA_real_, N, P_DIM)
  for (k in seq_len(K0)) {
    idx <- which(z_true == k)
    if (length(idx) > 0) {
      Y[idx, ] <- MASS::mvrnorm(length(idx), mu = MU0_LIST[[k]], Sigma = SIGMA0_LIST[[k]])
      Y[idx, ] <- Y[idx, ] + z2[idx]
    }
  }

  pi_raw <- x2 / sum(x2)
  list(Y = Y, x1 = x1, x2 = x2, psu_id = psu_id, stratum = stratum,
       z_true = z_true, pi_raw = pi_raw)
}

pps_sample <- function(pop, n_target, H = STRATA_H) {
  idx_all <- integer(0); w_all <- numeric(0); pi_all <- numeric(0)
  N_total <- length(pop$stratum)

  for (h in seq_len(H)) {
    idx_h <- which(pop$stratum == h)
    Nh <- length(idx_h)
    if (Nh == 0) next
    n_h_target <- max(1, round(n_target * Nh / N_total))

    pi_raw_h <- pop$x2[idx_h] / sum(pop$x2[idx_h])
    pi_i_h <- pmin(pi_raw_h * n_h_target, 0.999)
    included_h <- rbinom(length(pi_i_h), 1, pi_i_h) == 1

    idx_all <- c(idx_all, idx_h[included_h])
    w_all   <- c(w_all, 1 / pi_i_h[included_h])
    pi_all  <- c(pi_all, pi_i_h[included_h])
  }

  w_all <- w_all / sum(w_all) * length(idx_all)
  list(idx = idx_all, w = w_all, pi_i = pi_all,
       Y = pop$Y[idx_all, , drop = FALSE], psu_id = pop$psu_id[idx_all],
       stratum = pop$stratum[idx_all],
       z_true = pop$z_true[idx_all], x2 = pop$x2[idx_all])
}

get_true_prevalence <- function(n_mc = 200000) {
  x2 <- rgamma(n_mc, shape = 2, rate = 1)
  beta_k <- c(-0.5, 0.5, -0.25, 0.25)
  eta <- outer(x2 - 2, beta_k)
  probs <- exp(eta) / rowSums(exp(eta))
  colMeans(probs)
}

# -----------------------------------------------------------------------
# 2. SURVEY-WEIGHTED TRUNCATED STICK-BREAKING DPMM
# -----------------------------------------------------------------------

floor_eigen <- function(Sigma, floor = SIGMA_FLOOR) {
  eig <- eigen(Sigma, symmetric = TRUE)
  vals <- pmax(eig$values, floor)
  eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
}

fit_wdpmm <- function(Y, w, x2 = NULL, include_x2 = FALSE,
                       L = L_TRUNC, n_iter = N_MCMC_ITER, n_burnin = N_BURNIN,
                       a_alpha = 1, b_alpha = 4,
                       kappa0 = 0.5, nu0 = P_DIM + 2) {

  n <- nrow(Y); p <- ncol(Y)
  mu0 <- colMeans(Y)
  Psi0 <- diag(apply(Y, 2, var)) * (nu0 - p - 1)

  z <- sample.int(L, n, replace = TRUE)
  mu_k <- lapply(seq_len(L), function(k) mu0 + rnorm(p, 0, 1))
  Sigma_k <- lapply(seq_len(L), function(k) diag(p))
  alpha <- 1
  v <- rep(1 / L, L)

  keep <- (n_burnin + 1):n_iter
  z_draws <- matrix(NA_integer_, length(keep), n)
  min_eig_trace <- numeric(n_iter)
  draw_i <- 0

  for (iter in seq_len(n_iter)) {

    # --- stick-breaking weights pi_k ---
    if (!include_x2) {
      pi_k <- numeric(L)
      remaining <- 1
      for (k in seq_len(L - 1)) {
        n_k  <- sum(w[z == k])
        n_gt <- sum(w[z > k])
        v[k] <- rbeta(1, 1 + n_k, alpha + n_gt)
        pi_k[k] <- v[k] * remaining
        remaining <- remaining * (1 - v[k])
      }
      pi_k[L] <- remaining
      pi_mat <- matrix(pi_k, n, L, byrow = TRUE)
    } else {
      # Fast C++ 2-parameter logistic regression per stick
      pi_mat <- compute_stick_breaking_cpp(z, x2, w, as.integer(L))
    }

    # --- component indicators z_i: p(z_i=k) propto [pi_k * f(y_i|k)]^{w_i} ---
    logp <- matrix(-Inf, n, L)
    for (k in seq_len(L)) {
      dens <- tryCatch(
        mvtnorm::dmvnorm(Y, mean = mu_k[[k]], sigma = Sigma_k[[k]], log = TRUE),
        error = function(e) rep(-1e10, n)
      )
      logp[, k] <- w * (log(pmax(pi_mat[, k], 1e-12)) + dens)
    }
    logp <- logp - apply(logp, 1, max)
    probs <- exp(logp); probs <- probs / rowSums(probs)
    
    # Fast C++ multinomial sampling
    z <- sample_z_cpp(probs)

    # --- component parameters: weighted NIW conjugate update ---
    min_eig_iter <- Inf
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
      
      # Vectorized BLAS crossprod
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

      min_eig_iter <- min(min_eig_iter, min(eigen(Sigma_draw, only.values = TRUE)$values))
    }
    min_eig_trace[iter] <- min_eig_iter

    # --- concentration parameter alpha | v ---
    if (!include_x2) {
      log1mv <- log(pmax(1 - v[1:(L - 1)], 1e-12))
      alpha <- rgamma(1, a_alpha + (L - 1), b_alpha - sum(log1mv))
    }

    if (iter %in% keep) {
      draw_i <- draw_i + 1
      z_draws[draw_i, ] <- z
    }
  }

  list(z_draws = z_draws, min_eig_trace = min_eig_trace,
       mu_k = mu_k, Sigma_k = Sigma_k)
}

# -----------------------------------------------------------------------
# 3. LABEL-SWITCHING-SAFE PREVALENCE FUNCTIONAL (Stephens relabeling)
# -----------------------------------------------------------------------

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

  # --- CENTROID MATCHING TO TRUE PHENOTYPES ---
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

  # One-to-one match of the K0 best-matching fitted components to the K0
  # true phenotypes, fixing a canonical "anchor" component per phenotype.
  assigned_fitted <- integer(K0)
  avail_fitted <- seq_len(K_active)
  for (j in seq_len(K0)) {
    best_k <- avail_fitted[which.min(cost[avail_fitted, j])]
    assigned_fitted[j] <- best_k
    avail_fitted <- setdiff(avail_fitted, best_k)
  }

  # Merge every remaining fitted component (spurious fragments left active
  # by the L_TRUNC > K0 truncation) into whichever true phenotype it is
  # nearest to, rather than discarding its posterior mass. Without this,
  # fragmented phenotypes lose their split-off mass entirely and pi_hat
  # does not sum to 1. NOTE: this merge-to-known-truth step is only valid
  # because K0 and MU0_LIST are known ground truth in this validation
  # study; on real data with unknown K0, fragments would instead need to
  # be merged into each other (e.g. via posterior-similarity/co-clustering
  # or overlap-based component merging), not into a "true" centroid.
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

  # K_active here is the raw number of truncated stick-breaking components
  # kept active (a fragmentation diagnostic); K0 is used downstream for the
  # sandwich estimator's category count now that fragments are merged.
  list(prev_draws = prev_draws_aligned, z_point_estimate = as.integer(z_point_aligned),
       K_active = K0, K_raw_active = K_active)
}

# -----------------------------------------------------------------------
# 4. GODAMBE SANDWICH CORRECTION ON THE PREVALENCE FUNCTIONAL
# -----------------------------------------------------------------------

sandwich_correct_prevalence <- function(prev_draws, w, psu_id, stratum, z_point,
                                         K_active, B = B_BOOTSTRAP) {

  pi_bar <- colMeans(prev_draws)[-K_active]
  Kd <- length(pi_bar)

  H_hat <- (diag(pi_bar, nrow = Kd) - outer(pi_bar, pi_bar)) * sum(w)

  one_hot <- matrix(0, length(z_point), Kd)
  for (k in seq_len(Kd)) one_hot[, k] <- as.numeric(z_point == k)

  # Pre-aggregate PSU scores for the bootstrap
  unique_psus <- unique(psu_id)
  n_psus <- length(unique_psus)
  psu_map <- setNames(seq_len(n_psus) - 1L, as.character(unique_psus)) # 0-based for C++

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

  # Fast C++ bootstrap
  J_hat <- fast_sandwich_bootstrap_cpp(psu_scores, stratum_psu_indices, as.integer(B))

  H_inv <- solve(H_hat + diag(1e-8, Kd))
  Sigma_sandwich <- H_inv %*% J_hat %*% H_inv

  # Cholesky rotation applied to the (reduced) posterior draws
  draws_reduced <- prev_draws[, -K_active, drop = FALSE]
  draws_centered <- sweep(draws_reduced, 2, colMeans(draws_reduced))
  Sigma_naive <- cov(draws_reduced)
  R1 <- tryCatch(chol(Sigma_sandwich), error = function(e) chol(Sigma_sandwich + diag(1e-6, Kd)))
  R2 <- tryCatch(chol(Sigma_naive), error = function(e) chol(Sigma_naive + diag(1e-6, Kd)))
  rotated <- draws_centered %*% solve(R2) %*% R1
  adjusted_draws <- sweep(rotated, 2, colMeans(draws_reduced), "+")

  list(pi_hat = colMeans(prev_draws), adjusted_draws = adjusted_draws,
       Sigma_sandwich = Sigma_sandwich, Sigma_naive = Sigma_naive,
       H_hat = H_hat, J_hat = J_hat)
}

# -----------------------------------------------------------------------
# 5. PER-REPLICATE EXECUTION + REPRODUCIBILITY HOOK
# -----------------------------------------------------------------------

run_one_replicate <- function(include_x2_in_fit = TRUE, include_z2_in_truth = TRUE) {
  pop <- generate_population(include_z2 = include_z2_in_truth)
  samp <- pps_sample(pop, TARGET_N)

  fit <- fit_wdpmm(samp$Y, samp$w, x2 = samp$x2, include_x2 = include_x2_in_fit)

  K_active_min <- K0 + 2
  rl <- relabel_and_summarize(fit, samp$w, samp$Y, K_active_min)
  sw <- sandwich_correct_prevalence(rl$prev_draws, samp$w, samp$psu_id,
                                     samp$stratum, rl$z_point_estimate, rl$K_active)

  list(pi_hat = sw$pi_hat, adjusted_draws = sw$adjusted_draws,
       min_eig_trace = fit$min_eig_trace, K_active = rl$K_active,
       K_raw_active = rl$K_raw_active,
       Sigma_sandwich = sw$Sigma_sandwich, Sigma_naive = sw$Sigma_naive)
}

RESULTS_DIR <- "results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

BASE_SEED <- 123L

seed_for <- function(arm, task_id) {
  # 10,000 offset between arms ensures 1000 tasks per arm never share seeds
  offset <- switch(arm, informed = 0L, naive = 10000L, stop("unknown arm: ", arm))
  BASE_SEED + offset + as.integer(task_id)
}

run_and_save_replicate <- function(arm, task_id) {
  seed <- seed_for(arm, task_id)
  set.seed(seed)

  out_path <- file.path(RESULTS_DIR, sprintf("result_%s_%05d.rds", arm, task_id))
  err_path <- file.path(RESULTS_DIR, sprintf("ERROR_%s_%05d.rds", arm, task_id))

  t0 <- proc.time()
  result <- tryCatch(
    run_one_replicate(include_x2_in_fit = (arm == "informed")),
    error = function(e) {
      saveRDS(list(arm = arm, task_id = task_id, seed = seed,
                   error = conditionMessage(e)),
              err_path)
      NULL
    }
  )
  elapsed <- (proc.time() - t0)["elapsed"]

  if (is.null(result)) {
    cat(sprintf("[FAILED] arm=%s task_id=%d seed=%d elapsed=%.2fs -- see %s\n",
                arm, task_id, seed, elapsed, err_path))
    cat(sprintf("Reproduce with: reproduce_replicate(\"%s\", %d)\n", arm, task_id))
    quit(status = 1, save = "no")
  }

  result$arm <- arm; result$task_id <- task_id; result$seed <- seed
  result$elapsed_seconds <- elapsed
  saveRDS(result, out_path)
  cat(sprintf("[OK] arm=%s task_id=%d seed=%d elapsed=%.2fs -> %s\n",
              arm, task_id, seed, elapsed, out_path))
}

reproduce_replicate <- function(arm, task_id) {
  set.seed(seed_for(arm, task_id))
  run_one_replicate(include_x2_in_fit = (arm == "informed"))
}

aggregate_results <- function(arm) {
  files <- list.files(RESULTS_DIR, pattern = sprintf("^result_%s_.*\\.rds$", arm),
                      full.names = TRUE)
  if (length(files) == 0) stop("No completed result files found for arm: ", arm)

  true_prev <- get_true_prevalence()
  true_reduced <- true_prev[-length(true_prev)]

  n_files <- length(files)
  covered <- matrix(NA, n_files, length(true_reduced))
  K_active_vec <- integer(n_files)
  K_raw_active_vec <- rep(NA_integer_, n_files)
  min_eig_vec <- numeric(n_files)
  elapsed_vec <- numeric(n_files)

  for (i in seq_along(files)) {
    r <- readRDS(files[i])
    ci_lo <- apply(r$adjusted_draws, 2, quantile, probs = 0.025)
    ci_hi <- apply(r$adjusted_draws, 2, quantile, probs = 0.975)
    true_k <- true_reduced[seq_len(ncol(r$adjusted_draws))]
    covered[i, seq_along(true_k)] <- (true_k >= ci_lo[seq_along(true_k)]) &
                                     (true_k <= ci_hi[seq_along(true_k)])
    K_active_vec[i] <- r$K_active
    if (!is.null(r$K_raw_active)) K_raw_active_vec[i] <- r$K_raw_active
    min_eig_vec[i] <- min(r$min_eig_trace, na.rm = TRUE)
    if (!is.null(r$elapsed_seconds)) elapsed_vec[i] <- r$elapsed_seconds
  }

  list(n_completed = n_files,
       empirical_coverage = colMeans(covered, na.rm = TRUE),
       mean_K_active = mean(K_active_vec),
       mean_K_raw_active = mean(K_raw_active_vec, na.rm = TRUE),
       min_eig_summary = summary(min_eig_vec),
       mean_runtime_seconds = mean(elapsed_vec[elapsed_vec > 0]))
}

# -----------------------------------------------------------------------
# 6. ENTRY POINT
# -----------------------------------------------------------------------

slurm_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
slurm_arm     <- Sys.getenv("DPMM_ARM", unset = "informed")
run_local_demo <- Sys.getenv("RUN_LOCAL_DEMO", unset = "0") == "1"

if (!is.na(slurm_task_id)) {
  run_and_save_replicate(slurm_arm, as.integer(slurm_task_id))
} else if (run_local_demo) {
  cat("RUN_LOCAL_DEMO=1 -- running local demo\n")
  for (r in seq_len(min(5, R_REPLICATES))) run_and_save_replicate("informed", r)
  for (r in seq_len(min(5, R_REPLICATES))) run_and_save_replicate("naive", r)
  print(aggregate_results("informed"))
  print(aggregate_results("naive"))
} else {
  cat("Sourced without SLURM_ARRAY_TASK_ID -- defined functions only.\n")
}
