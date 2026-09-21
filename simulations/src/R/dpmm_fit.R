# =============================================================================
# Validation on a FIXED population with many independent samples -- RAW-
# MIXTURE + TOTAL-VARIANCE RECOMBINATION Godambe correction variant.
#
# dpmm_fixedpop_validate_rawmix.R (soft responsibilities from the true,
# unmerged raw fitted components) reduced the SE ratio (empirical/
# sandwich-implied) from 1.47-2.20 to 1.25-1.49, and a centered-coverage
# check (Test A: coverage of [pi_bar +/- 1.96*SE] around each replicate's
# own pi_hat, isolating variance calibration from location bias) still
# showed only 85-90% vs. the 95% target, confirming this closes PART but
# not all of the gap. This variant adds a second correction, law-of-
# total-variance recombination: Sigma_sandwich (design/classification
# variability given fixed mu_k, Sigma_k) and Sigma_naive (the raw MCMC
# posterior covariance, reflecting (mu_k,Sigma_k)-mixing/parameter
# uncertainty within one chain) are treated as additive, roughly
# independent variance components, so the Cholesky-rotation target
# becomes Sigma_sandwich + Sigma_naive rather than Sigma_sandwich alone.
# A back-of-envelope check on earlier (non-fixed-population) diagnostics
# suggested this could bring the SE ratio down to ~1.28 by itself; here
# it is combined with the already-improved raw-mixture Sigma_sandwich.
#
# dpmm_fixedpop_validate.R (the H_hat-fixed baseline) showed a stable,
# real undercoverage even with the H_hat bug fixed and even holding the
# population fixed: empirical SE exceeds the sandwich-implied SE by a
# consistent 1.5-2.2x across phenotypes. dpmm_fixedpop_validate_softresp.R
# tried fixing this by building the Godambe H/J from the MARGINAL mixture
# log-likelihood (soft responsibilities) instead of conditioning on the
# hard classification z_point, following Wu/Stephenson's (SWOLCA)
# adaptation of Williams & Savitsky confirmed from their Stan code
# (WSOLCA_main.stan: target += w_i * log_sum_exp(log_cond_c[i])) -- but
# using a single plug-in Gaussian per MERGED phenotype. That showed no
# consistent improvement (n=200 pilot), plausibly because a single
# Gaussian for an already-fragmented (K_raw_active ~ 6.4) phenotype is
# itself too "confident," understating uncertainty almost as much as the
# hard classification did.
#
# This variant reconstructs the TRUE sub-mixture density per phenotype
# from the ORIGINAL, unmerged raw fitted components (relabel_and_
# summarize's new a_mat output: for each phenotype j, a_ij is the
# properly-normalized mixture of every raw component mapped to j, using
# the last iteration's raw component weights and each raw component's own
# (mu_k, Sigma_k) -- not a single Gaussian). This is a more faithful
# analogue of SWOLCA's marginal likelihood, since their theta/pi
# parameters are the model's own fitted mixture components, not a
# post-hoc merge. The score/H_hat construction from a_mat is otherwise
# identical to the softresp variant:
#   m_i        = sum_j pi_bar_j * a_ij
#   score_i,k  = w_i*(a_ik - a_iK0)/m_i
#   H_kl       = sum_i w_i*(a_ik-a_iK0)*(a_il-a_iK0)/m_i^2
# Everything else -- relabeling/merge, the PSU/stratum bootstrap
# structure, and the Cholesky-rotation step -- is unchanged (that
# rotation already exactly implements Williams & Savitsky's
# Theta_m^a = (Theta_m - Theta_bar) R2^{-1} R1 + Theta_bar).
# =============================================================================

personal_lib <- "/nas/longleaf/home/aquijano/R/x86_64-pc-linux-gnu-library/4.5"
if (dir.exists(personal_lib) && !(personal_lib %in% .libPaths())) {
  .libPaths(c(personal_lib, .libPaths()))
}

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
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1])))
  sys.frame(1)$ofile  # fall back to the source()-invocation case
}
SCRIPT_DIR <- tryCatch(dirname(get_script_path()), error = function(e) getwd())
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

# Clustering-strength grid parameter (matches generate_population.R /
# sample_population.R): selects which population/samples suffix to read.
CLUSTER_C_STR <- Sys.getenv("CLUSTER_C", unset = "0.3")
CLUSTER_C     <- as.numeric(CLUSTER_C_STR)

# N_PSU grid parameter (matches sample_population.R): number of PSUs
# sampled per stratum. Default 4 matches the original (pre-grid) design and
# its folder/filenames; other values get an explicit _npsuN suffix.
N_PSU_STR <- Sys.getenv("N_PSU", unset = "4")
N_CLUSTERS_SAMPLE <- as.integer(N_PSU_STR)

# ARM_TAG: c<C>[_npsu<N>] -- N_PSU only appends a suffix when it differs
# from its original (pre-grid) default, so already-computed baselines keep
# their original folder/file names. Used for the pop/samples folder-per-arm
# layout (bulk du/tar/rm) -- samples don't depend on APPLY_CORR, only on
# CLUSTER_C/N_PSU, so this stays untagged by correction setting.
ARM_TAG_ENV <- Sys.getenv("ARM_TAG", unset = "")
if (nzchar(ARM_TAG_ENV)) {
  ARM_TAG <- ARM_TAG_ENV
} else {
  ARM_TAG <- sprintf("c%s", CLUSTER_C_STR)
  if (N_PSU_STR != "4") ARM_TAG <- paste0(ARM_TAG, "_npsu", N_PSU_STR)
}

# APPLY_CORR toggles the design-df sandwich scale correction (c_corr) and
# t-quantile width correction added in sandwich_correct_prevalence. Default
# "0" (corrections OFF) -- evaluate the new cluster-probability DGP on the
# base sandwich method first; set APPLY_CORR=1 to layer the correction on.
APPLY_CORR_STR <- Sys.getenv("APPLY_CORR", unset = "0")
APPLY_CORR <- APPLY_CORR_STR == "1"

# ORACLE toggle: replace the raw-mixture reconstruction of a_ik (plug-in,
# from the FITTED components) with the TRUE superpopulation Gaussian
# densities (MU0_LIST/SIGMA0_LIST) in relabel_and_summarize's a_mat. Tests
# whether the residual coverage gap is attributable to unpropagated
# posterior uncertainty in the fitted (mu_k, Sigma_k): if oracle coverage
# reaches nominal, that confirms it; if not, the gap has another source.
ORACLE_STR <- Sys.getenv("ORACLE", unset = "0")
ORACLE <- ORACLE_STR == "1"

# RESULT_TAG: like ARM_TAG but also carries a _corr1/_oracle1 suffix when
# APPLY_CORR/ORACLE deviate from their defaults (0), so different fits of
# the SAME samples land in distinctly-named results/posteriors rather than
# one overwriting/skipping (via resumability) the other. Respects an
# explicit ARM_TAG override exactly (no auto-suffix), matching that
# override's intent.
RESULT_TAG <- ARM_TAG
if (!nzchar(ARM_TAG_ENV)) {
  if (APPLY_CORR) RESULT_TAG <- paste0(RESULT_TAG, "_corr1")
  if (ORACLE) RESULT_TAG <- paste0(RESULT_TAG, "_oracle1")
}

# POP_PATH/SAMPLES_DIR env overrides let a one-off arm point at any
# population/samples location without touching the CLUSTER_C/N_PSU grid
# convention below.
POP_PATH_ENV <- Sys.getenv("POP_PATH", unset = "")
SAMPLES_DIR_ENV <- Sys.getenv("SAMPLES_DIR", unset = "")

if (nzchar(POP_PATH_ENV) && nzchar(SAMPLES_DIR_ENV)) {
  POP_PATH <- POP_PATH_ENV
  SAMPLES_DIR <- SAMPLES_DIR_ENV
} else {
  POP_PATH <- file.path(PROJECT_ROOT, "data", sprintf("pop_c%s", CLUSTER_C_STR), "population.rds")
  SAMPLES_DIR <- file.path(PROJECT_ROOT, "data", sprintf("samples_%s", ARM_TAG))
}

pop <- readRDS(POP_PATH)
K0 <- length(pop$MU0_LIST)
P_DIM <- ncol(pop$Y)
MU0_LIST <- pop$MU0_LIST
SIGMA0_LIST <- pop$SIGMA0_LIST
TRUE_PREVALENCE <- pop$true_prevalence

L_TRUNC       <- 10
N_MCMC_ITER   <- 1500
N_BURNIN      <- 500
SIGMA_FLOOR   <- 0.05
ALPHA0_SPARSE <- 1

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

relabel_and_summarize <- function(fit, w, Y, K_active_min, coverage = 0.99, oracle = ORACLE) {
  z_draws <- fit$z_draws
  occ_sorted <- sort(table(z_draws), decreasing = TRUE)
  cum_frac <- cumsum(occ_sorted) / sum(occ_sorted)
  K_active <- max(K_active_min, min(which(cum_frac >= coverage)))
  K_active <- min(K_active, length(occ_sorted))
  top_k <- as.integer(names(occ_sorted)[seq_len(K_active)])
  remap <- setNames(seq_along(top_k), top_k)

  # ZERO-WASTE RETENTION: rather than dropping an entire MCMC iteration
  # whenever ANY unit's raw label falls outside the top-K_active set
  # (previously discarded, on average, over half of the post-burnin
  # chain -- see ESS diagnostic that motivated this), remap every
  # inactive label to its NEAREST active label by Euclidean distance
  # between the final-iteration fitted centroids fit$mu_k. Every kept
  # iteration is retained, so Sigma_naive is estimated from the full
  # chain instead of a small, non-random surviving subset of it.
  all_labels <- as.integer(names(occ_sorted))
  inactive_labels <- setdiff(all_labels, top_k)
  if (length(inactive_labels) > 0) {
    top_k_mu <- do.call(rbind, fit$mu_k[top_k])
    nearest_top_k <- sapply(inactive_labels, function(ell) {
      d <- sqrt(rowSums(sweep(top_k_mu, 2, fit$mu_k[[ell]])^2))
      top_k[which.min(d)]
    })
    remap <- c(remap, setNames(remap[as.character(nearest_top_k)], inactive_labels))
  }
  z_collapsed <- matrix(remap[as.character(z_draws)], nrow(z_draws), ncol(z_draws))
  z_valid <- z_collapsed

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

  n_units <- ncol(z_valid)
  a_mat <- matrix(0, n_units, K0)

  if (oracle) {
    # ORACLE variant: use the TRUE superpopulation Gaussian density per
    # phenotype directly (MU0_LIST/SIGMA0_LIST), bypassing the fitted
    # (mu_k, Sigma_k) entirely -- isolates whether unpropagated posterior
    # uncertainty in the fitted components explains the residual coverage
    # gap (see ORACLE toggle above).
    for (j in seq_len(K0)) {
      a_mat[, j] <- mvtnorm::dmvnorm(Y, mean = MU0_LIST[[j]], sigma = SIGMA0_LIST[[j]])
    }
  } else {
    # RAW-MIXTURE phenotype densities a_ij, for the marginal-likelihood
    # sandwich correction (rawmix variant): rather than approximating each
    # merged phenotype as one Gaussian, reconstruct its true sub-mixture
    # density from the ORIGINAL (unmerged) raw fitted components mapped to
    # it, using the last iteration's raw component weights as a plug-in.
    z_raw_last <- class_mat[m, ]
    raw_w <- sapply(top_k, function(ell) sum(w[z_raw_last == ell]))
    for (j in seq_len(K0)) {
      raw_idx_j <- which(component_to_phenotype == j)
      raw_labels_j <- top_k[raw_idx_j]
      wj <- raw_w[raw_idx_j]
      if (sum(wj) <= 0) wj <- rep(1, length(wj)) # degenerate fallback
      wj <- wj / sum(wj)
      dens_j <- matrix(0, n_units, length(raw_labels_j))
      for (t in seq_along(raw_labels_j)) {
        ell <- raw_labels_j[t]
        dens_j[, t] <- mvtnorm::dmvnorm(Y, mean = fit$mu_k[[ell]], sigma = fit$Sigma_k[[ell]])
      }
      a_mat[, j] <- as.numeric(dens_j %*% wj)
    }
  }

  list(prev_draws = prev_draws_aligned, z_point_estimate = as.integer(z_point_aligned),
       K_active = K0, K_raw_active = K_active, a_mat = a_mat)
}

# RAW-MIXTURE Godambe sandwich correction: H and the (jackknife) score
# covariance both come from the MARGINAL mixture log-likelihood (see
# header), using a_mat (n x K0) -- the TRUE sub-mixture density per
# phenotype, reconstructed from the original unmerged raw fitted components
# (relabel_and_summarize's a_mat), not a single-Gaussian approximation.
sandwich_correct_prevalence <- function(prev_draws, w, psu_id, stratum, z_point,
                                         K_active, a_mat,
                                         apply_corr = APPLY_CORR) {
  pi_full <- colMeans(prev_draws)
  Kd <- K_active - 1

  m_i <- as.numeric(a_mat %*% pi_full)
  m_i <- pmax(m_i, 1e-300)
  a_diff <- a_mat[, seq_len(Kd), drop = FALSE] - a_mat[, K_active]

  # Soft score: w_i*(a_ik - a_iK)/m_i, replacing the hard one-hot score.
  unit_scores <- w * (a_diff / m_i)

  # H_kl = sum_i w_i*(a_ik-a_iK)*(a_il-a_iK)/m_i^2 -- consistent with the
  # score above (outer product identity).
  H_hat <- crossprod(a_diff * sqrt(w) / m_i)

  unique_psus <- unique(psu_id)
  n_psus <- length(unique_psus)
  psu_map <- setNames(seq_len(n_psus) - 1L, as.character(unique_psus))

  psu_scores <- matrix(0, n_psus, Kd)
  for (i in seq_along(psu_id)) {
    p_idx <- psu_map[[as.character(psu_id[i])]] + 1L
    psu_scores[p_idx, ] <- psu_scores[p_idx, ] + unit_scores[i, ]
  }

  strata_ids <- unique(stratum)
  stratum_psu_indices <- lapply(strata_ids, function(h) {
    h_psus <- unique(psu_id[stratum == h])
    as.integer(psu_map[as.character(h_psus)])
  })

  # Delete-one-PSU jackknife (closed-form Rao-Wu ultimate-cluster variance
  # for this linear score total -- see dpmm_helpers.cpp's derivation).
  J_hat <- jackknife_sandwich_cpp(psu_scores, stratum_psu_indices)

  H_inv <- solve(H_hat + diag(1e-8, Kd))
  Sigma_sandwich_raw <- H_inv %*% J_hat %*% H_inv

  # FINITE-SAMPLE DESIGN DF CORRECTION: with few PSUs relative to strata,
  # cluster-robust sandwich variance is known to be biased downward (see
  # e.g. Cameron & Miller 2015's review of few-cluster bias; this is the
  # same logic behind Stata's default cluster-robust scale factor). M is
  # the number of sampled PSUs and H the number of strata in THIS
  # replicate; df = M - H is the design degrees of freedom, used both to
  # scale up the sandwich piece (c_corr = M/df) and, downstream in
  # aggregate_fixedpop_results/test_a_centered_coverage, to widen the
  # quantile-based interval by t_{df,0.975}/z_0.975 instead of just z.
  M_psu <- length(unique(psu_id))
  H_strata <- length(unique(stratum))
  df_design <- max(M_psu - H_strata, 1)
  if (apply_corr) {
    c_corr <- M_psu / df_design
    t_crit <- qt(0.975, df_design)
  } else {
    c_corr <- 1
    t_crit <- qnorm(0.975)
  }
  Sigma_sandwich <- Sigma_sandwich_raw * c_corr

  draws_reduced <- prev_draws[, -K_active, drop = FALSE]
  draws_centered <- sweep(draws_reduced, 2, colMeans(draws_reduced))
  Sigma_naive <- cov(draws_reduced)

  # TOTAL-VARIANCE RECOMBINATION (law of total variance): the design-based
  # sandwich (classification/sampling variability given fixed mu_k,
  # Sigma_k) and the raw MCMC posterior covariance (which reflects
  # (mu_k,Sigma_k)-mixing/parameter uncertainty within one chain) are
  # treated as two additive, roughly independent variance components. The
  # rotation target becomes their SUM rather than Sigma_sandwich alone.
  # Sigma_sandwich here is already the df-corrected version above.
  Sigma_corrected <- Sigma_sandwich + Sigma_naive

  R1 <- tryCatch(chol(Sigma_corrected), error = function(e) chol(Sigma_corrected + diag(1e-6, Kd)))
  R2 <- tryCatch(chol(Sigma_naive), error = function(e) chol(Sigma_naive + diag(1e-6, Kd)))
  rotated <- draws_centered %*% solve(R2) %*% R1
  adjusted_draws <- sweep(rotated, 2, colMeans(draws_reduced), "+")

  list(pi_hat = colMeans(prev_draws), adjusted_draws = adjusted_draws,
       Sigma_sandwich = Sigma_sandwich, Sigma_sandwich_raw = Sigma_sandwich_raw,
       Sigma_naive = Sigma_naive, Sigma_corrected = Sigma_corrected,
       M_psu = M_psu, H_strata = H_strata, df_design = df_design,
       c_corr = c_corr, t_crit = t_crit)
}

run_one_sample <- function(sample_id) {
  samp_path <- file.path(SAMPLES_DIR, sprintf("sample_%04d.rds", sample_id))
  samp <- readRDS(samp_path)

  fit <- fit_wdpmm(samp$Y, samp$w)
  K_active_min <- K0 + 2
  rl <- relabel_and_summarize(fit, samp$w, samp$Y, K_active_min)
  sw <- sandwich_correct_prevalence(rl$prev_draws, samp$w, samp$psu_id,
                                     samp$stratum, rl$z_point_estimate, rl$K_active, rl$a_mat)

  list(pi_hat = sw$pi_hat, adjusted_draws = sw$adjusted_draws,
       K_active = rl$K_active, K_raw_active = rl$K_raw_active,
       true_prevalence = TRUE_PREVALENCE,
       M_psu = sw$M_psu, H_strata = sw$H_strata, df_design = sw$df_design,
       c_corr = sw$c_corr, t_crit = sw$t_crit,
       posterior = list(prev_draws_raw = rl$prev_draws, a_mat = rl$a_mat,
                         z_point_estimate = rl$z_point_estimate,
                         mu_k = fit$mu_k, Sigma_k = fit$Sigma_k,
                         Sigma_sandwich = sw$Sigma_sandwich, Sigma_sandwich_raw = sw$Sigma_sandwich_raw,
                         Sigma_naive = sw$Sigma_naive, Sigma_corrected = sw$Sigma_corrected,
                         M_psu = sw$M_psu, H_strata = sw$H_strata, df_design = sw$df_design,
                         c_corr = sw$c_corr, t_crit = sw$t_crit))
}

# Single shared results/ folder across the full grid (small per-replicate
# summaries); one posteriors_<ARM_TAG>/ folder per grid point for the
# heavier raw posterior objects (bulk du/tar/rm per arm), matching
# data/pop_c<v>/ and data/samples_<ARM_TAG>/.
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
POSTERIOR_DIR <- file.path(PROJECT_ROOT, sprintf("posteriors_%s", RESULT_TAG))
if (!dir.exists(POSTERIOR_DIR)) dir.create(POSTERIOR_DIR, recursive = TRUE)

run_and_save_sample <- function(sample_id) {
  out_path  <- file.path(RESULTS_DIR, sprintf("result_%s_%04d.rds", RESULT_TAG, sample_id))
  post_path <- file.path(POSTERIOR_DIR, sprintf("posterior_%04d.rds", sample_id))
  err_path  <- file.path(RESULTS_DIR, sprintf("ERROR_%s_%04d.rds", RESULT_TAG, sample_id))

  # Resumable: skip samples already completed (a prior ERROR is retried).
  if (file.exists(out_path)) {
    cat(sprintf("[SKIP] sample_id=%d already completed -> %s\n", sample_id, out_path))
    return(invisible(NULL))
  }

  t0 <- proc.time()
  result <- tryCatch(
    run_one_sample(sample_id),
    error = function(e) {
      saveRDS(list(sample_id = sample_id, cluster_c = CLUSTER_C, n_psu = N_CLUSTERS_SAMPLE,
                   sampling_frac = 0.01,
                   error = conditionMessage(e)), err_path)
      NULL
    }
  )
  elapsed <- (proc.time() - t0)["elapsed"]

  if (is.null(result)) {
    cat(sprintf("[FAILED] sample_id=%d elapsed=%.2fs -- see %s\n", sample_id, elapsed, err_path))
    quit(status = 1, save = "no")
  }

  posterior <- result$posterior
  posterior$sample_id <- sample_id
  posterior$cluster_c <- CLUSTER_C
  posterior$n_psu <- N_CLUSTERS_SAMPLE
  posterior$sampling_frac <- 0.01
  posterior$apply_corr <- APPLY_CORR
  saveRDS(posterior, post_path)
  result$posterior <- NULL

  result$sample_id <- sample_id
  result$cluster_c <- CLUSTER_C
  result$n_psu <- N_CLUSTERS_SAMPLE
  result$sampling_frac <- 0.01
  result$apply_corr <- APPLY_CORR
  result$elapsed_seconds <- elapsed
  saveRDS(result, out_path)
  cat(sprintf("[OK] sample_id=%d elapsed=%.2fs -> %s (posterior -> %s)\n",
              sample_id, elapsed, out_path, post_path))
}

aggregate_fixedpop_results <- function(cluster_c_str = CLUSTER_C_STR, n_psu_str = "4", arm_tag = NULL, apply_corr = FALSE) {
  if (is.null(arm_tag)) {
    if (nzchar(ARM_TAG_ENV)) {
      arm_tag <- ARM_TAG_ENV
    } else {
      arm_tag <- sprintf("c%s", cluster_c_str)
      if (n_psu_str != "4") arm_tag <- paste0(arm_tag, "_npsu", n_psu_str)
      if (apply_corr) arm_tag <- paste0(arm_tag, "_corr1")
    }
  }
  pattern <- sprintf("^result_%s_[0-9]+\\.rds$", gsub("\\.", "\\\\.", arm_tag))
  files <- list.files(RESULTS_DIR, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No completed result files found for arm ", arm_tag)

  n_files <- length(files)
  pi_mat <- matrix(NA, n_files, K0)
  true_prev_mat <- matrix(NA, n_files, K0)
  covered <- matrix(NA, n_files, K0)
  ci_widths <- matrix(NA, n_files, K0)
  K_raw_vec <- integer(n_files)
  df_vec <- numeric(n_files)
  c_corr_vec <- numeric(n_files)

  for (i in seq_along(files)) {
    r <- readRDS(files[i])
    pi_hat_i <- r$pi_hat[1:K0]
    pi_mat[i, ] <- pi_hat_i
    # Use THIS result's own stored true_prevalence (the population it was
    # actually fit against), not the session-global TRUE_PREVALENCE -- those
    # differ whenever this function is called for more than one CLUSTER_C
    # value within the same sourced R session.
    true_prev_i <- if (!is.null(r$true_prevalence)) r$true_prevalence[1:K0] else TRUE_PREVALENCE
    true_prev_mat[i, ] <- true_prev_i
    adj <- r$adjusted_draws[, seq_len(K0 - 1), drop = FALSE]
    last_col <- 1 - rowSums(adj)
    all_draws <- cbind(adj, last_col)
    ci_lo <- apply(all_draws, 2, quantile, probs = 0.025)
    ci_hi <- apply(all_draws, 2, quantile, probs = 0.975)

    # DESIGN-DF WIDTH CORRECTION: widen the empirical quantile interval by
    # t_{df,0.975}/z_0.975 (df = M_psu - H_strata for this replicate),
    # on top of the df-scaled Sigma_sandwich already baked into
    # adjusted_draws' spread via sandwich_correct_prevalence's c_corr.
    t_crit_i <- if (!is.null(r$t_crit)) r$t_crit else qnorm(0.975)
    width_scale <- t_crit_i / qnorm(0.975)
    ci_lo <- pi_hat_i - (pi_hat_i - ci_lo) * width_scale
    ci_hi <- pi_hat_i + (ci_hi - pi_hat_i) * width_scale

    ci_widths[i, ] <- ci_hi - ci_lo
    covered[i, ] <- (true_prev_i >= ci_lo) & (true_prev_i <= ci_hi)
    K_raw_vec[i] <- r$K_raw_active
    df_vec[i] <- if (!is.null(r$df_design)) r$df_design else NA_real_
    c_corr_vec[i] <- if (!is.null(r$c_corr)) r$c_corr else NA_real_
  }

  mean_pi <- colMeans(pi_mat)
  mean_true_prev <- colMeans(true_prev_mat)
  emp_se <- apply(pi_mat, 2, sd)
  mean_width <- colMeans(ci_widths)
  implied_se <- mean_width / (2 * 1.96)

  list(n_completed = n_files,
       true_prevalence = mean_true_prev,
       mean_pi_hat = mean_pi,
       bias = mean_pi - mean_true_prev,
       rel_bias_pct = 100 * (mean_pi - mean_true_prev) / mean_true_prev,
       empirical_se = emp_se,
       implied_se = implied_se,
       coverage_pct = colMeans(covered) * 100,
       mean_K_raw_active = mean(K_raw_vec),
       mean_df_design = mean(df_vec, na.rm = TRUE),
       mean_c_corr = mean(c_corr_vec, na.rm = TRUE))
}

# Test A (centered coverage): check whether the empirical replicate mean
# pi_bar (an essentially unbiased estimate of the truth, being an n=1000
# Monte Carlo average) falls inside each replicate's OWN quantile-based
# interval [ci_lo_m, ci_hi_m], instead of checking whether TRUE_PREVALENCE
# does. This tests purely whether each interval's WIDTH is calibrated to
# the actual replicate-to-replicate spread of pi_hat_m, independent of any
# systematic location bias of pi_hat_m relative to the true superpopulation
# value -- see rationale.tex's discussion of the coverage gap.
test_a_centered_coverage <- function(cluster_c_str = CLUSTER_C_STR, n_psu_str = "4", arm_tag = NULL, apply_corr = FALSE) {
  if (is.null(arm_tag)) {
    if (nzchar(ARM_TAG_ENV)) {
      arm_tag <- ARM_TAG_ENV
    } else {
      arm_tag <- sprintf("c%s", cluster_c_str)
      if (n_psu_str != "4") arm_tag <- paste0(arm_tag, "_npsu", n_psu_str)
      if (apply_corr) arm_tag <- paste0(arm_tag, "_corr1")
    }
  }
  pattern <- sprintf("^result_%s_[0-9]+\\.rds$", gsub("\\.", "\\\\.", arm_tag))
  files <- list.files(RESULTS_DIR, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No completed result files found for arm ", arm_tag)

  n_files <- length(files)
  pi_mat  <- matrix(NA, n_files, K0)
  true_prev_mat <- matrix(NA, n_files, K0)
  ci_lo_mat <- matrix(NA, n_files, K0)
  ci_hi_mat <- matrix(NA, n_files, K0)

  for (i in seq_along(files)) {
    r <- readRDS(files[i])
    pi_hat_i <- r$pi_hat[1:K0]
    pi_mat[i, ] <- pi_hat_i
    true_prev_mat[i, ] <- if (!is.null(r$true_prevalence)) r$true_prevalence[1:K0] else TRUE_PREVALENCE
    adj <- r$adjusted_draws[, seq_len(K0 - 1), drop = FALSE]
    last_col <- 1 - rowSums(adj)
    all_draws <- cbind(adj, last_col)
    ci_lo <- apply(all_draws, 2, quantile, probs = 0.025)
    ci_hi <- apply(all_draws, 2, quantile, probs = 0.975)

    # Same design-df width correction as aggregate_fixedpop_results.
    t_crit_i <- if (!is.null(r$t_crit)) r$t_crit else qnorm(0.975)
    width_scale <- t_crit_i / qnorm(0.975)
    ci_lo_mat[i, ] <- pi_hat_i - (pi_hat_i - ci_lo) * width_scale
    ci_hi_mat[i, ] <- pi_hat_i + (ci_hi - pi_hat_i) * width_scale
  }

  pi_bar <- colMeans(pi_mat)
  covered_a <- sweep(ci_lo_mat, 2, pi_bar, function(lo, pb) pb >= lo) &
               sweep(ci_hi_mat, 2, pi_bar, function(hi, pb) pb <= hi)

  list(n_completed = n_files,
       true_prevalence = colMeans(true_prev_mat),
       pi_bar = pi_bar,
       test_a_coverage_pct = colMeans(covered_a) * 100)
}

# Builds the per-phenotype summary table (estimate, bias, empirical/
# estimated SE and their ratio, standard + Test A coverage) for one grid
# arm and writes it to summaries/summary_<arm_tag>.csv, plus an RDS with
# the same table and top-level run metadata (n_completed, mean K_raw_active).
SUMMARIES_DIR <- file.path(PROJECT_ROOT, "summaries")

save_summary_table <- function(cluster_c_str = CLUSTER_C_STR, n_psu_str = "4", arm_tag = NULL, apply_corr = FALSE) {
  if (!dir.exists(SUMMARIES_DIR)) dir.create(SUMMARIES_DIR, recursive = TRUE)

  if (is.null(arm_tag)) {
    if (nzchar(ARM_TAG_ENV)) {
      arm_tag <- ARM_TAG_ENV
    } else {
      arm_tag <- sprintf("c%s", cluster_c_str)
      if (n_psu_str != "4") arm_tag <- paste0(arm_tag, "_npsu", n_psu_str)
      if (apply_corr) arm_tag <- paste0(arm_tag, "_corr1")
    }
  }

  agg <- aggregate_fixedpop_results(cluster_c_str, n_psu_str, arm_tag = arm_tag)
  ta  <- test_a_centered_coverage(cluster_c_str, n_psu_str, arm_tag = arm_tag)

  tbl <- data.frame(
    phenotype       = seq_len(K0),
    true_prevalence = agg$true_prevalence,
    estimate        = agg$mean_pi_hat,
    bias            = agg$bias,
    rel_bias_pct    = agg$rel_bias_pct,
    empirical_se    = agg$empirical_se,
    estimated_se    = agg$implied_se,
    se_ratio        = agg$empirical_se / agg$implied_se,
    coverage_pct    = agg$coverage_pct,
    test_a_coverage_pct = ta$test_a_coverage_pct
  )

  csv_path <- file.path(SUMMARIES_DIR, sprintf("summary_%s.csv", arm_tag))
  rds_path <- file.path(SUMMARIES_DIR, sprintf("summary_%s.rds", arm_tag))
  write.csv(tbl, csv_path, row.names = FALSE)
  saveRDS(list(arm_tag = arm_tag, n_completed = agg$n_completed,
               mean_K_raw_active = agg$mean_K_raw_active,
               mean_df_design = agg$mean_df_design, mean_c_corr = agg$mean_c_corr,
               table = tbl),
          rds_path)

  cat(sprintf("Summary for arm '%s' (n=%d, df=%.0f, c_corr=%.3f) written to %s\n",
              arm_tag, agg$n_completed, agg$mean_df_design, agg$mean_c_corr, csv_path))
  invisible(tbl)
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
