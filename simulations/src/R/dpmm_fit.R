# =============================================================================
# Fit the survey-weighted overfitted Gaussian mixture to one simulated sample
# (one SLURM array task) and save its calibrated phenotype prevalences.
#
#   1. fit_wdpmm              Gibbs sampler, L = 10 sparse components
#   2. summarize_phenotypes   merge components into phenotypes; prevalences
#   3. godambe_sandwich       design-based SEs for the prevalences
#
# Options (environment variables; non-default values add a RESULT_TAG suffix):
#   CLUSTER_C, N_PSU   scenario (0.3, 4)
#   MERGE              datadriven (default) | truth
#   THRESHOLD          1 (default) | 0   threshold consolidation
#   OUT_TAG            replaces RESULT_TAG (smoke tests)
#
# Summaries across replicates: summarize.R.
# =============================================================================

personal_lib <- "/nas/longleaf/home/aquijano/R/x86_64-pc-linux-gnu-library/4.5"
if (dir.exists(personal_lib) && !(personal_lib %in% .libPaths())) {
  .libPaths(c(personal_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
  library(Rcpp)
  library(RcppArmadillo)
})

# PROJECT_ROOT is simulations/, two levels above this file.
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1])))
  sys.frame(1)$ofile  # fall back to the source()-invocation case
}
SCRIPT_DIR <- tryCatch(dirname(get_script_path()), error = function(e) getwd())
PROJECT_ROOT <- dirname(dirname(SCRIPT_DIR))
cpp_file <- file.path(dirname(SCRIPT_DIR), "cpp", "dpmm_helpers.cpp")
cache_dir <- file.path(PROJECT_ROOT, ".rcpp_cache")
if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

# Array tasks share one compile cache; retry if another task holds it.
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

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------
CLUSTER_C_STR <- Sys.getenv("CLUSTER_C", unset = "0.3")
N_PSU_STR     <- Sys.getenv("N_PSU", unset = "4")
MERGE         <- match.arg(Sys.getenv("MERGE", unset = "datadriven"), c("datadriven", "truth"))
THRESHOLD     <- Sys.getenv("THRESHOLD", unset = "1") == "1" && MERGE == "datadriven"

# Scenario tag c<C>[_npsu<N>]; the default 4 PSUs adds no suffix.
SCENARIO <- sprintf("c%s", CLUSTER_C_STR)
if (N_PSU_STR != "4") SCENARIO <- paste0(SCENARIO, "_npsu", N_PSU_STR)

# Production is the bare scenario tag; a suffix names each departure from it.
RESULT_TAG <- paste0(SCENARIO,
                     if (MERGE == "truth") "_truth" else if (!THRESHOLD) "_nothresh")
OUT_TAG <- Sys.getenv("OUT_TAG", unset = "")
if (nzchar(OUT_TAG)) RESULT_TAG <- OUT_TAG

POP_PATH    <- file.path(PROJECT_ROOT, "data", sprintf("pop_c%s", CLUSTER_C_STR), "population.rds")
SAMPLES_DIR <- file.path(PROJECT_ROOT, "data", sprintf("samples_%s", SCENARIO))

# -----------------------------------------------------------------------------
# Fixed design constants
# -----------------------------------------------------------------------------
L_TRUNC       <- 10
ALPHA0_SPARSE <- 1
SIGMA_FLOOR   <- 0.05
N_MCMC_ITER   <- 1500
N_BURNIN      <- 500
FIT_BASE_SEED <- 20260930L   # replicate s uses FIT_BASE_SEED + s
MERGE_MIN_WEIGHT <- 0.02     # threshold: minimum share of total weight
MERGE_MIN_SEP    <- 1.5      # threshold: minimum centroid distance

pop <- readRDS(POP_PATH)
K0 <- length(pop$MU0_LIST)
P_DIM <- ncol(pop$Y)
MU0_LIST <- pop$MU0_LIST
TRUE_PREVALENCE <- pop$true_prevalence

# -----------------------------------------------------------------------------
# 1. Sampler
# -----------------------------------------------------------------------------
floor_eigen <- function(Sigma, floor = SIGMA_FLOOR) {
  eig <- eigen(Sigma, symmetric = TRUE)
  vals <- pmax(eig$values, floor)
  eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
}

# Gibbs sampler: unweighted allocation; weights (summing to n) enter the
# Dirichlet counts and the NIW sufficient statistics. Returns the kept
# allocations.
fit_wdpmm <- function(Y, w, L = L_TRUNC, n_iter = N_MCMC_ITER, n_burnin = N_BURNIN,
                       alpha0 = ALPHA0_SPARSE, kappa0 = 0.5, nu0 = P_DIM + 2) {
  n <- nrow(Y); p <- ncol(Y)
  mu0 <- colSums(Y * w) / sum(w)
  Psi0 <- diag(colSums(w * sweep(Y, 2, mu0)^2) / sum(w)) * (nu0 - p - 1)

  z <- sample.int(L, n, replace = TRUE)
  mu_k <- lapply(seq_len(L), function(k) mu0 + rnorm(p, 0, 1))
  Sigma_k <- lapply(seq_len(L), function(k) diag(p))

  keep <- (n_burnin + 1):n_iter
  z_draws <- matrix(NA_integer_, length(keep), n)
  draw_i <- 0
  n_dens_fallback <- 0L; n_sigma_fallback <- 0L

  for (iter in seq_len(n_iter)) {
    pi0 <- sparse_dirichlet_weights_cpp(z, w, as.integer(L), alpha0)

    logp <- matrix(-Inf, n, L)
    for (k in seq_len(L)) {
      dens <- tryCatch(
        mvtnorm::dmvnorm(Y, mean = mu_k[[k]], sigma = Sigma_k[[k]], log = TRUE),
        error = function(e) { n_dens_fallback <<- n_dens_fallback + 1L; rep(-1e10, n) }
      )
      logp[, k] <- log(max(pi0[k], 1e-12)) + dens
    }
    logp <- logp - apply(logp, 1, max)
    probs <- exp(logp); probs <- probs / rowSums(probs)
    z <- sample_z_cpp(probs)

    for (k in seq_len(L)) {
      idx <- which(z == k)
      if (length(idx) == 0) {
        # Empty component: draw from the NIW prior.
        Sigma_k[[k]] <- floor_eigen(solve(rWishart(1, nu0, solve(Psi0))[, , 1]))
        mu_k[[k]] <- MASS::mvrnorm(1, mu0, Sigma_k[[k]] / kappa0)
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
        error = function(e) { n_sigma_fallback <<- n_sigma_fallback + 1L; diag(p) }
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
  list(z_draws = z_draws,
       n_dens_fallback = n_dens_fallback, n_sigma_fallback = n_sigma_fallback)
}

# -----------------------------------------------------------------------------
# 2. Phenotypes and prevalences
# -----------------------------------------------------------------------------
weighted_centroid <- function(Y, w, idx) {
  if (length(idx) == 0) return(rep(1e6, ncol(Y)))  # empty: never the nearest
  colSums(Y[idx, , drop = FALSE] * w[idx]) / sum(w[idx])
}

# Threshold consolidation: repeatedly merge the lightest group into its
# nearest neighbour if it is below min_weight_frac of the total weight, else
# the closest pair of groups if they are less than min_sep apart.
merge_small_close_groups <- function(z, Y, w, min_weight_frac, min_sep) {
  total_w <- sum(w)
  repeat {
    groups <- sort(unique(z))
    K <- length(groups)
    if (K <= 1) break
    wt  <- sapply(groups, function(g) sum(w[z == g]))
    cen <- t(sapply(groups, function(g) weighted_centroid(Y, w, which(z == g))))
    dmat <- as.matrix(dist(cen))
    diag(dmat) <- Inf

    smallest_idx <- which.min(wt)
    if (wt[smallest_idx] / total_w < min_weight_frac) {
      g_small  <- groups[smallest_idx]
      g_target <- groups[which.min(dmat[smallest_idx, ])]
      z[z == g_small] <- g_target
      next
    }

    min_pair <- which(dmat == min(dmat), arr.ind = TRUE)[1, ]
    if (dmat[min_pair[1], min_pair[2]] < min_sep) {
      z[z == groups[min_pair[2]]] <- groups[min_pair[1]]
      next
    }
    break
  }
  match(z, sort(unique(z)))
}

# Label G estimated groups by true phenotype, for scoring: each phenotype
# takes its nearest unassigned group; leftover groups join their nearest
# phenotype.
align_groups_to_truth <- function(centroids) {
  G <- nrow(centroids)
  cost <- sapply(seq_len(K0), function(j) sqrt(colSums((t(centroids) - MU0_LIST[[j]])^2)))
  cost <- matrix(cost, G, K0)
  assigned <- integer(K0)
  avail <- seq_len(G)
  for (j in seq_len(K0)) {
    if (length(avail) == 0) avail <- seq_len(G)  # G < K0: allow reuse
    best <- avail[which.min(cost[avail, j])]
    assigned[j] <- best
    avail <- setdiff(avail, best)
  }
  group_to_pheno <- integer(G)
  group_to_pheno[assigned] <- seq_len(K0)
  leftover <- setdiff(seq_len(G), assigned)
  if (length(leftover) > 0) {
    group_to_pheno[leftover] <- apply(cost[leftover, , drop = FALSE], 1, which.min)
  }
  group_to_pheno
}

# Merge the raw components into phenotypes and summarize the draws.
#   datadriven: the units are partitioned into groups from the posterior
#               similarity matrix (plus threshold consolidation); in every
#               draw each component joins the group holding most of its weight.
#   truth:      in every draw each component joins the nearest true centroid.
# Each draw then gives every unit a phenotype, from which come the draw's
# prevalences and each phenotype's weighted mean and covariance. Returns the
# prevalence draws and the posterior means (pi, mu, Sigma) per phenotype.
summarize_phenotypes <- function(fit, w, Y) {
  z_draws <- fit$z_draws
  m <- nrow(z_draws); n <- ncol(z_draws); p <- ncol(Y)

  if (MERGE == "datadriven") {
    z_hat <- as.integer(mcclust::maxpear(mcclust::comp.psm(z_draws))$cl)
    z_hat <- match(z_hat, sort(unique(z_hat)))
    if (THRESHOLD) {
      z_hat <- merge_small_close_groups(z_hat, Y, w, MERGE_MIN_WEIGHT, MERGE_MIN_SEP)
    }
    K_hat <- max(z_hat)
    group_centroids <- t(sapply(seq_len(K_hat), function(g) weighted_centroid(Y, w, which(z_hat == g))))
    group_to_pheno <- align_groups_to_truth(matrix(group_centroids, K_hat))
    group_w <- w * outer(z_hat, seq_len(K_hat), "==")
  } else {
    K_hat <- NA_integer_
    truth_mu <- do.call(rbind, MU0_LIST)
  }

  wY  <- w * Y
  a   <- rep(seq_len(p), p); b <- rep(seq_len(p), each = p)
  wYY <- w * Y[, a] * Y[, b]
  prev_draws <- matrix(0, m, K0)
  pheno_draws <- matrix(0L, m, n)
  mu_sum <- matrix(0, K0, p); S_sum <- matrix(0, K0, p * p); n_present <- integer(K0)

  for (i in seq_len(m)) {
    z <- z_draws[i, ]
    # Phenotype of each component present in this draw (rows in label order).
    if (MERGE == "datadriven") {
      comp_pheno <- group_to_pheno[max.col(rowsum(group_w, z), ties.method = "first")]
    } else {
      cen <- rowsum(wY, z) / as.numeric(rowsum(w, z))
      comp_pheno <- apply(cen, 1, function(x) which.min(colSums((t(truth_mu) - x)^2)))
    }
    ph <- comp_pheno[match(z, sort(unique(z)))]
    pheno_draws[i, ] <- ph

    s0 <- rowsum(w, ph)
    js <- as.integer(rownames(s0)); s0 <- as.numeric(s0)
    mu <- rowsum(wY, ph) / s0
    S  <- rowsum(wYY, ph) / s0 - mu[, a, drop = FALSE] * mu[, b, drop = FALSE]
    prev_draws[i, js] <- s0 / sum(w)
    mu_sum[js, ] <- mu_sum[js, ] + mu
    S_sum[js, ]  <- S_sum[js, ] + S
    n_present[js] <- n_present[js] + 1L
  }

  present <- which(n_present > 0)
  occ <- sort(table(z_draws), decreasing = TRUE)
  list(prev_draws = prev_draws,
       pi = colMeans(prev_draws), present = present,
       mu = lapply(present, function(j) mu_sum[j, ] / n_present[j]),
       Sigma = lapply(present, function(j) matrix(S_sum[j, ] / n_present[j], p, p)),
       z_point_estimate = apply(pheno_draws, 2, function(x) which.max(tabulate(x, K0))),
       K_hat_estimated = K_hat,
       K_raw_active = min(which(cumsum(occ) / sum(occ) >= 0.99)))
}

# -----------------------------------------------------------------------------
# 3. Design-based variance
# -----------------------------------------------------------------------------
# eta = (pi_1..pi_{G-1}, then per g: mu_g, vech(Sigma_g)), pi_G = 1 - sum.
pack_eta <- function(pi, mu, Sigma) {
  p <- length(mu[[1]])
  lt <- which(lower.tri(diag(p), diag = TRUE))
  G <- length(pi)
  c(pi[-G], unlist(lapply(seq_len(G), function(g) c(mu[[g]], Sigma[[g]][lt]))))
}

# Unit scores w_i grad_eta log m_i(eta) (n x d):
#   pi_k:      (phi_ik - phi_iG)/m_i
#   mu_g:      r_ig Sigma_g^{-1}(Y_i - mu_g)
#   Sigma_ab:  r_ig 0.5[D_a D_b - (Sigma_g^{-1})_ab], doubled off-diagonal,
#              with D = Sigma_g^{-1}(Y_i - mu_g)
unit_scores <- function(eta, Y, w, G) {
  n <- nrow(Y); p <- ncol(Y)
  lt <- which(lower.tri(diag(p), diag = TRUE)); q <- length(lt)
  lt_row <- row(diag(p))[lt]; lt_col <- col(diag(p))[lt]
  vech_mult <- ifelse(lt_row == lt_col, 1, 2)

  pi_g <- c(eta[seq_len(G - 1)], 1 - sum(eta[seq_len(G - 1)]))
  pos <- G - 1
  mu <- Sig <- vector("list", G)
  for (g in seq_len(G)) {
    mu[[g]] <- eta[pos + seq_len(p)]; pos <- pos + p
    S <- matrix(0, p, p); S[lt] <- eta[pos + seq_len(q)]; pos <- pos + q
    Sig[[g]] <- S + t(S) - diag(diag(S), p)
  }

  logphi <- sapply(seq_len(G), function(g) mvtnorm::dmvnorm(Y, mu[[g]], Sig[[g]], log = TRUE))
  lp <- sweep(logphi, 2, log(pi_g), "+")
  lmax <- apply(lp, 1, max)
  log_m <- lmax + log(rowSums(exp(lp - lmax)))
  R <- exp(lp - log_m)
  A <- exp(logphi - log_m)

  blocks <- list(A[, -G, drop = FALSE] - A[, G])
  for (g in seq_len(G)) {
    Sinv <- solve(Sig[[g]])
    D <- sweep(Y, 2, mu[[g]]) %*% Sinv
    s_S <- 0.5 * (D[, lt_row, drop = FALSE] * D[, lt_col, drop = FALSE] -
                    matrix(Sinv[lt], n, q, byrow = TRUE))
    blocks[[length(blocks) + 1]] <- R[, g] * D
    blocks[[length(blocks) + 1]] <- R[, g] * sweep(s_S, 2, vech_mult, "*")
  }
  do.call(cbind, blocks) * w
}

# Godambe sandwich V = H^{-1} J H^{-1} for eta = (pi, mu, Sigma) of the
# weighted Gaussian mixture with one component per phenotype, evaluated at
# the posterior means in ph (from summarize_phenotypes). H is the numerical
# Jacobian of the score, J the delete-one-PSU jackknife of the PSU score
# totals. Returns the K x K prevalence covariance and SEs, and the t critical
# value on M - H degrees of freedom (sampled PSUs minus strata).
godambe_sandwich <- function(Y, w, ph, psu_id, stratum, K) {
  keep <- ph$present
  G <- length(keep)
  if (G < 2) stop("fewer than two phenotypes in the partition")

  eta <- pack_eta(ph$pi[keep], ph$mu, ph$Sigma)
  d <- length(eta)

  score_total <- function(e) colSums(unit_scores(e, Y, w, G))
  jac <- matrix(0, d, d)
  for (j in seq_len(d)) {
    h <- 1e-5 * max(1, abs(eta[j]))
    ep <- eta; ep[j] <- ep[j] + h
    em_ <- eta; em_[j] <- em_[j] - h
    jac[, j] <- (score_total(ep) - score_total(em_)) / (2 * h)
  }
  H_hat <- -(jac + t(jac)) / 2

  # PSU score totals, and each stratum's PSU rows (0-based).
  unique_psus <- unique(psu_id)
  psu_scores <- rowsum(unit_scores(eta, Y, w, G), psu_id, reorder = FALSE)
  stratum_psu_indices <- lapply(unique(stratum), function(h) {
    match(unique(psu_id[stratum == h]), unique_psus) - 1L
  })
  J_hat <- jackknife_sandwich_cpp(psu_scores, stratum_psu_indices)

  H_inv <- solve(H_hat)
  V_eta <- H_inv %*% J_hat %*% H_inv
  df_design <- max(length(unique_psus) - length(unique(stratum)), 1)

  # Prevalence block -> all G components (pi_G = 1 - sum) -> K phenotypes.
  B <- rbind(diag(G - 1), -1)
  V_group <- B %*% V_eta[seq_len(G - 1), seq_len(G - 1), drop = FALSE] %*% t(B)
  A_map <- matrix(0, K, G)
  A_map[cbind(keep, seq_len(G))] <- 1
  V_pi <- A_map %*% V_group %*% t(A_map)

  list(V_pi = V_pi, se = sqrt(diag(V_pi)), V_eta = V_eta,
       df_design = df_design, t_crit = qt(0.975, df_design))
}

# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------
run_one_sample <- function(sample_id) {
  samp <- readRDS(file.path(SAMPLES_DIR, sprintf("sample_%04d.rds", sample_id)))

  # Same seed in every arm, so arms share the chain for each sample.
  set.seed(FIT_BASE_SEED + sample_id)
  fit <- fit_wdpmm(samp$Y, samp$w)
  ph <- summarize_phenotypes(fit, samp$w, samp$Y)
  gs <- godambe_sandwich(samp$Y, samp$w, ph, samp$psu_id, samp$stratum, K0)

  list(pi_hat = ph$pi, se = gs$se, t_crit = gs$t_crit,
       V_pi = gs$V_pi, df_design = gs$df_design,
       K_raw_active = ph$K_raw_active, K_hat_estimated = ph$K_hat_estimated,
       true_prevalence = TRUE_PREVALENCE,
       n_dens_fallback = fit$n_dens_fallback, n_sigma_fallback = fit$n_sigma_fallback,
       posterior = list(prev_draws = ph$prev_draws, z_point_estimate = ph$z_point_estimate,
                        mu = ph$mu, Sigma = ph$Sigma, V_eta = gs$V_eta))
}

# results/: one small file per replicate and arm; posteriors_<tag>/: the
# larger posterior objects.
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
POSTERIOR_DIR <- file.path(PROJECT_ROOT, sprintf("posteriors_%s", RESULT_TAG))

run_and_save_sample <- function(sample_id) {
  dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(POSTERIOR_DIR, showWarnings = FALSE, recursive = TRUE)
  out_path  <- file.path(RESULTS_DIR, sprintf("result_%s_%04d.rds", RESULT_TAG, sample_id))
  post_path <- file.path(POSTERIOR_DIR, sprintf("posterior_%04d.rds", sample_id))
  err_path  <- file.path(RESULTS_DIR, sprintf("ERROR_%s_%04d.rds", RESULT_TAG, sample_id))

  # Skip completed samples; failed ones are retried.
  if (file.exists(out_path)) {
    cat(sprintf("[SKIP] sample_id=%d already completed -> %s\n", sample_id, out_path))
    return(invisible(NULL))
  }

  settings <- list(result_tag = RESULT_TAG, scenario = SCENARIO, merge = MERGE,
                   threshold = THRESHOLD,
                   seed = FIT_BASE_SEED + sample_id)
  t0 <- proc.time()
  result <- tryCatch(
    run_one_sample(sample_id),
    error = function(e) {
      saveRDS(c(list(sample_id = sample_id, error = conditionMessage(e),
                     call = deparse(conditionCall(e))), settings), err_path)
      NULL
    }
  )
  elapsed <- (proc.time() - t0)["elapsed"]

  if (is.null(result)) {
    cat(sprintf("[FAILED] sample_id=%d elapsed=%.2fs -- see %s\n", sample_id, elapsed, err_path))
    quit(status = 1, save = "no")
  }

  saveRDS(c(result$posterior, list(sample_id = sample_id), settings), post_path)
  result$posterior <- NULL
  saveRDS(c(result, list(sample_id = sample_id, elapsed_seconds = elapsed), settings), out_path)
  cat(sprintf("[OK] sample_id=%d elapsed=%.2fs -> %s\n", sample_id, elapsed, out_path))
}

slurm_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
if (!is.na(slurm_task_id)) {
  run_and_save_sample(as.integer(slurm_task_id))
} else {
  cat("Sourced without SLURM_ARRAY_TASK_ID -- defined functions only.\n")
}
