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
#
# Defaults below are set for a fast, runnable demo (n ~ 2000, R = 50
# replicates). For the real study, scale N_POP / N_PSU up and set
# R_REPLICATES = 500, N_MCMC_ITER = 5000+, per the design decided in
# the accompanying Methods write-up.
# =============================================================================

suppressPackageStartupMessages({
  library(MASS)             # mvrnorm
  library(mvtnorm)          # dmvnorm
  library(label.switching)  # Stephens (2000) relabeling algorithm
  library(coda)             # Rhat / effective sample size diagnostics
})

set.seed(20260918)

# -----------------------------------------------------------------------
# 0. GLOBAL SETTINGS
# -----------------------------------------------------------------------

N_POP          <- 2000000 # population size (before sampling) -- kept >> TARGET_N
                           # so the sampling fraction f = TARGET_N/N_POP stays
                           # near 0, consistent with the WOR/WR-equivalence
                           # simplification adopted for this method (no finite-
                           # population correction needed; ordinary PSU-level
                           # with-replacement bootstrap is used for J_hat below)
N_PSU          <- 200     # number of PSUs in the population (PPS3's z2 layer)
                           # -> ~10,000 population elements per PSU on average
STRATA_H       <- 10      # fixed number of strata; PSUs are nested within
                           # strata and sampled independently per stratum
                           # (proportional allocation). Purely a sampling-
                           # design partition -- does not enter the outcome
                           # model, consistent with strata being classification
                           # variables known at design time.
P_DIM          <- 3       # dimension of the CKM-marker-like outcome
K0             <- 4       # true number of phenotypes
TARGET_N       <- 2000    # target sample size after PPS selection
                           # f = TARGET_N / N_POP = 0.001 here. For the real
                           # study run, scale both up together (e.g. N_POP =
                           # 20,000,000 / TARGET_N = 16,415) to keep f this
                           # small rather than shrinking TARGET_N.
L_TRUNC        <- 10      # DPMM truncation level
N_MCMC_ITER    <- 1500
N_BURNIN       <- 500
SIGMA_FLOOR    <- 0.05    # minimum eigenvalue floor for Sigma_k
R_REPLICATES   <- 50      # Monte Carlo replicates (set to 500 for the real run)
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
# True marginal prevalence (population average over the x2-informative
# membership model, approximated once by simulation -- see get_true_prevalence())

# -----------------------------------------------------------------------
# 1. GENERATIVE MODEL: PPS1 (informative selection) + PPS3 (PSU effect)
#    extended to a K0-component mixture outcome
# -----------------------------------------------------------------------

#' Simulate one finite population.
#' x2 drives BOTH true phenotype membership (informativeness on the
#' outcome) and PPS selection probability -- the PPS1 mechanism.
#' z2 is an independent PSU-level effect entering the outcome only,
#' never the design -- the PPS3 mechanism.
generate_population <- function(N = N_POP, n_psu = N_PSU, K0 = K0,
                                 H = STRATA_H, beta_x2 = 1.2, z2_mean = 5,
                                 include_z2 = TRUE) {

  psu_id <- sample(seq_len(n_psu), N, replace = TRUE)
  # PSUs nested within strata: assign once at the PSU level, then broadcast
  # to elements. Roughly equal-sized strata by PSU count.
  psu_stratum <- rep(seq_len(H), length.out = n_psu)[sample(n_psu)]
  stratum <- psu_stratum[psu_id]

  z2_psu <- if (include_z2) rexp(n_psu, rate = 1 / z2_mean) else rep(0, n_psu)
  z2     <- z2_psu[psu_id] - z2_mean * include_z2  # centered

  x1 <- rnorm(N)
  x2 <- rgamma(N, shape = 2, rate = 1)             # positive "size" auxiliary

  # Multinomial logit membership on x2: informativeness acts on COMPOSITION
  alpha_k <- c(0, 0.5, -0.5, 0.2)                  # baseline = component 1
  beta_k  <- c(0, beta_x2, -beta_x2 / 2, beta_x2 / 3)
  eta     <- outer(x2, beta_k) + matrix(alpha_k, N, K0, byrow = TRUE)
  probs   <- exp(eta) / rowSums(exp(eta))
  z_true  <- apply(probs, 1, function(p) sample.int(K0, 1, prob = p))

  Y <- matrix(NA_real_, N, P_DIM)
  for (k in seq_len(K0)) {
    idx <- which(z_true == k)
    if (length(idx) > 0) {
      Y[idx, ] <- MASS::mvrnorm(length(idx), mu = MU0_LIST[[k]], Sigma = SIGMA0_LIST[[k]])
      Y[idx, ] <- Y[idx, ] + z2[idx]  # additive PSU-level shift, PPS3 mechanism
    }
  }

  # PPS inclusion probability proportional to x2 (informative selection)
  pi_raw <- x2 / sum(x2)

  list(Y = Y, x1 = x1, x2 = x2, psu_id = psu_id, stratum = stratum,
       z_true = z_true, pi_raw = pi_raw)
}

#' Stratified Poisson (approx PPS) sampling to target sample size n.
#' Allocation across strata is proportional to stratum population size
#' (standard proportional allocation); within each stratum, PPS selection
#' proceeds independently on that stratum's own x2-based inclusion
#' probabilities. Returns sampled indices, inclusion probabilities, and
#' normalized weights, along with element-level stratum and PSU labels
#' needed for the stratified PSU-bootstrap in the sandwich correction.
pps_sample <- function(pop, n_target, H = STRATA_H) {
  idx_all <- integer(0); w_all <- numeric(0); pi_all <- numeric(0)
  N_total <- length(pop$stratum)

  for (h in seq_len(H)) {
    idx_h <- which(pop$stratum == h)
    Nh <- length(idx_h)
    if (Nh == 0) next
    n_h_target <- max(1, round(n_target * Nh / N_total))  # proportional allocation

    pi_raw_h <- pop$x2[idx_h] / sum(pop$x2[idx_h])
    pi_i_h <- pmin(pi_raw_h * n_h_target, 0.999)
    included_h <- rbinom(length(pi_i_h), 1, pi_i_h) == 1

    idx_all <- c(idx_all, idx_h[included_h])
    w_all   <- c(w_all, 1 / pi_i_h[included_h])
    pi_all  <- c(pi_all, pi_i_h[included_h])
  }

  w_all <- w_all / sum(w_all) * length(idx_all)  # normalize to sum to sample size
  list(idx = idx_all, w = w_all, pi_i = pi_all,
       Y = pop$Y[idx_all, , drop = FALSE], psu_id = pop$psu_id[idx_all],
       stratum = pop$stratum[idx_all],
       z_true = pop$z_true[idx_all], x2 = pop$x2[idx_all])
}

#' Monte Carlo approximation of the true marginal prevalence pi_0,
#' i.e., E_{x2}[softmax membership probabilities], computed once.
get_true_prevalence <- function(n_mc = 200000) {
  x2 <- rgamma(n_mc, shape = 2, rate = 1)
  alpha_k <- c(0, 0.5, -0.5, 0.2)
  beta_k  <- c(0, 1.2, -0.6, 0.4)
  eta   <- outer(x2, beta_k) + matrix(alpha_k, n_mc, K0, byrow = TRUE)
  probs <- exp(eta) / rowSums(exp(eta))
  colMeans(probs)
}

# -----------------------------------------------------------------------
# 2. SURVEY-WEIGHTED TRUNCATED STICK-BREAKING DPMM (blocked Gibbs sampler)
# -----------------------------------------------------------------------

#' Enforce a minimum-eigenvalue floor on a covariance matrix.
floor_eigen <- function(Sigma, floor = SIGMA_FLOOR) {
  eig <- eigen(Sigma, symmetric = TRUE)
  vals <- pmax(eig$values, floor)
  eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
}

#' Fit the weighted truncated DPMM.
#' include_x2: if TRUE, membership logits condition on x2 (the "informed"
#' fitted model); if FALSE, weights are fit with a single global stick-
#' breaking prior only (the "naive," design-omitting fitted model) --
#' this is the direct analogue of PPS1's "analyst fits mu = f(x1) only."
fit_wdpmm <- function(Y, w, x2 = NULL, include_x2 = FALSE,
                       L = L_TRUNC, n_iter = N_MCMC_ITER, n_burnin = N_BURNIN,
                       a_alpha = 1, b_alpha = 1,
                       kappa0 = 0.5, nu0 = P_DIM + 2) {

  n <- nrow(Y); p <- ncol(Y)
  mu0 <- colMeans(Y)
  Psi0 <- diag(apply(Y, 2, var)) * (nu0 - p - 1)

  # init
  z <- sample.int(L, n, replace = TRUE)
  mu_k <- lapply(seq_len(L), function(k) mu0 + rnorm(p, 0, 1))
  Sigma_k <- lapply(seq_len(L), function(k) diag(p))
  alpha <- 1
  v <- rep(1 / L, L)

  # storage for post-burn-in draws
  keep <- (n_burnin + 1):n_iter
  z_draws <- matrix(NA_integer_, length(keep), n)
  min_eig_trace <- numeric(n_iter)
  draw_i <- 0

  for (iter in seq_len(n_iter)) {

    # --- stick-breaking weights pi_k (global, unless include_x2) ---
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
      # covariate-dependent stick-breaking: logistic regression of each
      # stick on x2, fit via a simple weighted GLM per stick (approx.
      # Gibbs update -- a practical stand-in for full Polya-Gamma
      # augmentation, adequate for this validation exercise)
      pi_mat <- matrix(0, n, L)
      remaining <- rep(1, n)
      for (k in seq_len(L - 1)) {
        ind_k  <- as.numeric(z == k)
        ind_ge <- as.numeric(z >= k)
        dat <- data.frame(y = ind_k[ind_ge == 1], x2 = x2[ind_ge == 1],
                           w = w[ind_ge == 1])
        if (length(unique(dat$y)) < 2 || nrow(dat) < 5) {
          nu_k <- rep(mean(ind_k[ind_ge == 1] + 1e-3), n)
        } else {
          fit <- suppressWarnings(glm(y ~ x2, data = dat, weights = w,
                                       family = binomial()))
          nu_k <- plogis(predict(fit, newdata = data.frame(x2 = x2)))
        }
        pi_mat[, k] <- nu_k * remaining
        remaining <- remaining * (1 - nu_k)
      }
      pi_mat[, L] <- remaining
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
    z <- apply(probs, 1, function(p) sample.int(L, 1, prob = p))

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
      Sk <- matrix(0, p, p)
      for (j in seq_along(idx)) {
        d <- Yk[j, ] - ybar_k
        Sk <- Sk + wk[j] * outer(d, d)
      }
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

    # --- concentration parameter alpha | v (conjugate Gamma update) ---
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

#' Relabel MCMC draws of z and return the relabeled weighted prevalence
#' vector for each retained draw, plus a single Dahl-style point estimate
#' of the hard assignment (posterior co-clustering matrix, least-squares).
relabel_and_summarize <- function(fit, w, K_active) {
  z_draws <- fit$z_draws
  # collapse to K_active most frequently occupied components for stable
  # relabeling (standard practical step before calling label.switching)
  occ_freq <- table(z_draws)
  top_k <- as.integer(names(sort(occ_freq, decreasing = TRUE))[seq_len(K_active)])
  remap <- setNames(seq_along(top_k), top_k)
  z_collapsed <- matrix(sapply(z_draws, function(x) {
    if (as.character(x) %in% names(remap)) remap[as.character(x)] else NA
  }), nrow(z_draws), ncol(z_draws))

  valid_rows <- which(rowSums(is.na(z_collapsed)) == 0)
  z_collapsed <- z_collapsed[valid_rows, , drop = FALSE]

  # Stephens (2000) relabeling via label.switching package needs
  # per-draw, per-unit allocation probabilities; approximate with a
  # one-hot encoding (adequate for hard-assignment relabeling here)
  m <- nrow(z_collapsed); n <- ncol(z_collapsed)
  p_array <- array(0, dim = c(m, n, K_active))
  for (i in seq_len(m)) for (k in seq_len(K_active))
    p_array[i, , k] <- as.numeric(z_collapsed[i, ] == k)

  ls <- tryCatch(
    label.switching::stephens(p_array),
    error = function(e) NULL
  )
  if (is.null(ls)) {
    perm <- matrix(rep(seq_len(K_active), m), m, K_active, byrow = TRUE)
  } else {
    perm <- ls$permutations
  }

  prev_draws <- matrix(NA_real_, m, K_active)
  z_relabeled_last <- NULL
  for (i in seq_len(m)) {
    map_i <- perm[i, ]
    z_i <- z_collapsed[i, ]
    z_relab <- match(z_i, map_i)
    tab <- sapply(seq_len(K_active), function(k) sum(w[z_relab == k]))
    prev_draws[i, ] <- tab / sum(w)
    if (i == m) z_relabeled_last <- z_relab
  }

  list(prev_draws = prev_draws, z_point_estimate = z_relabeled_last)
}

# -----------------------------------------------------------------------
# 4. GODAMBE SANDWICH CORRECTION ON THE PREVALENCE FUNCTIONAL
#    (PSU-level bootstrap for J_hat; multinomial Fisher info for H_hat)
# -----------------------------------------------------------------------

sandwich_correct_prevalence <- function(prev_draws, w, psu_id, stratum, z_point,
                                         K_active, B = B_BOOTSTRAP) {

  pi_bar <- colMeans(prev_draws)[-K_active]   # drop last category, identifiability
  Kd <- length(pi_bar)

  # H_hat: multinomial Fisher information, scaled by total weight
  H_hat <- (diag(pi_bar) - outer(pi_bar, pi_bar)) * sum(w)

  # J_hat via STRATIFIED PSU-level bootstrap of the weighted score total:
  # resample PSUs with replacement independently WITHIN each stratum
  # (the standard ultimate-cluster / Rao-Wu-type multistage bootstrap
  # structure), then pool the resampled units' scores across strata for
  # each bootstrap replicate. Pooling PSUs across strata before resampling
  # would misrepresent a stratified design's variance.
  one_hot <- sapply(seq_len(Kd), function(k) as.numeric(z_point == k))
  strata_ids <- unique(stratum)
  psus_by_stratum <- lapply(strata_ids, function(h) unique(psu_id[stratum == h]))

  score_boot <- matrix(NA_real_, B, Kd)
  for (b in seq_len(B)) {
    sel_all <- integer(0)
    for (s in seq_along(strata_ids)) {
      psus_h <- psus_by_stratum[[s]]
      if (length(psus_h) == 0) next
      boot_psu_h <- sample(psus_h, length(psus_h), replace = TRUE)
      sel_h <- unlist(lapply(boot_psu_h, function(pu) {
        which(psu_id == pu & stratum == strata_ids[s])
      }))
      sel_all <- c(sel_all, sel_h)
    }
    score_boot[b, ] <- colSums(w[sel_all] * (one_hot[sel_all, , drop = FALSE] -
                                 matrix(pi_bar, length(sel_all), Kd, byrow = TRUE)))
  }
  J_hat <- cov(score_boot)

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

  list(pi_hat = pi_bar, Sigma_sandwich = Sigma_sandwich,
       adjusted_draws = adjusted_draws)
}

# -----------------------------------------------------------------------
# 5. ONE FULL REPLICATE: generate -> sample -> fit -> relabel -> correct
# -----------------------------------------------------------------------

run_one_replicate <- function(include_x2_in_fit = TRUE, include_z2_in_truth = TRUE) {
  pop    <- generate_population(include_z2 = include_z2_in_truth)
  samp   <- pps_sample(pop, TARGET_N)

  fit <- fit_wdpmm(samp$Y, samp$w, x2 = samp$x2, include_x2 = include_x2_in_fit)

  occ <- sort(table(fit$z_draws[nrow(fit$z_draws), ]), decreasing = TRUE)
  K_active <- min(K0 + 2, length(occ))  # allow slight over-clustering room

  rl <- relabel_and_summarize(fit, samp$w, K_active)
  sw <- sandwich_correct_prevalence(rl$prev_draws, samp$w, samp$psu_id,
                                     samp$stratum, rl$z_point_estimate, K_active)

  list(pi_hat = sw$pi_hat, adjusted_draws = sw$adjusted_draws,
       min_eig_trace = fit$min_eig_trace, K_active = K_active)
}

# -----------------------------------------------------------------------
# 6. PER-REPLICATE EXECUTION (SLURM-array-safe) + REPRODUCIBILITY HOOK
# -----------------------------------------------------------------------
#
# Designed to run as a SLURM job array: one task = one Monte Carlo
# replicate. The seed is a deterministic function of (arm, task_id), so
# ANY replicate -- including a failed one -- can be reproduced exactly
# later via reproduce_replicate(), regardless of cluster scheduling
# order or which other tasks ran. Nothing is saved except this seed
# mapping plus each replicate's small result/error object; the full
# population and sample are regenerated on demand rather than stored.
#
# Submit two array jobs, one per arm:
#   sbatch --array=1-500 --export=ALL,DPMM_ARM=informed run_dpmm_array.sbatch
#   sbatch --array=1-500 --export=ALL,DPMM_ARM=naive    run_dpmm_array.sbatch

RESULTS_DIR <- "results"
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

#' Map (arm, task_id) to a reproducible seed. Distinct offsets per arm
#' keep the two arms' RNG streams from ever coinciding for the same
#' task_id.
seed_for <- function(arm, task_id) {
  offset <- switch(arm, informed = 0L, naive = 200000L, stop("unknown arm: ", arm))
  offset + as.integer(task_id)
}

#' Run exactly one replicate under (arm, task_id): set the seed, run,
#' and save the result -- or, on failure, save the error plus the exact
#' seed needed to reproduce it -- keyed by task_id so tasks are
#' independently inspectable.
run_and_save_replicate <- function(arm, task_id) {
  seed <- seed_for(arm, task_id)
  set.seed(seed)

  out_path <- file.path(RESULTS_DIR, sprintf("result_%s_%05d.rds", arm, task_id))
  err_path <- file.path(RESULTS_DIR, sprintf("ERROR_%s_%05d.rds", arm, task_id))

  result <- tryCatch(
    run_one_replicate(include_x2_in_fit = (arm == "informed")),
    error = function(e) {
      saveRDS(list(arm = arm, task_id = task_id, seed = seed,
                   error = conditionMessage(e)),
              err_path)
      NULL
    }
  )

  if (is.null(result)) {
    cat(sprintf("[FAILED] arm=%s task_id=%d seed=%d -- see %s\n",
                arm, task_id, seed, err_path))
    cat(sprintf("Reproduce with: reproduce_replicate(\"%s\", %d)\n", arm, task_id))
    quit(status = 1, save = "no")  # nonzero exit -> SLURM marks task FAILED
  }

  result$arm <- arm; result$task_id <- task_id; result$seed <- seed
  saveRDS(result, out_path)
  cat(sprintf("[OK] arm=%s task_id=%d seed=%d -> %s\n", arm, task_id, seed, out_path))
}

#' Reproduce a specific (arm, task_id) replicate exactly -- e.g. to
#' debug a task that failed or looked odd on the cluster. Source this
#' script interactively (SLURM_ARRAY_TASK_ID unset, so nothing auto-runs
#' beyond the local demo guard below), then call:
#'   rep <- reproduce_replicate("naive", 340)
#' `rep` is bitwise-identical to what task_id=340 produced on SLURM.
reproduce_replicate <- function(arm, task_id) {
  set.seed(seed_for(arm, task_id))
  run_one_replicate(include_x2_in_fit = (arm == "informed"))
}

#' Aggregate all completed replicate files for one arm into the
#' coverage / K*-recovery / min-eigenvalue summary, as a post-hoc step
#' over saved per-task .rds files (run after the SLURM array completes).
aggregate_results <- function(arm) {
  files <- list.files(RESULTS_DIR, pattern = sprintf("^result_%s_.*\\.rds$", arm),
                       full.names = TRUE)
  if (length(files) == 0) stop("No completed result files found for arm: ", arm)

  true_prev <- get_true_prevalence()
  true_reduced <- true_prev[-length(true_prev)]

  n_files <- length(files)
  covered <- matrix(NA, n_files, length(true_reduced))
  K_active_vec <- integer(n_files)
  min_eig_vec <- numeric(n_files)

  for (i in seq_along(files)) {
    r <- readRDS(files[i])
    ci_lo <- apply(r$adjusted_draws, 2, quantile, probs = 0.025)
    ci_hi <- apply(r$adjusted_draws, 2, quantile, probs = 0.975)
    true_k <- true_reduced[seq_len(ncol(r$adjusted_draws))]
    covered[i, seq_along(true_k)] <- (true_k >= ci_lo[seq_along(true_k)]) &
                                       (true_k <= ci_hi[seq_along(true_k)])
    K_active_vec[i] <- r$K_active
    min_eig_vec[i] <- min(r$min_eig_trace, na.rm = TRUE)
  }

  list(n_completed = n_files,
       empirical_coverage = colMeans(covered, na.rm = TRUE),
       mean_K_active = mean(K_active_vec),
       min_eig_summary = summary(min_eig_vec))
}

# -----------------------------------------------------------------------
# 7. ENTRY POINT: one SLURM array task, or a small local interactive demo
# -----------------------------------------------------------------------

slurm_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
slurm_arm     <- Sys.getenv("DPMM_ARM", unset = "informed")

if (!is.na(slurm_task_id)) {
  # One SLURM array task = exactly one replicate: run it, save it, exit.
  run_and_save_replicate(slurm_arm, as.integer(slurm_task_id))
} else {
  # Local/interactive fallback -- small sequential demo, same functions
  # the cluster uses, so behavior matches exactly (just serialized).
  cat("No SLURM_ARRAY_TASK_ID set -- running a small local demo (R_REPLICATES\n")
  cat("per arm) instead of a single array task. See run_dpmm_array.sbatch\n")
  cat("for cluster submission.\n\n")

  for (r in seq_len(R_REPLICATES)) run_and_save_replicate("informed", r)
  for (r in seq_len(R_REPLICATES)) run_and_save_replicate("naive", r)

  cat("\n=== Informed arm (x2 included in fitted model) ===\n")
  print(aggregate_results("informed"))
  cat("\n=== Naive arm (x2 excluded, PPS1 informativeness case) ===\n")
  print(aggregate_results("naive"))
}
