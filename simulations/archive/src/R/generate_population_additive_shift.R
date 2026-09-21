# =============================================================================
# Generate Finite Population for DPMM Survey Phenotyping Study
# Based on Kim, Rao, & Wang (2024, JASA Supplementary Section S9)
# =============================================================================

suppressPackageStartupMessages({
  library(MASS)
  library(data.table)
})

# POP_ID selects a superpopulation replicate for the multi-population
# variance-decomposition study (25 populations x 100 samples each, to
# separate superpopulation-regeneration variance from survey-sampling
# variance under a FIXED population). Unset (default) reproduces the
# original single, fixed population exactly -- same seed, same output path
# -- so the already-running single-population jobs are unaffected.
POP_ID_STR <- Sys.getenv("POP_ID", unset = "")
IS_MULTIPOP <- nzchar(POP_ID_STR)
POP_ID <- if (IS_MULTIPOP) as.integer(POP_ID_STR) else NA_integer_

set.seed(if (IS_MULTIPOP) 20260919L + 137L * POP_ID else 20260919L)

# Parameters
H_STRATA     <- 50     # Number of strata
CN_CLUSTERS  <- 20     # Constant offset for clusters per stratum: N_h = 5*Pois(a_h) + CN
CM_ELEMENTS  <- 150    # Constant offset for elements per cluster: M_hi = 5*Pois(a_h + b_hi) + CM
K0           <- 4      # Number of true phenotypes
P_DIM        <- 3      # Outcome dimension

# Clustering-strength grid parameter: a single shared coefficient applied to
# BOTH the stratum-level effect a_h and the cluster-level effect b_hi in the
# outcome shift below (previously 0.2 and 0.3 respectively -- now equal, so
# one knob controls the overall degree of clustering/design effect in the
# population). Passed via env var so the same script can be run once per
# point on the {0.3, 1} grid without code duplication.
CLUSTER_C_STR <- Sys.getenv("CLUSTER_C", unset = "0.3")
CLUSTER_C     <- as.numeric(CLUSTER_C_STR)
cat(sprintf("Clustering-strength coefficient CLUSTER_C = %s\n", CLUSTER_C_STR))

# Phenotype Superpopulation Centroids & Covariances
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

cat("========================================================================\n")
cat(" GENERATING HIERARCHICAL POPULATION (Kim, Rao, & Wang 2024 JASA S9)\n")
cat("========================================================================\n")

# 1. Stratum-level effects and cluster counts
a_h <- rexp(H_STRATA, rate = 1)
N_h <- 5 * rpois(H_STRATA, lambda = a_h) + CN_CLUSTERS
total_clusters <- sum(N_h)

cat(sprintf("Strata: %d | Clusters per stratum: min = %d, mean = %.1f, max = %d\n",
            H_STRATA, min(N_h), mean(N_h), max(N_h)))
cat(sprintf("Total Clusters (PSUs) in population: %d\n", total_clusters))

# 2. Cluster-level effects and cluster sizes (data.table, vectorized)
cluster_dt <- data.table(stratum = rep(seq_len(H_STRATA), times = N_h))
cluster_dt[, cluster_in_stratum := seq_len(.N), by = stratum]
cluster_dt[, cluster_id := .I]
cluster_dt[, a_h := a_h[stratum]]
cluster_dt[, b_hi := rexp(.N, rate = 1)]
cluster_dt[, M_hi := as.integer(5L * rpois(.N, lambda = a_h + b_hi) + CM_ELEMENTS)]
setkey(cluster_dt, cluster_id)

N_total <- sum(cluster_dt$M_hi)

cat(sprintf("Elements per cluster: min = %d, mean = %.1f, max = %d\n",
            min(cluster_dt$M_hi), mean(cluster_dt$M_hi), max(cluster_dt$M_hi)))
cat(sprintf("Total Population Size (N): %d elements\n", N_total))

# 3. Element-level attributes and outcome generation (data.table, vectorized)
t0 <- proc.time()

# Expand cluster-level rows to element level by M_hi (one row per element),
# carrying stratum/cluster_id/a_h/b_hi down to every element in that cluster.
pop_dt <- cluster_dt[rep(seq_len(.N), M_hi),
                      .(stratum, psu_id = cluster_id, a_h, b_hi)]

# Covariates (vectorized over the full population at once)
pop_dt[, x1 := rnorm(.N)]
pop_dt[, x2 := rgamma(.N, shape = 2, rate = 1)]

# Latent phenotype membership: multinomial logit on x2
beta_k <- c(-0.5, 0.5, -0.25, 0.25)
eta <- outer(pop_dt$x2 - 2, beta_k)
probs <- exp(eta) / rowSums(exp(eta))
pop_dt[, z_true := apply(probs, 1, function(p) sample.int(K0, 1, prob = p))]

# Centered cluster shift (PPS3 mechanism); CLUSTER_C weights the stratum-
# level (a_h) and cluster-level (b_hi) effects EQUALLY, unlike the original
# fixed 0.2/0.3 split.
pop_dt[, cluster_shift := CLUSTER_C * (a_h - 1) + CLUSTER_C * (b_hi - 1)]

# Outcome Y = g(X'\beta + a_h + b_{hi}), drawn per true-phenotype group then
# shifted by each element's own cluster_shift.
Y_mat <- matrix(NA_real_, N_total, P_DIM)
for (k in seq_len(K0)) {
  idx_k <- which(pop_dt$z_true == k)
  if (length(idx_k) > 0) {
    Y_mat[idx_k, ] <- MASS::mvrnorm(length(idx_k), mu = MU0_LIST[[k]], Sigma = SIGMA0_LIST[[k]])
  }
}
Y_mat <- Y_mat + pop_dt$cluster_shift  # recycles column-wise, i.e. per-row

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("Elements generated in %.2f seconds.\n", elapsed))

# True population prevalences
true_prev <- prop.table(table(pop_dt$z_true))
cat("\nTrue Population Phenotype Prevalences:\n")
for (k in seq_len(K0)) {
  cat(sprintf("  Phenotype %d: %.4f (%.2f%%)\n", k, true_prev[k], true_prev[k] * 100))
}

# 4. Save Population Object -- cluster_info is reconstructed as a list of
# lists (rather than kept as cluster_dt) so sample_population.R and the
# fitting scripts, which index it as pop$cluster_info[[cid]]$field, are
# unaffected by this internal data.table refactor.
cluster_info <- lapply(seq_len(nrow(cluster_dt)), function(i) {
  list(stratum = cluster_dt$stratum[i],
       cluster_in_stratum = cluster_dt$cluster_in_stratum[i],
       cluster_id = cluster_dt$cluster_id[i],
       a_h = cluster_dt$a_h[i],
       b_hi = cluster_dt$b_hi[i],
       M_hi = cluster_dt$M_hi[i])
})

pop_obj <- list(
  cluster_c = CLUSTER_C,
  H_STRATA = H_STRATA,
  N_h = N_h,
  total_clusters = nrow(cluster_dt),
  N_total = N_total,
  cluster_info = cluster_info,
  stratum = pop_dt$stratum,
  psu_id = pop_dt$psu_id,
  x1 = pop_dt$x1,
  x2 = pop_dt$x2,
  z_true = pop_dt$z_true,
  Y = Y_mat,
  MU0_LIST = MU0_LIST,
  SIGMA0_LIST = SIGMA0_LIST,
  true_prevalence = as.numeric(true_prev)
)

pop_obj$pop_id <- POP_ID

# One folder per CLUSTER_C grid point (bulk du/tar/rm per arm); pop_obj also
# carries cluster_c internally so the value is traceable from the file
# alone. Multi-population replicates nest under a separate
# multipop_c<C>/pop_<PP>/ tree so they never collide with the single fixed
# population's data/pop_c<C>/population.rds.
out_dir <- if (IS_MULTIPOP) {
  file.path(sprintf("data/multipop_c%s", CLUSTER_C_STR), sprintf("pop_%02d", POP_ID))
} else {
  sprintf("data/pop_c%s", CLUSTER_C_STR)
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
pop_file <- file.path(out_dir, "population.rds")
saveRDS(pop_obj, pop_file)

cat(sprintf("\nPopulation successfully saved to: %s (%.1f MB)\n",
            pop_file, file.size(pop_file) / (1024^2)))
cat("========================================================================\n")
