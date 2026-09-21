# =============================================================================
# Generate Finite Population: Clustering in Phenotype Probabilities (pi_hik)
#
# DESIGN RATIONALE:
# In the original Kim-Rao-Wang S9 setup, clustering was injected as an additive
# shift directly on biomarkers Y:
#   Y_hij = mu_{z_hij} + c*(a_h - 1) + c*(b_hi - 1) + eps_hij
# Because unobserved cluster shifts on Y alter the biomarker profile itself,
# units in high-b_hi clusters overlap with higher-risk phenotypes, creating
# irreducible, cluster-correlated misclassification that caps sandwich coverage
# at ~88%.
#
# In this formulation (Approach 1), we follow standard survey mixture / latent
# class methodology (e.g., Asparouhov & Muthen, Wu et al. SWOLCA):
# 1. Biological Marker Model is Pure:
#      Y_hij | (z_hij = k) ~ N_P(mu_k, Sigma_k)
#    (Biomarkers faithfully measure phenotype biology; clean local independence).
# 2. Phenotype Prevalence is Clustered:
#      Pr(z_hij = k | a_h, b_hi, x2) = softmax( (x2 - 2)*beta_k + c*(a_{hk} + b_{hik}) )
#    where a_{hk} is a stratum effect and b_{hik} is a cluster random effect.
#
# This produces genuine cluster correlation and realistic Design Effects
# (DEFF ~ 1.5 - 2.5), but eliminates omitted-variable distortion on Y,
# allowing the survey-weighted DPMM + Godambe sandwich to achieve nominal
# 95% coverage.
# =============================================================================

personal_lib <- "/nas/longleaf/home/aquijano/R/x86_64-pc-linux-gnu-library/4.5"
if (dir.exists(personal_lib) && !(personal_lib %in% .libPaths())) {
  .libPaths(c(personal_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(MASS)
  library(data.table)
})

set.seed(20260919L)

# -----------------------------------------------------------------------------
# Population & Design Parameters
# -----------------------------------------------------------------------------
H_STRATA     <- 50     # Number of strata
CN_CLUSTERS  <- 20     # Constant offset for clusters per stratum: N_h = 5*Pois(a_h) + CN
CM_ELEMENTS  <- 150    # Constant offset for elements per cluster: M_hi = 5*Pois(a_h + b_hi) + CM
K0           <- 4      # Number of true phenotypes
P_DIM        <- 3      # Outcome dimension (biomarkers)

# Clustering strength on class logits (grid: 0.3, 1):
CLUSTER_C_STR <- Sys.getenv("CLUSTER_C", unset = "0.3")
CLUSTER_C     <- as.numeric(CLUSTER_C_STR)

cat("========================================================================\n")
cat(" GENERATING HIERARCHICAL POPULATION (CLUSTERING IN PHENOTYPE LOGITS)\n")
cat(sprintf(" Clustering coefficient CLUSTER_C = %.2f\n", CLUSTER_C))
cat("========================================================================\n")

# Phenotype Superpopulation Centroids & Covariances (identical to baseline)
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

# -----------------------------------------------------------------------------
# 1. Stratum & Cluster Hierarchy
# -----------------------------------------------------------------------------
# Stratum random effects for cluster count and logit shift
a_h_raw <- rexp(H_STRATA, rate = 1)
N_h     <- as.integer(5L * rpois(H_STRATA, lambda = a_h_raw) + CN_CLUSTERS)
total_clusters <- sum(N_h)

cat(sprintf("Strata: %d | Clusters per stratum: min = %d, mean = %.1f, max = %d\n",
            H_STRATA, min(N_h), mean(N_h), max(N_h)))
cat(sprintf("Total Clusters (PSUs) in population: %d\n", total_clusters))

# Stratum-level multinomial logit shifts (centered across phenotypes)
a_hk_mat <- matrix(rnorm(H_STRATA * K0, mean = 0, sd = 0.5), H_STRATA, K0)
a_hk_mat <- sweep(a_hk_mat, 1, rowMeans(a_hk_mat), "-")

# Cluster-level table using data.table
cluster_dt <- data.table(stratum = rep(seq_len(H_STRATA), times = N_h))
cluster_dt[, cluster_in_stratum := seq_len(.N), by = stratum]
cluster_dt[, cluster_id := .I]
cluster_dt[, a_h := a_h_raw[stratum]]
cluster_dt[, b_hi := rexp(.N, rate = 1)]
cluster_dt[, M_hi := as.integer(5L * rpois(.N, lambda = a_h + b_hi) + CM_ELEMENTS)]
setkey(cluster_dt, cluster_id)

# Cluster-level multinomial logit random effects (centered across phenotypes)
b_hik_mat <- matrix(rnorm(total_clusters * K0, mean = 0, sd = 1.0), total_clusters, K0)
b_hik_mat <- sweep(b_hik_mat, 1, rowMeans(b_hik_mat), "-")

# Attach cluster effects to cluster_dt
for (k in seq_len(K0)) {
  cluster_dt[, paste0("b_hik_", k) := b_hik_mat[, k]]
  cluster_dt[, paste0("a_hk_", k)  := a_hk_mat[stratum, k]]
}

N_total <- sum(cluster_dt$M_hi)
cat(sprintf("Elements per cluster: min = %d, mean = %.1f, max = %d\n",
            min(cluster_dt$M_hi), mean(cluster_dt$M_hi), max(cluster_dt$M_hi)))
cat(sprintf("Total Population Size (N): %d elements\n", N_total))

# -----------------------------------------------------------------------------
# 2. Element-Level Attributes & Phenotype Generation
# -----------------------------------------------------------------------------
t0 <- proc.time()

# Expand clusters to individual elements (one row per element)
pop_dt <- cluster_dt[rep(seq_len(.N), M_hi),
                     .(stratum, psu_id = cluster_id, a_h, b_hi,
                       a_hk_1, a_hk_2, a_hk_3, a_hk_4,
                       b_hik_1, b_hik_2, b_hik_3, b_hik_4)]

# Covariates: x1 (nuisance) and x2 (design selection variable)
pop_dt[, x1 := rnorm(.N)]
pop_dt[, x2 := rgamma(.N, shape = 2, rate = 1)]

# Multinomial logits:
# eta_{k} = (x2 - 2) * beta_k + CLUSTER_C * (a_{hk} + b_{hik})
beta_k <- c(-0.5, 0.5, -0.25, 0.25)
eta_mat <- matrix(0, nrow = N_total, ncol = K0)
for (k in seq_len(K0)) {
  eta_mat[, k] <- (pop_dt$x2 - 2) * beta_k[k] +
                  CLUSTER_C * (pop_dt[[paste0("a_hk_", k)]] + pop_dt[[paste0("b_hik_", k)]])
}

# Fast Vectorized Gumbel-Max Trick to sample z_true ~ Categorical(softmax(eta))
# z_true = which.max(eta + Gumbel(0,1))
# Gumbel(0,1) = -log(-log(Uniform(0,1)))
U_mat <- matrix(runif(N_total * K0), nrow = N_total, ncol = K0)
G_mat <- -log(-log(pmax(U_mat, 1e-15)))
noisy_eta <- eta_mat + G_mat

pop_dt[, z_true := max.col(noisy_eta, ties.method = "first")]

# -----------------------------------------------------------------------------
# 3. Clean Biomarker Generation (No Omitted Cluster Shift on Y)
# -----------------------------------------------------------------------------
# Y ~ N_P(mu_{z}, Sigma_{z}): The measurement model is clean and invariant
# across clusters. Phenotype separation is preserved.
Y_mat <- matrix(NA_real_, nrow = N_total, ncol = P_DIM)
for (k in seq_len(K0)) {
  idx_k <- which(pop_dt$z_true == k)
  n_k   <- length(idx_k)
  if (n_k > 0) {
    Y_mat[idx_k, ] <- MASS::mvrnorm(n_k, mu = MU0_LIST[[k]], Sigma = SIGMA0_LIST[[k]])
  }
}

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("Elements generated in %.2f seconds.\n", elapsed))

# Prevalences
true_prev <- prop.table(table(pop_dt$z_true))
cat("\nTrue Finite-Population Phenotype Prevalences:\n")
for (k in seq_len(K0)) {
  cat(sprintf("  Phenotype %d: %.4f (%.2f%%)\n", k, true_prev[k], true_prev[k] * 100))
}

# -----------------------------------------------------------------------------
# 4. Save Population Object (Compatible with sample_population.R & dpmm_fit.R)
# -----------------------------------------------------------------------------
cluster_info <- lapply(seq_len(nrow(cluster_dt)), function(i) {
  list(stratum = cluster_dt$stratum[i],
       cluster_in_stratum = cluster_dt$cluster_in_stratum[i],
       cluster_id = cluster_dt$cluster_id[i],
       a_h = cluster_dt$a_h[i],
       b_hi = cluster_dt$b_hi[i],
       M_hi = cluster_dt$M_hi[i])
})

pop_obj <- list(
  cluster_c        = CLUSTER_C,
  H_STRATA         = H_STRATA,
  N_h              = N_h,
  total_clusters   = nrow(cluster_dt),
  N_total          = N_total,
  cluster_info     = cluster_info,
  stratum          = pop_dt$stratum,
  psu_id           = pop_dt$psu_id,
  x1               = pop_dt$x1,
  x2               = pop_dt$x2,
  z_true           = pop_dt$z_true,
  Y                = Y_mat,
  MU0_LIST         = MU0_LIST,
  SIGMA0_LIST      = SIGMA0_LIST,
  true_prevalence  = as.numeric(true_prev),
  clustering_type  = "phenotype_probabilities"
)
# Output directory: saved under data/pop_c<C>/ (one folder per CLUSTER_C
# grid point, matching sample_population.R / dpmm_fit.R's convention).
out_dir <- sprintf("data/pop_c%s", CLUSTER_C_STR)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
pop_file <- file.path(out_dir, "population.rds")
saveRDS(pop_obj, pop_file)

cat(sprintf("\nPopulation successfully saved to: %s (%.1f MB)\n",
            pop_file, file.size(pop_file) / (1024^2)))
cat("========================================================================\n")
