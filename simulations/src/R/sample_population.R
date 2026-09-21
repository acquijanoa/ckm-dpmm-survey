# =============================================================================
# Two-Stage PPS Sampling for DPMM Survey Phenotyping Study
# Based on Kim, Rao, & Wang (2024, JASA Supplementary Section S9)
# =============================================================================

personal_lib <- "/nas/longleaf/home/aquijano/R/x86_64-pc-linux-gnu-library/4.5"
if (dir.exists(personal_lib) && !(personal_lib %in% .libPaths())) {
  .libPaths(c(personal_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(MASS)
})

# Same CLUSTER_C grid parameter as generate_population.R. One folder per
# grid point for pop/samples (bulk du/tar/rm per arm); each saved object
# also carries a cluster_c field so the value is traceable from the file
# alone.
CLUSTER_C_STR <- Sys.getenv("CLUSTER_C", unset = "0.3")

# N_PSU grid parameter: number of PSUs (clusters) sampled per stratum, at
# fixed overall sampling fraction f -- more PSUs means each one contributes
# fewer elements (m_target below), spreading the same total sample across
# more, smaller clusters. Default 4 matches the original (pre-grid) design
# and its folder/filenames, so the already-computed N_PSU=4 baseline is
# unaffected; other values get an explicit _npsuN suffix.
N_PSU_STR <- Sys.getenv("N_PSU", unset = "4")

# arm_tag: c<C>[_npsu<N>] -- N_PSU only appends a suffix when it differs
# from its original (pre-grid) default, so already-computed baselines keep
# their original folder/file names.
arm_tag <- sprintf("c%s", CLUSTER_C_STR)
if (N_PSU_STR != "4") arm_tag <- paste0(arm_tag, "_npsu", N_PSU_STR)

POP_FILE    <- sprintf("data/pop_c%s/population.rds", CLUSTER_C_STR)
SAMPLES_DIR <- sprintf("data/samples_%s", arm_tag)
if (!dir.exists(SAMPLES_DIR)) dir.create(SAMPLES_DIR, recursive = TRUE)

SAMPLING_FRACTION <- 0.01  # f = 0.01 (WR ~ WOR)
N_CLUSTERS_SAMPLE <- as.integer(N_PSU_STR)  # n_h PSUs sampled per stratum
BASE_SEED         <- 42000L

draw_one_sample <- function(pop, sample_id) {
  set.seed(BASE_SEED + sample_id)

  H <- pop$H_STRATA
  sampled_indices <- integer(0)
  pi_1_vec <- numeric(0)
  pi_2_vec <- numeric(0)

  # Pre-split population indices by cluster
  # pop$psu_id is 1-based cluster_id
  for (h in seq_len(H)) {
    # Clusters in stratum h
    h_clusters <- which(sapply(pop$cluster_info, function(cl) cl$stratum == h))
    M_h_clusters <- sapply(h_clusters, function(cid) pop$cluster_info[[cid]]$M_hi)
    M_h_total <- sum(M_h_clusters)

    # Stage 1: Sample n_h clusters with PPS proportional to cluster size M_hi
    n_h <- min(N_CLUSTERS_SAMPLE, length(h_clusters))
    # Inclusion prob at Stage 1: pi_{1, hi} = n_h * M_{hi} / M_h
    pi_1_h <- pmin(1.0, n_h * M_h_clusters / M_h_total)

    # PPS selection without replacement (WR ~ WOR under small f)
    sel_cluster_idx <- sample(seq_along(h_clusters), size = n_h, replace = FALSE, prob = M_h_clusters)
    selected_cids <- h_clusters[sel_cluster_idx]
    selected_pi1  <- pi_1_h[sel_cluster_idx]

    # Stage 2: Within each selected cluster, sample elements with PPS on x2
    # Target elements per cluster to achieve overall sampling fraction f
    # m_hi = round(f * M_h_total / n_h)
    m_target <- max(2, round(SAMPLING_FRACTION * M_h_total / n_h))

    for (k in seq_along(selected_cids)) {
      cid <- selected_cids[k]
      pi1 <- selected_pi1[k]

      unit_indices <- which(pop$psu_id == cid)
      M_curr <- length(unit_indices)
      m_k <- min(m_target, M_curr)

      x2_vals <- pop$x2[unit_indices]
      prob_x2 <- x2_vals / sum(x2_vals)

      # Stage 2 inclusion prob: pi_{2, j | hi} = m_k * x2_j / sum(x2)
      pi_2_units <- pmin(0.999, m_k * prob_x2)

      # Sample elements (WOR ~ WR)
      sel_units <- sample(unit_indices, size = m_k, replace = FALSE, prob = prob_x2)

      sampled_indices <- c(sampled_indices, sel_units)
      pi_1_vec <- c(pi_1_vec, rep(pi1, m_k))
      pi_2_vec <- c(pi_2_vec, pi_2_units[match(sel_units, unit_indices)])
    }
  }

  # Joint inclusion probabilities: pi_{hij} = pi_{1, hi} * pi_{2, j | hi}
  pi_joint <- pi_1_vec * pi_2_vec
  w_raw <- 1 / pi_joint

  # Standard normalization: scale weights to sum to realized sample size n
  n_realized <- length(sampled_indices)
  w_norm <- w_raw / sum(w_raw) * n_realized

  list(
    sample_id = sample_id,
    cluster_c = pop$cluster_c,
    n_psu = N_CLUSTERS_SAMPLE,
    sampling_frac = SAMPLING_FRACTION,
    n = n_realized,
    idx = sampled_indices,
    Y = pop$Y[sampled_indices, , drop = FALSE],
    x1 = pop$x1[sampled_indices],
    x2 = pop$x2[sampled_indices],
    z_true = pop$z_true[sampled_indices],
    stratum = pop$stratum[sampled_indices],
    psu_id = pop$psu_id[sampled_indices],
    pi_i = pi_joint,
    w = w_norm
  )
}

# Entry point: process specified samples
args <- commandArgs(trailingOnly = TRUE)
if (!file.exists(POP_FILE)) stop("Population file not found: ", POP_FILE)

pop <- readRDS(POP_FILE)

if (length(args) >= 2) {
  start_id <- as.integer(args[1])
  end_id   <- as.integer(args[2])
} else if (length(args) == 1) {
  start_id <- as.integer(args[1])
  end_id   <- as.integer(args[1])
} else {
  slurm_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
  if (!is.na(slurm_id)) {
    start_id <- as.integer(slurm_id)
    end_id   <- as.integer(slurm_id)
  } else {
    # Default: run all 1000 samples
    start_id <- 1L
    end_id   <- 1000L
  }
}

cat(sprintf("Drawing samples %d to %d from population (N = %d)...\n",
            start_id, end_id, pop$N_total))

t0 <- proc.time()
for (s_id in start_id:end_id) {
  s_obj <- draw_one_sample(pop, s_id)
  out_path <- file.path(SAMPLES_DIR, sprintf("sample_%04d.rds", s_id))
  saveRDS(s_obj, out_path)
  if (s_id %% 100 == 0 || s_id == end_id) {
    cat(sprintf("  [OK] Sample %04d saved: n = %d units, %d PSUs, %d strata\n",
                s_id, s_obj$n, length(unique(s_obj$psu_id)), length(unique(s_obj$stratum))))
  }
}

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("Completed %d sample(s) in %.2f seconds (avg %.3f s/sample).\n",
            end_id - start_id + 1, elapsed, elapsed / (end_id - start_id + 1)))
