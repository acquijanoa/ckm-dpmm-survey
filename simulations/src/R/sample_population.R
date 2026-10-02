# =============================================================================
# Draw informative two-stage PPS samples from a population (Kim, Rao & Wang
# 2024, Suppl. S9). Stage 1: N_PSU clusters per stratum, PPS to size M_hi.
# Stage 2: elements PPS to x2, for an overall sampling fraction f = 0.01.
# Weights 1/(pi_1 pi_2), normalized to sum to n.
#
#   Rscript sample_population.R <from> <to>   (or one SLURM array task each)
# Output: data/samples_<SCENARIO>/sample_NNNN.rds.
# =============================================================================

personal_lib <- "/nas/longleaf/home/aquijano/R/x86_64-pc-linux-gnu-library/4.5"
if (dir.exists(personal_lib) && !(personal_lib %in% .libPaths())) {
  .libPaths(c(personal_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(MASS)
})

CLUSTER_C_STR <- Sys.getenv("CLUSTER_C", unset = "0.3")

N_PSU_STR <- Sys.getenv("N_PSU", unset = "4")

# Scenario tag c<C>[_<scenario>][_npsu<N>]; defaults add no suffix.
POP_SCENARIO <- Sys.getenv("POP_SCENARIO", unset = "base")
pop_tag <- paste0("c", CLUSTER_C_STR, if (POP_SCENARIO != "base") paste0("_", POP_SCENARIO))
arm_tag <- pop_tag
if (N_PSU_STR != "4") arm_tag <- paste0(arm_tag, "_npsu", N_PSU_STR)

POP_FILE    <- sprintf("data/pop_%s/population.rds", pop_tag)
SAMPLES_DIR <- sprintf("data/samples_%s", arm_tag)
if (!dir.exists(SAMPLES_DIR)) dir.create(SAMPLES_DIR, recursive = TRUE)

SAMPLING_FRACTION <- 0.01
N_CLUSTERS_SAMPLE <- as.integer(N_PSU_STR)  # PSUs sampled per stratum
BASE_SEED         <- 42000L

draw_one_sample <- function(pop, sample_id) {
  set.seed(BASE_SEED + sample_id)

  H <- pop$H_STRATA
  sampled_indices <- integer(0)
  pi_1_vec <- numeric(0)
  pi_2_vec <- numeric(0)

  for (h in seq_len(H)) {
    h_clusters <- which(sapply(pop$cluster_info, function(cl) cl$stratum == h))
    M_h_clusters <- sapply(h_clusters, function(cid) pop$cluster_info[[cid]]$M_hi)
    M_h_total <- sum(M_h_clusters)

    # Stage 1: PPS to cluster size, pi_1 = n_h M_hi / M_h
    n_h <- min(N_CLUSTERS_SAMPLE, length(h_clusters))
    pi_1_h <- pmin(1.0, n_h * M_h_clusters / M_h_total)

    sel_cluster_idx <- sample(seq_along(h_clusters), size = n_h, replace = FALSE, prob = M_h_clusters)
    selected_cids <- h_clusters[sel_cluster_idx]
    selected_pi1  <- pi_1_h[sel_cluster_idx]

    # Stage 2: PPS to x2 within each cluster, m elements per cluster so the
    # overall sampling fraction is f; pi_2 = m x2_j / sum(x2)
    m_target <- max(2, round(SAMPLING_FRACTION * M_h_total / n_h))

    for (k in seq_along(selected_cids)) {
      cid <- selected_cids[k]
      pi1 <- selected_pi1[k]

      unit_indices <- which(pop$psu_id == cid)
      M_curr <- length(unit_indices)
      m_k <- min(m_target, M_curr)

      x2_vals <- pop$x2[unit_indices]
      prob_x2 <- x2_vals / sum(x2_vals)

      pi_2_units <- pmin(0.999, m_k * prob_x2)

      sel_units <- sample(unit_indices, size = m_k, replace = FALSE, prob = prob_x2)

      sampled_indices <- c(sampled_indices, sel_units)
      pi_1_vec <- c(pi_1_vec, rep(pi1, m_k))
      pi_2_vec <- c(pi_2_vec, pi_2_units[match(sel_units, unit_indices)])
    }
  }

  pi_joint <- pi_1_vec * pi_2_vec
  w_raw <- 1 / pi_joint

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

# Sample range: command-line <from> <to>, else the SLURM array task, else 1-1000
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
