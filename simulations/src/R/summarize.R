# =============================================================================
# Summarize dpmm_fit.R replicates per arm into summaries/summary_<tag>.csv/.rds.
#
#   Rscript src/R/summarize.R                 # every arm in results/
#   Rscript src/R/summarize.R <tag> [<tag>]   # named arms
#
# Per phenotype: bias, empirical vs. estimated SE, coverage of the truth and
# of the across-replicate mean (centered coverage). Intervals are
# pi_hat +/- t_crit * se. Replicates whose SE was flagged as unusable
# (se_flag) are counted and left out of the SE and coverage columns.
# =============================================================================

get_project_root <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))))
  getwd()
}
PROJECT_ROOT  <- get_project_root()
RESULTS_DIR   <- file.path(PROJECT_ROOT, "results")
SUMMARIES_DIR <- file.path(PROJECT_ROOT, "summaries")

result_files <- function(tag) {
  pattern <- sprintf("^result_%s_[0-9]+\\.rds$", gsub("\\.", "\\\\.", tag))
  list.files(RESULTS_DIR, pattern = pattern, full.names = TRUE)
}

list_tags <- function() {
  f <- list.files(RESULTS_DIR, pattern = "^result_.*_[0-9]+\\.rds$")
  unique(sub("^result_(.*)_[0-9]+\\.rds$", "\\1", f))
}

summarize_arm <- function(tag) {
  files <- result_files(tag)
  if (length(files) == 0) stop("No result files for arm ", tag)
  rs <- lapply(files, readRDS)
  K <- length(rs[[1]]$pi_hat)

  pi_mat <- do.call(rbind, lapply(rs, `[[`, "pi_hat"))
  truth  <- colMeans(do.call(rbind, lapply(rs, `[[`, "true_prevalence")))
  se_mat <- do.call(rbind, lapply(rs, `[[`, "se"))
  half   <- se_mat * vapply(rs, `[[`, 0, "t_crit")
  lo     <- pi_mat - half
  hi     <- pi_mat + half
  pi_bar <- colMeans(pi_mat)

  covers <- function(target) {
    colMeans(sweep(lo, 2, target, "<=") & sweep(hi, 2, target, ">="), na.rm = TRUE)
  }
  emp_se <- apply(pi_mat, 2, sd)
  est_se <- colMeans(se_mat, na.rm = TRUE)

  tbl <- data.frame(
    phenotype           = seq_len(K),
    true_prevalence     = truth,
    estimate            = pi_bar,
    bias                = pi_bar - truth,
    rel_bias_pct        = 100 * (pi_bar - truth) / truth,
    empirical_se        = emp_se,
    estimated_se        = est_se,
    se_ratio            = emp_se / est_se,
    coverage_pct        = 100 * covers(truth),
    centered_coverage_pct = 100 * covers(pi_bar)
  )

  get_num <- function(name) vapply(rs, function(r) if (is.null(r[[name]])) NA_real_ else as.numeric(r[[name]]), 0)
  k_hat <- get_num("K_hat_estimated")
  worst_dev <- apply(abs(sweep(pi_mat, 2, truth)), 1, max)
  diagnostics <- list(
    arm_tag              = tag,
    n_completed          = length(rs),
    n_errors             = length(list.files(RESULTS_DIR, sprintf("^ERROR_%s_[0-9]+\\.rds$", gsub("\\.", "\\\\.", tag)))),
    k0_recovery_pct      = if (all(is.na(k_hat))) NA_real_ else 100 * mean(k_hat == K, na.rm = TRUE),
    mean_K_raw_active    = mean(get_num("K_raw_active")),
    mean_df_design       = mean(get_num("df_design")),
    n_fallback_reps      = sum(get_num("n_dens_fallback") + get_num("n_sigma_fallback") > 0, na.rm = TRUE),
    n_outliers_0.10      = sum(worst_dev > 0.10),
    n_se_flagged         = sum(!is.finite(rowSums(se_mat))),
    mean_coverage_pct    = mean(tbl$coverage_pct)
  )

  dir.create(SUMMARIES_DIR, showWarnings = FALSE, recursive = TRUE)
  write.csv(tbl, file.path(SUMMARIES_DIR, sprintf("summary_%s.csv", tag)), row.names = FALSE)
  saveRDS(c(diagnostics, list(table = tbl)), file.path(SUMMARIES_DIR, sprintf("summary_%s.rds", tag)))
  invisible(c(diagnostics, list(table = tbl)))
}

print_summary <- function(s) {
  cat(sprintf("\n== %s  (n=%d, errors=%d, K0 recovery=%s, mean coverage=%.1f%%) ==\n",
              s$arm_tag, s$n_completed, s$n_errors,
              if (is.na(s$k0_recovery_pct)) "n/a" else sprintf("%.1f%%", s$k0_recovery_pct),
              s$mean_coverage_pct))
  cat(sprintf("   outliers(>0.10)=%d, SE flagged=%d, reps with fallbacks=%d\n",
              s$n_outliers_0.10, s$n_se_flagged, s$n_fallback_reps))
  print(format(s$table, digits = 3), row.names = FALSE)
}

if (sys.nframe() == 0L) {
  tags <- commandArgs(trailingOnly = TRUE)
  if (length(tags) == 0) tags <- list_tags()
  for (tag in tags) print_summary(summarize_arm(tag))
}
