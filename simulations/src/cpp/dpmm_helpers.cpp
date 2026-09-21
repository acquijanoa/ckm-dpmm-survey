#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//
// Sparse finite-mixture weights (Rousseau & Mengersen, 2011; Malsiner-Walli,
// Fruehwirth-Schnatter & Gruen, 2016 "sparse finite mixtures"; and the same
// principle used in practice by the WOLCA/SWOLCA survey-weighted overfitted
// latent class models -- Stephenson, Wu & Dominici 2024; Wu, Williams,
// Savitsky & Stephenson 2024 -- which report that empty components "drop
// out during Gibbs sampling" under exactly this construction).
//
// Motivation: the quasi-Bernoulli stick-breaking variant (dpmm_helpers_qb.cpp)
// was pilot-tested and showed *zero* reduction in raw active-component count
// (K_raw_active stayed at 10.0/10 across 80 replicates, two priors) --
// because it only prunes components with near-zero occupancy, and this
// model's spurious extra components have genuinely nontrivial occupancy
// (real over-splitting from the likelihood), not idle noise. This sparse
// Dirichlet approach targets that directly: instead of asking "is this one
// component empty," it puts a joint sparsity-inducing prior on the *entire*
// L-vector of weights, Dirichlet(alpha0/L, ..., alpha0/L) with alpha0 small,
// which is known (Rousseau & Mengersen 2011) to drive the posterior weight
// of superfluous components to (asymptotically) zero even when none of them
// is individually empty at any given iteration.
//
// Only the naive (covariate-free) arm's helpers live here -- the x2-informed
// stick-breaking mechanism and the bootstrap-based sandwich it was paired
// with belong to the archived informed/naive-arm comparison and have moved
// to archive/src/cpp/dpmm_helpers_informed_arm.cpp, since neither is called
// by the active src/R/dpmm_fit.R pipeline.
//
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// ---------------------------------------------------------------------------
// 1. Vectorized sampling of cluster indicators z_i from probabilities
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerVector sample_z_cpp(const NumericMatrix& probs) {
    int n = probs.nrow();
    int L = probs.ncol();
    IntegerVector z(n);

    for (int i = 0; i < n; ++i) {
        double u = R::runif(0.0, 1.0);
        double cum = 0.0;
        int chosen = L;
        for (int k = 0; k < L; ++k) {
            cum += probs(i, k);
            if (u <= cum || k == L - 1) {
                chosen = k + 1; // 1-based index for R
                break;
            }
        }
        z[i] = chosen;
    }
    return z;
}

// ---------------------------------------------------------------------------
// 2. Sparse Dirichlet mixture weights (naive arm)
// ---------------------------------------------------------------------------
// pi ~ Dirichlet(alpha0/L + n_1, ..., alpha0/L + n_L), sampled via
// independent Gamma(shape, 1) draws normalized to sum to 1.
// [[Rcpp::export]]
NumericVector sparse_dirichlet_weights_cpp(const IntegerVector& z,
                                            const NumericVector& w,
                                            int L, double alpha0) {
    int n = z.size();
    std::vector<double> n_k(L + 1, 0.0); // 1-based
    for (int i = 0; i < n; ++i) n_k[z[i]] += w[i];

    NumericVector pi0(L);
    double total = 0.0;
    double shape0 = alpha0 / (double)L;
    for (int k = 1; k <= L; ++k) {
        double shape = shape0 + n_k[k];
        double g = R::rgamma(shape, 1.0);
        g = std::max(g, 1e-300);
        pi0[k - 1] = g;
        total += g;
    }
    for (int k = 0; k < L; ++k) pi0[k] /= total;
    return pi0;
}

// ---------------------------------------------------------------------------
// 3. Stratified Delete-One-PSU Jackknife Covariance (J_hat)
// ---------------------------------------------------------------------------
// For a LINEAR (additive-over-PSU) score total T = sum_h sum_i t_hi, the
// delete-one-PSU jackknife variance estimator reduces algebraically to the
// closed-form Rao-Wu "ultimate cluster" with-replacement variance:
//   J_hat = sum_h [n_h/(n_h-1)] * sum_{i in h} (t_hi - tbar_h)(t_hi - tbar_h)'
// This is exactly the B -> infinity limit of a stratified PSU bootstrap's
// Monte Carlo estimate of the same quantity (same expectation, zero
// resampling noise, and no B hyperparameter). Strata with a single sampled
// PSU (n_h < 2) cannot be jackknifed and contribute zero -- same degenerate
// behavior a bootstrap has there (always resampling the one available PSU).
// [[Rcpp::export]]
arma::mat jackknife_sandwich_cpp(const arma::mat& psu_scores,
                                  const Rcpp::List& stratum_psu_indices) {
    int Kd = psu_scores.n_cols;
    int H = stratum_psu_indices.size();
    arma::mat J_hat(Kd, Kd, arma::fill::zeros);

    for (int h = 0; h < H; ++h) {
        Rcpp::IntegerVector psus_in_h = stratum_psu_indices[h];
        int nh = psus_in_h.size();
        if (nh < 2) continue;

        arma::rowvec tbar_h(Kd, arma::fill::zeros);
        for (int i = 0; i < nh; ++i) tbar_h += psu_scores.row(psus_in_h[i]);
        tbar_h /= (double)nh;

        arma::mat stratum_contrib(Kd, Kd, arma::fill::zeros);
        for (int i = 0; i < nh; ++i) {
            arma::rowvec dev = psu_scores.row(psus_in_h[i]) - tbar_h;
            stratum_contrib += dev.t() * dev;
        }
        J_hat += (nh / (double)(nh - 1)) * stratum_contrib;
    }

    return J_hat;
}
