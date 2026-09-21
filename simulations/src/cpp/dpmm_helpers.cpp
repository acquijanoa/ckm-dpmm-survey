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
// Two components:
//   1. sparse_dirichlet_weights_cpp(): exact conjugate Gibbs update for the
//      covariate-free (naive-arm) mixture weights pi ~ Dirichlet(alpha0/L +
//      n_1, ..., alpha0/L + n_L), n_k = weighted count assigned to k.
//   2. compute_stick_breaking_sparse_cpp(): for the x2-informed arm, keeps
//      the existing per-stick weighted logistic-regression construction
//      (nu_k(x2_i) = P(z_i = k | z_i >= k, x2_i)) but adds an L2 shrinkage
//      prior pulling each stick's intercept toward the value implied by a
//      sparse Dirichlet draw pi0 (converted to the equivalent
//      stick-breaking proportion v0_k), so the x2-adjusted weights cannot
//      drift away from the shared sparse baseline without real evidence in
//      the data. When a stick has too little data to estimate a slope, it
//      now falls back to the sparse baseline v0_k directly (rather than an
//      unconditional empirical mean), which is itself pulled toward zero by
//      alpha0's sparsity -- the desired behavior.
//
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// ---------------------------------------------------------------------------
// 1. Vectorized sampling of cluster indicators z_i from probabilities
//    (unchanged from dpmm_helpers.cpp)
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
// 2. Sparse Dirichlet mixture weights (naive arm; also the shared baseline
//    used by the x2-informed arm's shrinkage prior below)
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
// 3. x2-informed per-stick logistic regression, shrunk toward the sparse
//    Dirichlet baseline pi0
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix compute_stick_breaking_sparse_cpp(const IntegerVector& z,
                                                 const NumericVector& x2,
                                                 const NumericVector& w,
                                                 int L,
                                                 const NumericVector& pi0,
                                                 double tau) {
    int n = z.size();
    NumericMatrix pi_mat(n, L);
    std::vector<double> remaining(n, 1.0);

    // Deterministic stick-breaking decomposition of the sparse Dirichlet
    // draw pi0, giving the shrinkage target v0_k for each stick.
    std::vector<double> v0(L, 0.0);
    double remaining0 = 1.0;
    for (int k = 0; k < L - 1; ++k) {
        double vk = pi0[k] / std::max(remaining0, 1e-300);
        vk = std::min(std::max(vk, 1e-6), 1.0 - 1e-6);
        v0[k] = vk;
        remaining0 *= (1.0 - vk);
    }

    double inv_tau2 = 1.0 / (tau * tau);

    for (int k = 1; k < L; ++k) { // k is 1-based stick index: 1 to L-1
        double v0k = v0[k - 1];
        double prior_mean0 = std::log(v0k / (1.0 - v0k));

        int m = 0;
        double sum_y = 0.0;
        for (int i = 0; i < n; ++i) {
            if (z[i] >= k) {
                m++;
                if (z[i] == k) sum_y += 1.0;
            }
        }

        std::vector<double> nu_k(n);
        if (sum_y < 0.5 || sum_y > m - 0.5 || m < 5) {
            // Not enough evidence to estimate an x2 slope: fall back to the
            // sparse baseline directly (itself pulled toward 0 by alpha0).
            for (int i = 0; i < n; ++i) nu_k[i] = v0k;
        } else {
            std::vector<double> sub_y(m), sub_x(m), sub_w(m);
            int idx = 0;
            for (int i = 0; i < n; ++i) {
                if (z[i] >= k) {
                    sub_y[idx] = (z[i] == k) ? 1.0 : 0.0;
                    sub_x[idx] = x2[i];
                    sub_w[idx] = w[i];
                    idx++;
                }
            }

            double b0 = prior_mean0; // start at the shrinkage target
            double b1 = 0.0;

            for (int iter = 0; iter < 12; ++iter) {
                double g0 = 0.0, g1 = 0.0;
                double h00 = 0.0, h01 = 0.0, h11 = 0.0;

                for (int j = 0; j < m; ++j) {
                    double eta = b0 + b1 * sub_x[j];
                    double mu = 1.0 / (1.0 + std::exp(-eta));
                    mu = std::max(1e-12, std::min(1.0 - 1e-12, mu));
                    double W = sub_w[j] * mu * (1.0 - mu);
                    double r = sub_w[j] * (sub_y[j] - mu);

                    g0 += r;
                    g1 += r * sub_x[j];
                    h00 += W;
                    h01 += W * sub_x[j];
                    h11 += W * sub_x[j] * sub_x[j];
                }

                // L2 shrinkage of the intercept toward prior_mean0 (the
                // sparse-Dirichlet-implied baseline); slope is unpenalized.
                g0 -= (b0 - prior_mean0) * inv_tau2;
                h00 += inv_tau2;

                double det = h00 * h11 - h01 * h01;
                if (std::abs(det) < 1e-12 || !std::isfinite(det)) break;

                double delta0 = (h11 * g0 - h01 * g1) / det;
                double delta1 = (-h01 * g0 + h00 * g1) / det;

                b0 += delta0;
                b1 += delta1;

                if (std::abs(delta0) + std::abs(delta1) < 1e-6) break;
            }

            for (int i = 0; i < n; ++i) {
                double eta = b0 + b1 * x2[i];
                nu_k[i] = 1.0 / (1.0 + std::exp(-eta));
            }
        }

        for (int i = 0; i < n; ++i) {
            pi_mat(i, k - 1) = nu_k[i] * remaining[i];
            remaining[i] *= (1.0 - nu_k[i]);
        }
    }

    for (int i = 0; i < n; ++i) pi_mat(i, L - 1) = remaining[i];

    return pi_mat;
}

// ---------------------------------------------------------------------------
// 4. Stratified PSU-level Sandwich Bootstrap Covariance (J_hat)
//    (unchanged from dpmm_helpers.cpp)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat fast_sandwich_bootstrap_cpp(const arma::mat& psu_scores,
                                      const Rcpp::List& stratum_psu_indices,
                                      int B) {
    int Kd = psu_scores.n_cols;
    int H = stratum_psu_indices.size();
    arma::mat score_boot(B, Kd, arma::fill::zeros);

    for (int b = 0; b < B; ++b) {
        arma::rowvec b_score(Kd, arma::fill::zeros);
        for (int h = 0; h < H; ++h) {
            Rcpp::IntegerVector psus_in_h = stratum_psu_indices[h];
            int nh = psus_in_h.size();
            if (nh == 0) continue;
            for (int r = 0; r < nh; ++r) {
                int pick = (int)(R::runif(0.0, 1.0) * nh);
                if (pick >= nh) pick = nh - 1;
                int psu_row = psus_in_h[pick];
                b_score += psu_scores.row(psu_row);
            }
        }
        score_boot.row(b) = b_score;
    }

    arma::mat J_hat = arma::cov(score_boot);
    return J_hat;
}

// ---------------------------------------------------------------------------
// 5. Stratified Delete-One-PSU Jackknife Covariance (J_hat)
// ---------------------------------------------------------------------------
// For a LINEAR (additive-over-PSU) score total T = sum_h sum_i t_hi, the
// delete-one-PSU jackknife variance estimator reduces algebraically to the
// closed-form Rao-Wu "ultimate cluster" with-replacement variance:
//   J_hat = sum_h [n_h/(n_h-1)] * sum_{i in h} (t_hi - tbar_h)(t_hi - tbar_h)'
// This is exactly the B -> infinity limit of fast_sandwich_bootstrap_cpp's
// Monte Carlo bootstrap above (same expectation, zero resampling noise, and
// no B hyperparameter). Strata with a single sampled PSU (n_h < 2) cannot be
// jackknifed and contribute zero -- same degenerate behavior the bootstrap
// has there (always resampling the one available PSU).
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
