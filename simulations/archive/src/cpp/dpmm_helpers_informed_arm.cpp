#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
//
// Archived: the x2-informed "informed arm" stick-breaking mechanism and the
// bootstrap-based sandwich it was originally paired with. Both belong to the
// original PPS1/PPS3 validation design (informed vs. naive arm comparison),
// which lives in the now-archived simulations/archive/src/R/
// dpmm_survey_simulation.R. Neither function is called by the currently
// active src/R/dpmm_fit.R pipeline (which uses only the naive-arm sparse
// Dirichlet weights, sparse_dirichlet_weights_cpp in src/cpp/dpmm_helpers.cpp,
// and the closed-form jackknife_sandwich_cpp in place of the bootstrap below).
// Split out of src/cpp/dpmm_helpers.cpp to keep the active helper file down
// to only what the active pipeline actually calls.
//
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// ---------------------------------------------------------------------------
// x2-informed per-stick logistic regression, shrunk toward the sparse
// Dirichlet baseline pi0
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
// Standard Unpenalized Stick-Breaking IRLS Logistic Regression
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix compute_stick_breaking_cpp(const IntegerVector& z,
                                         const NumericVector& x2,
                                         const NumericVector& w,
                                         int L) {
    int n = z.size();
    NumericMatrix pi_mat(n, L);
    std::vector<double> remaining(n, 1.0);

    for (int k = 1; k < L; ++k) {
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
            double mean_val = (sum_y + 1e-3 * m) / (double)m;
            for (int i = 0; i < n; ++i) nu_k[i] = mean_val;
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

            double p_init = std::max(1e-4, std::min(1.0 - 1e-4, sum_y / (double)m));
            double b0 = std::log(p_init / (1.0 - p_init));
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
// Stratified PSU-level Sandwich Bootstrap Covariance (J_hat)
// ---------------------------------------------------------------------------
// Superseded by the closed-form jackknife_sandwich_cpp in
// src/cpp/dpmm_helpers.cpp, which is the B -> infinity limit of this Monte
// Carlo bootstrap for the linear score total used here (same expectation,
// zero resampling noise, no B hyperparameter).
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
