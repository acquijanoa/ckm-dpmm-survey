#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

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
// 2. 2-parameter Weighted Logistic Regression for Stick Breaking
// ---------------------------------------------------------------------------
// Fits logistic regression y ~ x2 with weights w using 2x2 analytic IRLS
// [[Rcpp::export]]
NumericMatrix compute_stick_breaking_cpp(const IntegerVector& z, 
                                         const NumericVector& x2, 
                                         const NumericVector& w, 
                                         int L) {
    int n = z.size();
    NumericMatrix pi_mat(n, L);
    std::vector<double> remaining(n, 1.0);

    for (int k = 1; k < L; ++k) { // k is 1-based stick index: 1 to L-1
        // Identify subset where z >= k
        int m = 0;
        double sum_y = 0.0;
        double sum_w = 0.0;
        for (int i = 0; i < n; ++i) {
            if (z[i] >= k) {
                m++;
                if (z[i] == k) sum_y += 1.0;
                sum_w += w[i];
            }
        }

        std::vector<double> nu_k(n);
        // Degenerate checks
        if (sum_y < 0.5 || sum_y > m - 0.5 || m < 5) {
            double mean_val = (sum_y + 1e-3 * m) / (double)m;
            for (int i = 0; i < n; ++i) {
                nu_k[i] = mean_val;
            }
        } else {
            // Extract subset
            std::vector<double> sub_y(m);
            std::vector<double> sub_x(m);
            std::vector<double> sub_w(m);
            int idx = 0;
            for (int i = 0; i < n; ++i) {
                if (z[i] >= k) {
                    sub_y[idx] = (z[i] == k) ? 1.0 : 0.0;
                    sub_x[idx] = x2[i];
                    sub_w[idx] = w[i];
                    idx++;
                }
            }

            // IRLS for 2-parameter logistic regression (intercept + slope)
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
                if (std::abs(det) < 1e-12 || !std::isfinite(det)) {
                    break;
                }

                double delta0 = (h11 * g0 - h01 * g1) / det;
                double delta1 = (-h01 * g0 + h00 * g1) / det;

                b0 += delta0;
                b1 += delta1;

                if (std::abs(delta0) + std::abs(delta1) < 1e-6) {
                    break;
                }
            }

            // Predict for all n
            for (int i = 0; i < n; ++i) {
                double eta = b0 + b1 * x2[i];
                nu_k[i] = 1.0 / (1.0 + std::exp(-eta));
            }
        }

        // Update pi_mat column k-1 (0-based)
        for (int i = 0; i < n; ++i) {
            pi_mat(i, k - 1) = nu_k[i] * remaining[i];
            remaining[i] *= (1.0 - nu_k[i]);
        }
    }

    // Last stick
    for (int i = 0; i < n; ++i) {
        pi_mat(i, L - 1) = remaining[i];
    }

    return pi_mat;
}

// ---------------------------------------------------------------------------
// 3. Stratified PSU-level Sandwich Bootstrap Covariance (J_hat)
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
