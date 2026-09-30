#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
//
// Rcpp kernels for dpmm_fit.R: allocation sampling, sparse Dirichlet
// mixture weights, and the stratified delete-one-PSU jackknife.

using namespace Rcpp;

// z_i ~ Categorical(probs[i, ]) by inverse CDF; returns 1-based labels.
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
                chosen = k + 1;
                break;
            }
        }
        z[i] = chosen;
    }
    return z;
}

// pi ~ Dirichlet(alpha0/L + n_1^w, ..., alpha0/L + n_L^w), with weighted
// counts n_k^w = sum_{z_i = k} w_i, via normalized Gamma draws.
// [[Rcpp::export]]
NumericVector sparse_dirichlet_weights_cpp(const IntegerVector& z,
                                            const NumericVector& w,
                                            int L, double alpha0) {
    int n = z.size();
    std::vector<double> n_k(L + 1, 0.0);  // indexed by 1-based label
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

// Delete-one-PSU jackknife covariance of a score total. For a total that is
// linear in the PSU totals t_hi this equals the ultimate-cluster estimator
//   J = sum_h n_h/(n_h - 1) sum_i (t_hi - tbar_h)(t_hi - tbar_h)'.
// psu_scores: one row per PSU; stratum_psu_indices: 0-based rows per stratum.
// Strata with a single PSU contribute zero.
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
