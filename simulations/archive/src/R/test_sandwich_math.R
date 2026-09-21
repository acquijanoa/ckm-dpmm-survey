# =============================================================================
# Diagnostic: Inspect Sandwich Matrix and Hessian Formulation
# =============================================================================

suppressPackageStartupMessages({
  library(MASS)
})

res <- readRDS("results/result_informed_00001.rds")

cat("========================================================================\n")
cat(" DIAGNOSTIC AUDIT OF GODAMBE SANDWICH FORMULATION\n")
cat("========================================================================\n")

cat("\n1. Current (Stored) Sandwich Covariance vs Naive Covariance:\n")
cat("Diag of Sigma_sandwich (stored in RDS):\n")
print(diag(res$Sigma_sandwich))
cat("Implied Sandwich SEs:\n")
print(sqrt(diag(res$Sigma_sandwich)))

cat("\nDiag of Sigma_naive (stored in RDS):\n")
print(diag(res$Sigma_naive))
cat("Implied Naive SEs:\n")
print(sqrt(diag(res$Sigma_naive)))

cat("\nRatio of Sandwich SE to Naive SE:\n")
print(sqrt(diag(res$Sigma_sandwich)) / sqrt(diag(res$Sigma_naive)))

pi_hat <- res$pi_hat
Kd <- length(pi_hat) - 1
pi_sub <- pi_hat[1:Kd]

# Current formulation of H:
# H = (diag(pi) - outer(pi, pi)) * sum(w)
# In our simulations, sum(w) = n = 2000
sum_w <- 2000
V_multinomial <- diag(pi_sub) - outer(pi_sub, pi_sub)
H_current <- V_multinomial * sum_w
H_inv_current <- solve(H_current)

cat("\n2. Comparing H matrices:\n")
cat("Eigenvalues of V_multinomial:\n")
print(eigen(V_multinomial)$values)

cat("\nDiagonal of H_inv_current * sum_w (Inflation factor on scores):\n")
print(diag(H_inv_current) * sum_w)

cat("\nSquared inflation factor (Variance multiplier):\n")
print((diag(H_inv_current) * sum_w)^2)

cat("\n========================================================================\n")
