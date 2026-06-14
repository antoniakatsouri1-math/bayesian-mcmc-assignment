# ============================================================
# ΑΣΚΗΣΗ 1 — Πολλαπλή Γραμμική Παλινδρόμηση
# Μέρος Δ: Gibbs Sampler (υλοποίηση από μηδέν)
# ============================================================

library(mvtnorm)

gibbs_regression <- function(X, Y, M_diag = 0.0001,
                              n_iter = 20000, burn_in = 5000) {
  n        <- nrow(X)
  p_plus_1 <- ncol(X)

  X_mat <- as.matrix(X)
  Y_vec <- as.numeric(Y)

  M        <- diag(M_diag, p_plus_1)
  a        <- 0.0001
  b        <- 0.0001
  XTX      <- t(X_mat) %*% X_mat
  XTY      <- t(X_mat) %*% Y_vec
  V_beta   <- solve(M + XTX)
  beta_hat <- as.numeric(V_beta %*% XTY)

  beta_samples   <- matrix(0, nrow = n_iter, ncol = p_plus_1)
  sigma2_samples <- numeric(n_iter)
  current_sigma2 <- 1.0

  for (t in 1:n_iter) {
    curr_Sigma   <- current_sigma2 * V_beta
    current_beta <- as.numeric(rmvnorm(1, mean = beta_hat,
                                       sigma = curr_Sigma))
    beta_samples[t, ] <- current_beta

    resid          <- Y_vec - as.numeric(X_mat %*% current_beta)
    shape_sig      <- a + n / 2
    scale_sig      <- b + 0.5 * sum(resid^2)
    current_sigma2 <- 1 / rgamma(1, shape = shape_sig,
                                  rate  = scale_sig)
    sigma2_samples[t] <- current_sigma2
  }

  keep_idx <- (burn_in + 1):n_iter
  return(list(
    beta   = beta_samples[keep_idx, ],
    sigma2 = sigma2_samples[keep_idx]
  ))
}

# --- Εκτέλεση (χρησιμοποιεί X, Y, true_betas από exercise1_jags.R) ---
set.seed(123)
gibbs_results <- gibbs_regression(X, Y, M_diag = 0.0001,
                                   n_iter = 20000, burn_in = 5000)

gibbs_beta_mean   <- colMeans(gibbs_results$beta)
gibbs_beta_sd     <- apply(gibbs_results$beta, 2, sd)
gibbs_sigma2_mean <- mean(gibbs_results$sigma2)
gibbs_sigma2_sd   <- sd(gibbs_results$sigma2)

comparison_table <- data.frame(
  True_Value   = c(true_betas,           2.5^2),
  OLS_Estimate = c(coef(ols_model),      ols_summary$sigma^2),
  OLS_SE       = c(ols_summary$coefficients[, 2], NA),
  Gibbs_Mean   = c(gibbs_beta_mean,      gibbs_sigma2_mean),
  Gibbs_SD     = c(gibbs_beta_sd,        gibbs_sigma2_sd)
)
rownames(comparison_table) <- c(paste0("beta_", 0:15), "sigma2")

cat("--- Gibbs Sampler Results ---\n")
print(round(comparison_table, 4))
