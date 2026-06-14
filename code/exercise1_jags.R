# ============================================================
# ΑΣΚΗΣΗ 1 — Πολλαπλή Γραμμική Παλινδρόμηση
# Μέρος Γ: Μπεϋζιανή Εκτίμηση μέσω JAGS
# ============================================================

library(MASS)
library(rjags)
library(coda)

# --- Προσομοίωση Δεδομένων ---
set.seed(123)
n <- 50
p <- 15

X_1_10 <- mvrnorm(n = n, mu = rep(0, 10), Sigma = diag(10))

X <- matrix(NA, nrow = n, ncol = p + 1)
X[, 1] <- 1
X[, 2:11] <- X_1_10

for(i in 1:n) {
  mu_x <- 0.3*X[i,2] + 0.5*X[i,3] + 0.7*X[i,4] +
          0.9*X[i,5] + 1.5*X[i,6]
  for(j in 12:16) {
    X[i, j] <- rnorm(1, mean = mu_x, sd = 1)
  }
}

Y <- rep(NA, n)
for(i in 1:n) {
  mu_y <- 4 + 2*X[i,2] - X[i,6] + 1.5*X[i,8] +
          X[i,12] + 0.5*X[i,14]
  Y[i] <- rnorm(1, mean = mu_y, sd = 2.5)
}

true_betas <- c(4, 2, 0, 0, 0, -1, 0, 1.5, 0, 0, 0, 1, 0, 0.5, 0, 0)

# --- OLS ---
ols_model   <- lm(Y ~ X - 1)
ols_summary <- summary(ols_model)

# --- JAGS Model ---
jags_model_string <- "
model {
  for (i in 1:n) {
    Y[i] ~ dnorm(mu[i], tau)
    mu[i] <- inprod(X[i,], beta[])
  }
  for (j in 1:16) {
    beta[j] ~ dnorm(0, tau * 0.0001)
  }
  tau    ~ dgamma(0.0001, 0.0001)
  sigma2 <- 1 / tau
}
"

model <- jags.model(
  textConnection(jags_model_string),
  data     = list(Y = Y, X = X, n = n),
  n.chains = 2
)
update(model, 2000)

samples <- coda.samples(
  model          = model,
  variable.names = c("beta", "sigma2"),
  n.iter         = 10000
)

samples_matrix <- as.matrix(samples)
param_names    <- colnames(samples_matrix)
beta_idx       <- grep("^beta",   param_names)
sigma2_idx     <- grep("^sigma2", param_names)

# --- Αποτελέσματα ---
bayes_mean <- colMeans(samples_matrix)
bayes_sd   <- apply(samples_matrix, 2, sd)

comparison_table <- data.frame(
  True_Value   = c(true_betas,           2.5^2),
  OLS_Estimate = c(coef(ols_model),      ols_summary$sigma^2),
  Bayes_Mean   = c(bayes_mean[beta_idx], bayes_mean[sigma2_idx]),
  Bayes_SD     = c(bayes_sd[beta_idx],   bayes_sd[sigma2_idx])
)
rownames(comparison_table) <- c(paste0("beta_", 0:15), "sigma2")

cat("--- Gelman-Rubin Diagnostic ---\n")
print(gelman.diag(samples, multivariate = FALSE))
cat("\n--- Comparison Table ---\n")
print(round(comparison_table, 4))

# --- Γραφήματα ---
n_params <- length(param_names)
params_per_page <- 4

# Trace Plots
for(start in seq(1, n_params, by = params_per_page)) {
  end <- min(start + params_per_page - 1, n_params)
  idx <- start:end
  png(sprintf("plots/trace_plot_%02d_%02d.png", start, end),
      width = 1200, height = 800, res = 120)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for(k in idx) {
    chain1 <- samples[[1]][, k]
    chain2 <- samples[[2]][, k]
    ylim   <- range(c(chain1, chain2))
    plot(chain1, type = "l", col = "steelblue",
         main = paste("Trace:", param_names[k]),
         xlab = "Iteration", ylab = "Value", ylim = ylim)
    lines(chain2, col = "tomato")
    legend("topright", legend = c("Chain 1","Chain 2"),
           col = c("steelblue","tomato"), lty = 1, cex = 0.7)
  }
  dev.off()
}

# Density Plots
for(start in seq(1, n_params, by = params_per_page)) {
  end <- min(start + params_per_page - 1, n_params)
  idx <- start:end
  png(sprintf("plots/density_plot_%02d_%02d.png", start, end),
      width = 1200, height = 800, res = 120)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for(k in idx) {
    chain1 <- density(samples[[1]][, k])
    chain2 <- density(samples[[2]][, k])
    xlim   <- range(c(chain1$x, chain2$x))
    ylim   <- range(c(chain1$y, chain2$y))
    plot(chain1, col = "steelblue", lwd = 2,
         main = paste("Density:", param_names[k]),
         xlab = "Value", xlim = xlim, ylim = ylim)
    lines(chain2, col = "tomato", lwd = 2)
    legend("topright", legend = c("Chain 1","Chain 2"),
           col = c("steelblue","tomato"), lty = 1, cex = 0.7)
  }
  dev.off()
}

# Boxplots
for(start in seq(1, n_params, by = params_per_page)) {
  end <- min(start + params_per_page - 1, n_params)
  idx <- start:end
  png(sprintf("plots/boxplot_%02d_%02d.png", start, end),
      width = 1200, height = 800, res = 120)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  for(k in idx) {
    boxplot(samples_matrix[, k],
            main   = paste("Boxplot:", param_names[k]),
            ylab   = "Τιμή",
            col    = "lightblue",
            border = "navy")
    abline(h = mean(samples_matrix[, k]),
           col = "red", lwd = 2, lty = 2)
  }
  dev.off()
}

cat("Γραφήματα αποθηκεύτηκαν στον φάκελο plots/\n")
