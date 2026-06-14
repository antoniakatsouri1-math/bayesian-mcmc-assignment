# ============================================================
# ΑΣΚΗΣΗ 4 — Jeffreys Prior & Metropolis-Hastings (Poisson)
# ============================================================

set.seed(456)
n     <- 10
sum_y <- 35

# --- Log-posterior υπό Jeffreys prior ---
# p(lambda) ∝ lambda^(-1/2)
# Ύστερη: Gamma(35.5, 10)
log_posterior <- function(lambda, n, sum_y) {
  if (lambda <= 0) return(-Inf)
  return(-n * lambda + (sum_y - 0.5) * log(lambda))
}

# --- Random Walk Metropolis-Hastings ---
N_iter          <- 20000
lambda_chain    <- numeric(N_iter)
lambda_chain[1] <- 1.0
sigma_prop      <- 1.65
accepted        <- 0

for (t in 2:N_iter) {
  current  <- lambda_chain[t-1]
  proposed <- rnorm(1, mean = current, sd = sigma_prop)

  log_alpha <- log_posterior(proposed, n, sum_y) -
               log_posterior(current,  n, sum_y)

  if (log(runif(1)) < log_alpha) {
    lambda_chain[t] <- proposed
    accepted        <- accepted + 1
  } else {
    lambda_chain[t] <- current
  }
}

acc_rate <- accepted / (N_iter - 1) * 100
cat("Ποσοστό Αποδοχής:", round(acc_rate, 2), "%\n")

# --- Burn-in ---
burn_in     <- 2000
final_chain <- lambda_chain[(burn_in + 1):N_iter]

# --- Γραφήματα ---
png("plots/MH_Plots_Corrected_1.png", width = 900, height = 450)
par(mfrow = c(1, 2))
hist(final_chain, breaks = 40, col = "lightblue",
     main = "Ιστόγραμμα του λ", xlab = "λ", probability = TRUE)
lines(density(final_chain), col = "darkblue", lwd = 2)
plot(lambda_chain, type = "l", col = "gray40",
     main = "Trace Plot του λ",
     xlab = "Επαναλήψεις", ylab = "λ")
abline(v = burn_in, col = "red", lwd = 2, lty = 2)
dev.off()

ergodic_mean <- cumsum(final_chain) / seq_along(final_chain)

png("plots/MH_Plots_Corrected_2.png", width = 900, height = 450)
par(mfrow = c(1, 2))
plot(ergodic_mean, type = "l", col = "darkgreen", lwd = 2,
     main = "Ergodic Mean Plot",
     xlab = "Επαναλήψεις (μετά burn-in)", ylab = "Τρεχούσα Μέση")
acf(final_chain, main = "Διάγραμμα Αυτοσυσχετίσεων (ACF)",
    lag.max = 50)
dev.off()

# --- Στατιστικά ---
cat("\n=== Περιγραφικά Στατιστικά ===\n")
print(summary(final_chain))
cat("SD   =", round(sd(final_chain), 4), "\n")
ci <- quantile(final_chain, c(0.025, 0.975))
cat("95% CI: [", round(ci[1], 4), ",", round(ci[2], 4), "]\n")
