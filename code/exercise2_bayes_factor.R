# ============================================================
# ΑΣΚΗΣΗ 2 — Μπεϋζιανός Έλεγχος Υπόθεσης (Poisson)
# ============================================================

# --- Δεδομένα ---
n   <- 10
y   <- 35

# --- Ερώτημα Α: Bayes Factor ---
# H0: lambda = 2
# H1: lambda ~ Gamma(2, 1)

# Υπολογισμός σε λογαριθμική κλίμακα
log_BF10 <- lfactorial(36) + 20 - 37 * log(11) - 35 * log(2)
BF10     <- exp(log_BF10)

cat("=== Bayes Factor ===\n")
cat("ln(BF10) =", round(log_BF10, 4), "\n")
cat("BF10     =", round(BF10, 4), "\n")
cat("Ερμηνεία: Ισχυρή ένδειξη υπέρ H1 (κλίμακα Jeffreys)\n\n")

# --- Ερώτημα Β: Ύστερη Προβλεπτική Κατανομή ---
# Ύστερη: lambda|y ~ Gamma(37, 11)
alpha_post <- y + 2      # = 37
beta_post  <- n + 1      # = 11

cat("=== Ύστερη Κατανομή ===\n")
cat("lambda | y ~ Gamma(", alpha_post, ",", beta_post, ")\n\n")

# Μπεϋζιανή προβλεπτική: x11|y ~ NegBin(37, 11/12)
E_bayes   <- alpha_post / beta_post
Var_bayes <- E_bayes + alpha_post / beta_post^2

cat("=== Προβλεπτική Κατανομή ===\n")
cat("x11|y ~ NegBin(r=37, p=11/12)\n")
cat("E(x11|y)   =", round(E_bayes,   4), "\n")
cat("Var(x11|y) =", round(Var_bayes, 4), "\n\n")

# Κλασική
lambda_mle <- y / n
cat("=== Κλασική Προσέγγιση ===\n")
cat("lambda_MLE =", lambda_mle, "\n")
cat("E_class    =", lambda_mle, "\n")
cat("Var_class  =", lambda_mle, "\n\n")

# P(x11 > 5 | y)
prob_le5 <- pnbinom(5, size = alpha_post,
                    prob = beta_post / (beta_post + 1))
prob_gt5 <- 1 - prob_le5

cat("=== P(x11 > 5 | y) ===\n")
cat("P(x11 <= 5 | y) =", round(prob_le5, 4), "\n")
cat("P(x11 >  5 | y) =", round(prob_gt5, 4), "\n")
