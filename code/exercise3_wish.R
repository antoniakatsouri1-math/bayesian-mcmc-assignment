# ============================================================
# ΑΣΚΗΣΗ 3 — Wish Dataset (Kaggle) — Binomial Model
# ============================================================

library(dplyr)
library(R2jags)
library(coda)

# --- Φόρτωση Δεδομένων ---
wish_data <- read.csv("data/summer-products-with-rating-and-performance_2020-08.csv")

# --- Προεπεξεργασία ---
summary_table <- wish_data %>%
  mutate(
    badge_category = case_when(
      badge_product_quality == 1 ~ "Excellent Quality",
      badge_local_product   == 1 ~ "Top Merchant",
      TRUE                       ~ "No Badge"
    )
  ) %>%
  group_by(badge_category) %>%
  summarise(
    n_total = n(),
    y_ad    = sum(uses_ad_boosts == 1, na.rm = TRUE)
  ) %>%
  arrange(badge_category)

print(summary_table)

# --- JAGS Model ---
# j=1: No Badge         ~ Beta(1,1)
# j=2: Excellent Quality ~ Beta(1,1)
# j=3: Top Merchant     ~ Beta(40,10)

model_code <- "
model {
  for (i in 1:3) {
    y[i] ~ dbin(theta[i], n[i])
  }
  theta[1] ~ dbeta(1,  1)
  theta[2] ~ dbeta(1,  1)
  theta[3] ~ dbeta(40, 10)
}
"
writeLines(model_code, "wish_model.jags")

jags_data <- list(
  n = c(1433, 117, 23),
  y = c(616,   51, 14)
)

set.seed(123)
jags_fit <- jags(
  data               = jags_data,
  inits              = function() list(theta = c(0.5, 0.5, 0.5)),
  parameters.to.save = c("theta"),
  model.file         = "wish_model.jags",
  n.chains           = 3,
  n.iter             = 20000,
  n.burnin           = 5000,
  n.thin             = 2
)
print(jags_fit)

# --- Διαγνωστικά ---
mcmc_samples <- as.mcmc(jags_fit)

png("plots/1_Trace_and_Density_Plots.png", width = 800, height = 600)
plot(mcmc_samples)
dev.off()

png("plots/2_Geweke_Diagnostic_Plot.png", width = 800, height = 500)
geweke.plot(mcmc_samples)
dev.off()

png("plots/3_Gelman_Rubin_Plot.png", width = 800, height = 500)
gelman.plot(mcmc_samples)
dev.off()

png("plots/4_Autocorrelation_Plots.png", width = 800, height = 600)
autocorr.plot(mcmc_samples)
dev.off()

cat("Διαγνωστικά γραφήματα αποθηκεύτηκαν.\n")
