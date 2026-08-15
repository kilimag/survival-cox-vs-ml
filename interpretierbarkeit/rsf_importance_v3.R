# ===================================================================
# BAUSTEIN 3: RSF Variable Importance (Hazard 0.3)
# -------------------------------------------------------------------
# Zeigt: RSF erkennt x2 als wichtig - genau die Variable, die das
# naive Cox mit HR ~ 1 übersieht. Permutation Importance über
# mehrere Iterationen gemittelt.
# ===================================================================

library(survival)
library(ranger)
library(ggplot2)
library(dplyr)

HAZARD_FACTOR <- 0.3

simulate_censoring <- function(n = 1000, cens_rate) {
  x1 <- runif(n, -2, 2); x2 <- rnorm(n, 0, 1)
  x3 <- rbinom(n, 1, 0.5); x4 <- rnorm(n, 0, 1)
  lp_base <- 1.5 * sin(pi * x1) + 1.2 * (x2^2) + (1.5 * x1 * x3)
  true_time <- ( -log(runif(n)) / (HAZARD_FACTOR * exp(lp_base + 1.5 * x3)) )^(1/(1.1 + 0.5 * x3))
  cens_time <- rexp(n, rate = cens_rate)
  time   <- pmin(true_time, cens_time, 5)
  status <- ifelse(true_time <= pmin(cens_time, 5), 1, 0)
  data.frame(time, status, x1, x2, x3, x4)
}

calibrate_rate <- function(target_cens, n_calib = 20000, seed = 1) {
  f <- function(rate) {
    set.seed(seed)
    d <- simulate_censoring(n = n_calib, cens_rate = rate)
    mean(d$status == 0) - target_cens
  }
  uniroot(f, interval = c(1e-4, 50), tol = 1e-4)$root
}

form <- Surv(time, status) ~ x1 + x2 + x3 + x4

cens_target <- 0.20
n_iter      <- 100

cat("Kalibriere Rate fuer", cens_target*100, "% ...\n")
rate <- calibrate_rate(cens_target)
cat("Rate:", round(rate, 4), "\n")

# -------------------------------------------------------------------
# Permutation Importance über Iterationen
# -------------------------------------------------------------------
imp_mat <- matrix(NA, nrow = n_iter, ncol = 4,
                  dimnames = list(NULL, c("x1","x2","x3","x4")))

cat("\nBerechne Importance über", n_iter, "Iterationen...\n")
for (i in 1:n_iter) {
  if (i %% 10 == 0) cat("  Iteration", i, "\n")
  set.seed(1000 + i)
  d <- simulate_censoring(n = 1000, cens_rate = rate)
  m <- ranger(form, data = d, num.trees = 500,
              importance = "permutation")
  imp_mat[i, ] <- importance(m)[c("x1","x2","x3","x4")]
}

# Zusammenfassen
imp_data <- data.frame(
  Variable = c("x1","x2","x3","x4"),
  Importance = colMeans(imp_mat, na.rm = TRUE),
  SD = apply(imp_mat, 2, sd, na.rm = TRUE),
  Rolle = c("Signal","Signal","Signal","Rauschen")
)

cat("\n=== RSF Permutation Importance (gemittelt) ===\n")
print(imp_data, row.names = FALSE)

# -------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------
imp_data$Variable <- factor(imp_data$Variable,
                            levels = imp_data$Variable[order(imp_data$Importance)])
plot_imp <- ggplot(imp_data, aes(x = Importance, y = Variable, fill = Rolle)) +
  geom_col(width = 0.6) +
  geom_errorbarh(aes(xmin = pmax(0, Importance - SD), xmax = Importance + SD),
                 height = 0.2, color = "grey30") +
  geom_vline(xintercept = 0, color = "grey50") +
  scale_fill_manual(values = c("Signal"="#4DAF4A", "Rauschen"="#999999")) +
  theme_minimal() +
  labs(title = "Random Survival Forest: Permutation Variable Importance",
       subtitle = "RSF erkennt x2 als wichtig - genau die Variable, die das naive Cox mit HR ~ 1 übersieht",
       x = "Permutation Importance (gemittelt, +/- SD)", y = "Kovariate",
       fill = "Wahre Rolle") +
  theme(text = element_text(size = 12), legend.position = "bottom")

ggsave("Plot_v3_RSF_Importance:withoutcap.png", plot = plot_imp, width = 10, height = 6, dpi = 300)
print(plot_imp)

write.csv(imp_data, "Ergebnisse_v3_RSF_Importance.csv", row.names = FALSE)
cat("\nGespeichert: Plot_v3_RSF_Importance.png\n")
