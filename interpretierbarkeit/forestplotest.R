# ===================================================================
# BAUSTEIN 1: Forest-Plot der Cox-HRs (Hazard 0.3) - BASIS-VERSION
# -------------------------------------------------------------------
# Zeigt die vom naiven Cox geschätzten Hazard Ratios, gemittelt
# über Iterationen. Kernbefund: x2 (starkes quadratisches Signal)
# bekommt HR ~ 1.00 -> Cox übersieht es.
#
# HINWEIS: Diese Basis-Version zeigt NUR die geschätzten HRs.
# Das Einzeichnen der "Wahrheit" (x3, x4) kommt als Erweiterung
# dazu im Skript "forest_plot_hr_v5_untereinander". 
#
# Zwei Zensierungsstufen (20% und 50%) - 50% war Prof-Wunsch.
# ===================================================================

library(survival)
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

cens_levels <- c(0.20, 0.50)   # 50% war Prof-Wunsch
n_iter      <- 100

cat("Kalibriere Raten...\n")
rates <- sapply(cens_levels, calibrate_rate)
names(rates) <- paste0(cens_levels*100, "%")
print(round(rates, 4))

# -------------------------------------------------------------------
# HRs ueber Iterationen sammeln
# -------------------------------------------------------------------
hr_data <- data.frame()

for (k in seq_along(cens_levels)) {
  cens <- cens_levels[k]; rate <- rates[k]
  cat("\nZensierung", cens*100, "% ...\n")
  
  coefs <- matrix(NA, nrow = n_iter, ncol = 4,
                  dimnames = list(NULL, c("x1","x2","x3","x4")))
  
  for (i in 1:n_iter) {
    set.seed(round(cens*1000) + i)
    d <- simulate_censoring(n = 1000, cens_rate = rate)
    m <- tryCatch(coxph(form, data = d), error = function(e) NULL)
    if (!is.null(m)) coefs[i, ] <- coef(m)
  }
  
  # Mittelwerte und CIs der HRs (exp der log-HR)
  for (v in c("x1","x2","x3","x4")) {
    log_hr <- coefs[, v]
    hr_data <- rbind(hr_data, data.frame(
      Zensierung = paste0(cens*100, "% Zensierung"),
      Variable = v,
      HR = exp(mean(log_hr, na.rm = TRUE)),
      HR_low = exp(quantile(log_hr, 0.025, na.rm = TRUE)),
      HR_high = exp(quantile(log_hr, 0.975, na.rm = TRUE)),
      Rolle = ifelse(v == "x4", "Rauschen", "Signal")
    ))
  }
}

cat("\n=== Geschaetzte Hazard Ratios (gemittelt) ===\n")
print(hr_data[, c("Zensierung","Variable","HR","HR_low","HR_high")], row.names = FALSE)

# -------------------------------------------------------------------
# Forest-Plot
# -------------------------------------------------------------------
hr_data$Variable <- factor(hr_data$Variable, levels = c("x4","x3","x2","x1"))

plot_forest <- ggplot(hr_data, aes(x = HR, y = Variable, color = Rolle)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3.5) +
  geom_errorbarh(aes(xmin = HR_low, xmax = HR_high), height = 0.2, linewidth = 0.8) +
  facet_wrap(~ Zensierung) +
  scale_color_manual(values = c("Signal"="#377EB8", "Rauschen"="#999999")) +
  theme_minimal() +
  labs(title = "Vom naiven Cox geschaetzte Hazard Ratios",
       subtitle = "Kernbefund: x2 (starkes quadratisches Signal) wird mit HR ~ 1 als irrelevant eingestuft",
       x = "Hazard Ratio (gemittelt, mit 95%-Bereich)", y = "Kovariate",
       color = "Wahre Rolle",
       caption = "Gestrichelt: HR = 1 (kein Effekt). x1/x2 wirken nichtlinear - ihr Effekt ist als HR nicht darstellbar.") +
  theme(text = element_text(size = 12), legend.position = "bottom")

ggsave("Plot_v3_ForestPlot.png", plot = plot_forest, width = 11, height = 6, dpi = 300)
print(plot_forest)

write.csv(hr_data, "Ergebnisse_v3_ForestHR.csv", row.names = FALSE)
cat("\nGespeichert: Plot_v3_ForestPlot.png\n")
