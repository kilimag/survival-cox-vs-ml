# ===================================================================
# Schoenfeld-Residuen & PH-Detection (Hazard 0.3)
# -------------------------------------------------------------------
# Teil A: Einzelner Schoenfeld-Fit (Illustration, nur x1)
# Teil B: PH-Detection-Quote über 100 Iterationen, alle Variablen,
#         20% vs 90% Zensierung -> der 100%->x% Befund
# ===================================================================

library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)

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

# -------------------------------------------------------------------
# TEIL B zuerst: PH-Detection-Quote über Iterationen
# -------------------------------------------------------------------
cens_levels <- c(0.20, 0.90)
n_iter      <- 100

cat("Kalibriere Raten...\n")
rates <- sapply(cens_levels, calibrate_rate)
names(rates) <- paste0(cens_levels*100, "%")
print(round(rates, 4))

detection <- data.frame()

for (k in seq_along(cens_levels)) {
  cens <- cens_levels[k]; rate <- rates[k]
  cat("\nZensierung", cens*100, "% ...\n")
  
  for (i in 1:n_iter) {
    set.seed(round(cens*1000) + i)
    d <- simulate_censoring(n = 1000, cens_rate = rate)
    
    res <- tryCatch({
      m <- coxph(form, data = d)
      zph <- cox.zph(m)
      # p-Werte pro Variable + GLOBAL
      pvals <- zph$table[, "p"]
      data.frame(
        Zensierung = cens,
        Variable = rownames(zph$table),
        p = pvals,
        verletzt = pvals < 0.05
      )
    }, error = function(e) NULL)
    
    if (!is.null(res)) detection <- rbind(detection, res)
  }
}

# Detection-Quote je Variable und Zensierung
quote <- detection %>%
  group_by(Zensierung, Variable) %>%
  summarise(Detection_Rate = mean(verletzt, na.rm = TRUE),
            n_valid = sum(!is.na(verletzt)), .groups = "drop")

cat("\n=== PH-Detection-Quote (Anteil Simulationen mit p < 0.05) ===\n")
quote_wide <- quote %>%
  select(Variable, Zensierung, Detection_Rate) %>%
  pivot_wider(names_from = Zensierung, values_from = Detection_Rate)
print(as.data.frame(quote_wide), row.names = FALSE)

# -------------------------------------------------------------------
# Plot B: Detection-Quote
# -------------------------------------------------------------------
quote$Variable <- factor(quote$Variable,
                         levels = c("x1","x2","x3","x4","GLOBAL"))
quote$Zens_label <- paste0(quote$Zensierung*100, "%")

plot_det <- ggplot(quote, aes(x = Variable, y = Detection_Rate,
                              fill = Zens_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "grey40") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0,1)) +
  scale_fill_manual(values = c("20%"="#377EB8", "90%"="#E41A1C")) +
  theme_minimal() +
  labs(title = "Erkennung der PH-Verletzung durch cox.zph()",
       subtitle = "Anteil der Simulationen mit signifikanter Verletzung (p < 0.05)",
       x = "Kovariate", y = "Erkennungsquote", fill = "Zensierung") +
  theme(text = element_text(size = 12), legend.position = "bottom")

print(plot_det)

ggsave("Plot_v3_PH_Detection.png", plot = plot_det, width = 10, height = 6, dpi = 300)
print(plot_det)

# -------------------------------------------------------------------
# TEIL A: Einzelner Schoenfeld-Fit zur Illustration (x1 vs x3)
# -------------------------------------------------------------------
set.seed(202)
rate20 <- rates[1]
d_ill <- simulate_censoring(n = 1000, cens_rate = rate20)
m_ill <- coxph(form, data = d_ill)
zph_ill <- cox.zph(m_ill)

cat("\n=== Einzelner Fit (20% Zensierung, Seed 202) ===\n")
print(zph_ill$table)

# Schoenfeld-Plot fuer x1 speichern
png("Plot_v3_Schoenfeld_x1.png", width = 900, height = 500, res = 110)
plot(zph_ill[1], main = paste0("Skaliertes Schoenfeld-Residuum: x1 (p = ",
                               round(zph_ill$table["x1","p"], 4), ")"))
abline(h = 0, lty = 3, col = "grey50")
dev.off()

cat("\nGespeichert:\n  Plot_v3_PH_Detection.png\n  Plot_v3_Schoenfeld_x1.png\n")

write.csv(quote, "Ergebnisse_v3_PH_Detection.csv", row.names = FALSE)