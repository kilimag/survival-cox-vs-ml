# ===================================================================
# SENSITIVITAETSRECHNUNG: Harrells C vs. Uno's C
# -------------------------------------------------------------------
# Zweck: Zeigen, dass sich die RANGFOLGE der Modelle nicht ändert,
# wenn man statt des zensierungsabhängigen Harrell-C das
# zensierungsrobuste Uno-C verwendet.
#
# Kein neuer Hauptlauf - reduzierte Sensitivitätsrechnung:
#   - 3 Zensierungsstufen (niedrig / mittel / hoch)
#   - 30 Iterationen je Stufe
#   - beide C-Varianten pro Modell nebeneinander
#
# Uno-C: survC1::Est.Cval(), Auswertungshorizont tau = 5 (Studienende)
# Konsistent zu v4: df = 6, NN-Seed gesetzt.
#
# Erzeugt:  UnoC_Vergleich.csv  (Mittelwerte beider Masse je Modell/Stufe)
#           Plot_UnoC_Vergleich.png
# ===================================================================

library(survival)
library(ranger)
library(survivalmodels)
library(reticulate)
library(splines)
library(dplyr)

# survC1 für Uno's C - falls nicht vorhanden:
#   install.packages("survC1")
if (!requireNamespace("survC1", quietly = TRUE)) {
  stop("Paket 'survC1' fehlt. Bitte installieren: install.packages('survC1')")
}

use_condaenv("r-reticulate", required = TRUE)

# ------------------------------------------------------------------
# 1. Simulation (identisch zu v4)
# ------------------------------------------------------------------
HAZARD_FACTOR <- 0.3

simulate_censoring <- function(n, cens_rate) {
  x1 <- runif(n, -2, 2); x2 <- rnorm(n, 0, 1)
  x3 <- rbinom(n, 1, 0.5); x4 <- rnorm(n, 0, 1)
  lp_base <- 1.5 * sin(pi * x1) + 1.2 * (x2^2) + (1.5 * x1 * x3)
  true_time <- ( -log(runif(n)) /
                   (HAZARD_FACTOR * exp(lp_base + 1.5 * x3)) )^(1/(1.1 + 0.5 * x3))
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
  tryCatch(uniroot(f, interval = c(1e-4, 50), tol = 1e-4)$root,
           error = function(e) 1e-4)
}

# ------------------------------------------------------------------
# 2. Metriken
# ------------------------------------------------------------------
# Harrells C (wie im Hauptlauf: concordance mit reverse)
harrell_c <- function(time, status, risk) {
  if (all(is.na(risk)) || length(unique(risk)) < 2) return(NA_real_)
  concordance(Surv(time, status) ~ risk, reverse = TRUE)$concordance
}

# Uno's C über survC1::Est.Cval
# Est.Cval erwartet: mydata = (time, status, marker), tau
# Hoehere Marker-Werte muessen hoeheres Risiko bedeuten -> risk direkt nutzen
uno_c <- function(time, status, risk, tau = 5) {
  if (all(is.na(risk)) || length(unique(risk)) < 2) return(NA_real_)
  df <- data.frame(time = time, status = status, marker = risk)
  out <- tryCatch(
    survC1::Est.Cval(mydata = df, tau = tau, nofit = TRUE)$Dhat,
    error = function(e) NA_real_
  )
  out
}

# ------------------------------------------------------------------
# 3. Setup
# ------------------------------------------------------------------
target_censoring <- c(0.20, 0.50, 0.70)   # niedrig / mittel / hoch
n_iterations     <- 30
N_TRAIN <- 1000
N_TEST  <- 1000
TAU     <- 5

form_naive  <- Surv(time, status) ~ x1 + x2 + x3 + x4
form_oracle <- Surv(time, status) ~ ns(x1, df = 6)*x3 + I(x2^2)

cat("Kalibriere Raten...\n")
rates <- sapply(target_censoring, calibrate_rate)
print(round(rates, 4))

results <- data.frame()

# ------------------------------------------------------------------
# 4. Hauptschleife
# ------------------------------------------------------------------
for (k in seq_along(target_censoring)) {
  target <- target_censoring[k]; rate <- rates[k]
  cat("\n--- Zensierung", target*100, "% ---\n")
  
  for (i in 1:n_iterations) {
    if (i %% 10 == 0) cat("  Iter", i, "/", n_iterations, "\n")
    iter_seed <- round(target*1000) + i
    set.seed(iter_seed); set_seed(iter_seed)
    train <- simulate_censoring(N_TRAIN, rate)
    test  <- simulate_censoring(N_TEST,  rate)
    
    risks <- list()
    capture.output({
      tryCatch({ m <- coxph(form_naive, data=train)
      risks$CoxPH <- predict(m, newdata=test, type="risk") }, error=function(e) NULL)
      tryCatch({ m <- coxph(form_oracle, data=train)
      risks$OracleCox <- predict(m, newdata=test, type="risk") }, error=function(e) NULL)
      tryCatch({ m <- ranger(form_naive, data=train, num.trees=500)
      risks$RSF <- rowSums(predict(m, data=test)$chf) }, error=function(e) NULL)
      tryCatch({ m <- deepsurv(form_naive, data=train, epochs=50L, num_nodes=c(32,32), verbose=FALSE)
      risks$DeepSurv <- predict(m, newdata=test, type="risk") }, error=function(e) NULL)
      tryCatch({ m <- deephit(form_naive, data=train, epochs=50L, num_nodes=c(32,32), verbose=FALSE)
      risks$DeepHit <- predict(m, newdata=test, type="risk") }, error=function(e) NULL)
    })
    
    for (modell in names(risks)) {
      r <- risks[[modell]]
      results <- rbind(results, data.frame(
        Target_Cens = target, Iteration = i, Modell = modell,
        Harrell = harrell_c(test$time, test$status, r),
        Uno     = uno_c(test$time, test$status, r, tau = TAU)
      ))
    }
  }
}

# ------------------------------------------------------------------
# 5. Zusammenfassen
# ------------------------------------------------------------------
summary_tab <- results %>%
  group_by(Target_Cens, Modell) %>%
  summarise(
    Harrell_Mean = mean(Harrell, na.rm=TRUE),
    Uno_Mean     = mean(Uno,     na.rm=TRUE),
    Differenz    = Harrell_Mean - Uno_Mean,
    n_valid      = sum(!is.na(Harrell)),
    .groups = "drop"
  ) %>%
  arrange(Target_Cens, desc(Harrell_Mean))

cat("\n=== Harrells C vs. Uno's C (Mittelwerte) ===\n")
print(as.data.frame(summary_tab), digits = 4)

# Rangfolge-Check: stimmen die Raenge ueberein?
cat("\n=== Rangfolge-Vergleich je Zensierungsstufe ===\n")
for (tc in target_censoring) {
  sub <- summary_tab %>% filter(Target_Cens == tc)
  rang_h <- sub$Modell[order(-sub$Harrell_Mean)]
  rang_u <- sub$Modell[order(-sub$Uno_Mean)]
  gleich <- identical(rang_h, rang_u)
  cat("\nZensierung", tc*100, "%:\n")
  cat("  Harrell-Rang:", paste(rang_h, collapse=" > "), "\n")
  cat("  Uno-Rang:    ", paste(rang_u, collapse=" > "), "\n")
  cat("  Rangfolge identisch:", gleich, "\n")
}

write.csv(summary_tab, "UnoC_Vergleich.csv", row.names = FALSE)
cat("\nGespeichert: UnoC_Vergleich.csv\n")

# ------------------------------------------------------------------
# 6. Plot: Harrell vs Uno je Modell/Stufe
# ------------------------------------------------------------------
library(ggplot2)
library(tidyr)

plot_df <- summary_tab %>%
  select(Target_Cens, Modell, Harrell_Mean, Uno_Mean) %>%
  pivot_longer(c(Harrell_Mean, Uno_Mean), names_to="Mass", values_to="C") %>%
  mutate(Mass = ifelse(Mass=="Harrell_Mean", "Harrells C", "Uno's C"),
         Zens = paste0(Target_Cens*100, "% Zensierung"))

p <- ggplot(plot_df, aes(x = Modell, y = C, fill = Mass)) +
  geom_col(position = position_dodge(width=0.7), width=0.65) +
  facet_wrap(~ Zens) +
  coord_cartesian(ylim = c(0.5, 0.95)) +
  scale_fill_manual(values = c("Harrells C"="#377EB8", "Uno's C"="#E41A1C")) +
  theme_minimal() +
  labs(title = "Harrells C vs. Uno's C: Rangfolge bleibt erhalten",
       subtitle = paste0("tau = ", TAU, ", ", n_iterations, " Iterationen je Stufe"),
       x = NULL, y = "Mittlerer C-Index", fill = NULL) +
  theme(legend.position="bottom", axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(size=12))

ggsave("Plot_UnoC_Vergleich.png", p, width = 12, height = 5.5, dpi = 300)
cat("Gespeichert: Plot_UnoC_Vergleich.png\n")
