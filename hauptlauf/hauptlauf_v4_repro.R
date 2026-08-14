# ===================================================================
# HAUPTLAUF v4 - Performanz-Vergleich (reproduzierbar + df=6)
# -------------------------------------------------------------------
# Änderungen gegenüber v3:
#   A) Oracle-Spline df = 6 statt df = 4
#      (df=4 approximiert 1.5*sin(pi*x1) über zwei Perioden zu grob;
#       MaxAbs-Fehler 1.78 -> mit df=6 nur noch 0.10)
#   B) NN-Reproduzierbarkeit: survivalmodels::set_seed() wird vor
#      jeder Iteration gesetzt und steuert damit den PyTorch-RNG,
#      den set.seed() allein NICHT erreicht. Dadurch sind auch
#      DeepSurv und DeepHit reproduzierbar.
# ===================================================================

library(survival)
library(ranger)
library(survivalmodels)
library(reticulate)
library(ggplot2)
library(dplyr)
library(reshape2)
library(splines)

use_condaenv("r-reticulate", required = TRUE)

# -------------------------------------------------------------------
# 1. Simulationsfunktion mit Hazard-Faktor 0.3
# -------------------------------------------------------------------
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

# -------------------------------------------------------------------
# 2. Kalibrierung der Zensierungsraten
# -------------------------------------------------------------------
calibrate_rate <- function(target_cens, n_calib = 20000, seed = 1) {
  f <- function(rate) {
    set.seed(seed)
    d <- simulate_censoring(n = n_calib, cens_rate = rate)
    mean(d$status == 0) - target_cens
  }
  tryCatch(
    uniroot(f, interval = c(1e-4, 50), tol = 1e-4)$root,
    error = function(e) {
      warning(sprintf("Ziel %.0f%% evtl. nicht erreichbar: %s", target_cens*100, e$message))
      1e-4
    }
  )
}

# -------------------------------------------------------------------
# 3. Metriken
# -------------------------------------------------------------------
c_idx <- function(time, status, risk_score) {
  if (all(is.na(risk_score))) return(NA_real_)
  concordance(Surv(time, status) ~ risk_score, reverse = TRUE)$concordance
}

get_surv_matrix <- function(model_type, model, test_data, t_grid) {
  if (model_type %in% c("cox", "oracle_cox")) {
    sf <- survfit(model, newdata = test_data)
    surv_mat <- sapply(1:nrow(test_data), function(i) {
      stepfun(sf$time, c(1, sf$surv[, i]))(t_grid)
    })
    return(t(surv_mat))
  }
  if (model_type == "rsf") {
    pred <- predict(model, data = test_data)
    surv_mat <- sapply(1:nrow(test_data), function(i) {
      stepfun(pred$unique.death.times, c(1, pred$survival[i, ]))(t_grid)
    })
    return(t(surv_mat))
  }
  if (model_type %in% c("deepsurv", "deephit")) {
    surv_pred <- predict(model, newdata = test_data, type = "survival")
    pred_times <- as.numeric(colnames(surv_pred))
    surv_mat <- t(sapply(1:nrow(surv_pred), function(i) {
      approx(x = pred_times, y = surv_pred[i, ], xout = t_grid,
             method = "linear", rule = 2)$y
    }))
    return(surv_mat)
  }
  stop("Unbekannter Modelltyp: ", model_type)
}

ibs_ipcw <- function(surv_matrix, test_time, test_status,
                     train_time, train_status, t_grid) {
  if (all(is.na(surv_matrix))) return(NA_real_)
  km_cens <- survfit(Surv(train_time, 1 - train_status) ~ 1)
  G_fun <- stepfun(km_cens$time, c(1, km_cens$surv))
  G_min <- 0.01
  n_test <- length(test_time)
  brier_t <- numeric(length(t_grid))
  for (k in seq_along(t_grid)) {
    t_star <- t_grid[k]
    G_tstar <- max(G_fun(t_star), G_min)
    G_Ti    <- pmax(G_fun(test_time), G_min)
    case1 <- (test_time <= t_star) & (test_status == 1)
    case2 <- (test_time > t_star)
    weights <- numeric(n_test)
    weights[case1] <- 1 / G_Ti[case1]
    weights[case2] <- 1 / G_tstar
    y_true <- numeric(n_test)
    y_true[case2] <- 1
    s_pred <- surv_matrix[, k]
    brier_t[k] <- mean(weights * (s_pred - y_true)^2)
  }
  delta <- diff(t_grid)
  sum(0.5 * (brier_t[-1] + brier_t[-length(brier_t)]) * delta) /
    (max(t_grid) - min(t_grid))
}

# -------------------------------------------------------------------
# 4. Setup
# -------------------------------------------------------------------
target_censoring <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.70, 0.90)
n_iterations     <- 100
t_grid <- seq(0.1, 4.5, length.out = 50)

form_naive  <- Surv(time, status) ~ x1 + x2 + x3 + x4
form_oracle <- Surv(time, status) ~ ns(x1, df = 6)*x3 + I(x2^2)   # <-- df = 6 (war 4)

cat("Kalibriere Zensierungs-Raten (Hazard-Faktor =", HAZARD_FACTOR, ")...\n")
calibrated_rates <- sapply(target_censoring, calibrate_rate)
names(calibrated_rates) <- paste0(target_censoring * 100, "%")
print(round(calibrated_rates, 4))

all_results <- data.frame()

# -------------------------------------------------------------------
# 5. Hauptschleife
# -------------------------------------------------------------------
cat("\nStarte Production Run (", length(target_censoring) * n_iterations,
    "Iterationen x 5 Modelle )...\n")
cat("Reproduzierbar: R-RNG + PyTorch-RNG (survivalmodels::set_seed) synchron.\n")

for (k in seq_along(target_censoring)) {
  target <- target_censoring[k]
  rate   <- calibrated_rates[k]
  
  cat("\n======================================================\n")
  cat("ZIEL-ZENSIERUNG:", target * 100, "% (Rate =", round(rate, 4), ")\n")
  cat("======================================================\n")
  
  for (i in 1:n_iterations) {
    if (i %% 20 == 0) cat("  Iteration", i, "von", n_iterations, "...\n")
    
    iter_seed <- round(target * 1000) + i
    set.seed(iter_seed)                 # steuert R: Simulation, ranger
    set_seed(iter_seed)                 # steuert PyTorch: DeepSurv, DeepHit
    train_data <- simulate_censoring(n = 1000, cens_rate = rate)
    test_data  <- simulate_censoring(n = 1000, cens_rate = rate)
    emp_cens   <- mean(test_data$status == 0)
    
    c_cox <- c_ora <- c_rsf <- c_ds <- c_hit <- NA_real_
    ibs_cox <- ibs_ora <- ibs_rsf <- ibs_ds <- ibs_hit <- NA_real_
    
    capture.output({
      tryCatch({
        m_cox <- coxph(form_naive, data = train_data, x = TRUE)
        c_cox   <- c_idx(test_data$time, test_data$status,
                         predict(m_cox, newdata = test_data, type = "risk"))
        ibs_cox <- ibs_ipcw(get_surv_matrix("cox", m_cox, test_data, t_grid),
                            test_data$time, test_data$status,
                            train_data$time, train_data$status, t_grid)
      }, error = function(e) cat("  Cox-Fehler:", e$message, "\n"))
      
      tryCatch({
        m_ora <- coxph(form_oracle, data = train_data, x = TRUE)
        c_ora   <- c_idx(test_data$time, test_data$status,
                         predict(m_ora, newdata = test_data, type = "risk"))
        ibs_ora <- ibs_ipcw(get_surv_matrix("oracle_cox", m_ora, test_data, t_grid),
                            test_data$time, test_data$status,
                            train_data$time, train_data$status, t_grid)
      }, error = function(e) cat("  OracleCox-Fehler:", e$message, "\n"))
      
      tryCatch({
        m_rsf <- ranger(form_naive, data = train_data, num.trees = 500, importance = "none")
        c_rsf   <- c_idx(test_data$time, test_data$status,
                         rowSums(predict(m_rsf, data = test_data)$chf))
        ibs_rsf <- ibs_ipcw(get_surv_matrix("rsf", m_rsf, test_data, t_grid),
                            test_data$time, test_data$status,
                            train_data$time, train_data$status, t_grid)
      }, error = function(e) cat("  RSF-Fehler:", e$message, "\n"))
      
      tryCatch({
        m_ds <- deepsurv(form_naive, data = train_data, epochs = 50L,
                         num_nodes = c(32, 32), verbose = FALSE)
        c_ds   <- c_idx(test_data$time, test_data$status,
                        predict(m_ds, newdata = test_data, type = "risk"))
        ibs_ds <- ibs_ipcw(get_surv_matrix("deepsurv", m_ds, test_data, t_grid),
                           test_data$time, test_data$status,
                           train_data$time, train_data$status, t_grid)
      }, error = function(e) cat("  DeepSurv-Fehler:", e$message, "\n"))
      
      tryCatch({
        m_hit <- deephit(form_naive, data = train_data, epochs = 50L,
                         num_nodes = c(32, 32), verbose = FALSE)
        c_hit   <- c_idx(test_data$time, test_data$status,
                         predict(m_hit, newdata = test_data, type = "risk"))
        ibs_hit <- ibs_ipcw(get_surv_matrix("deephit", m_hit, test_data, t_grid),
                            test_data$time, test_data$status,
                            train_data$time, train_data$status, t_grid)
      }, error = function(e) cat("  DeepHit-Fehler:", e$message, "\n"))
    })
    
    all_results <- rbind(all_results, data.frame(
      Target_Cens = target, Empirische_Cens = emp_cens, Iteration = i,
      CoxPH_C=c_cox, OracleCox_C=c_ora, RSF_C=c_rsf, DeepSurv_C=c_ds, DeepHit_C=c_hit,
      CoxPH_IBS=ibs_cox, OracleCox_IBS=ibs_ora, RSF_IBS=ibs_rsf, DeepSurv_IBS=ibs_ds, DeepHit_IBS=ibs_hit
    ))
    
    if (i %% 20 == 0) write.csv(all_results, "Ergebnisse_v4_Rohdaten.csv", row.names = FALSE)
  }
}

# -------------------------------------------------------------------
# 6. Aufbereiten
# -------------------------------------------------------------------
cat("\nBerechnung fertig! Erstelle Zusammenfassung...\n")

long_c <- all_results %>%
  select(Target_Cens, Empirische_Cens, Iteration, ends_with("_C")) %>%
  rename_with(~ sub("_C$", "", .x), ends_with("_C")) %>%
  melt(id.vars = c("Target_Cens","Empirische_Cens","Iteration"),
       variable.name = "Modell", value.name = "C_Index")

long_ibs <- all_results %>%
  select(Target_Cens, Empirische_Cens, Iteration, ends_with("_IBS")) %>%
  rename_with(~ sub("_IBS$", "", .x), ends_with("_IBS")) %>%
  melt(id.vars = c("Target_Cens","Empirische_Cens","Iteration"),
       variable.name = "Modell", value.name = "IBS")

summary_c <- long_c %>%
  group_by(Target_Cens, Modell) %>%
  summarise(Mean=mean(C_Index,na.rm=TRUE), SD=sd(C_Index,na.rm=TRUE),
            Mean_Cens=mean(Empirische_Cens,na.rm=TRUE),
            n_valid=sum(!is.na(C_Index)), .groups='drop')

summary_ibs <- long_ibs %>%
  group_by(Target_Cens, Modell) %>%
  summarise(Mean=mean(IBS,na.rm=TRUE), SD=sd(IBS,na.rm=TRUE),
            Mean_Cens=mean(Empirische_Cens,na.rm=TRUE),
            n_valid=sum(!is.na(IBS)), .groups='drop')

cat("\n--- C-Index ---\n"); print(summary_c, n=Inf)
cat("\n--- IBS ---\n"); print(summary_ibs, n=Inf)

# -------------------------------------------------------------------
# 7. Speichern v4-Dateiname
# -------------------------------------------------------------------
write.csv(all_results, "Ergebnisse_v4_Rohdaten.csv", row.names=FALSE)
write.csv(summary_c,   "Ergebnisse_v4_CIndex.csv", row.names=FALSE)
write.csv(summary_ibs, "Ergebnisse_v4_IBS.csv", row.names=FALSE)

cat("\nFertig! Dateien: Ergebnisse_v4_Rohdaten.csv, _CIndex.csv, _IBS.csv\n")
cat("Zum Plotten das bestehende Plot-Skript auf die v4-CSVs zeigen lassen.\n")
