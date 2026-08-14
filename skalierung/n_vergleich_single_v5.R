# ===================================================================
# n-VERGLEICH - EINZELLAUF v5 (Harrell UND Uno's C)
# -------------------------------------------------------------------
# Wie v4 (df=6, NN-Seed), aber schreibt pro Modell ZWEI Maße:
#   - Harrell C  (Spalten *_H)  -> für Facetten-Plot (innerhalb Stufe)
#   - Uno's C    (Spalten *_U)  -> für Ereigniszahl-Plot (über Stufen)
#
# Uno's C über survC1::Est.Cval, fester Horizont tau = 5 (Studienende).
# Damit wird der Ereigniszahl-Plot über Zensierungsstufen hinweg
# vergleichbar und ist nicht mehr durch die Harrell-Verzerrung belastet.
#
# survC1 nötig:  install.packages("survC1")
#
# ANLEITUNG - Skript FÜNFMAL laufen lassen, Konfig-Block anpassen:
#   Lauf 1:  N_THIS_RUN <- 100     INCLUDE_RSF <- TRUE
#   Lauf 2:  N_THIS_RUN <- 300     INCLUDE_RSF <- TRUE
#   Lauf 3:  N_THIS_RUN <- 1000    INCLUDE_RSF <- TRUE
#   Lauf 4:  N_THIS_RUN <- 3000    INCLUDE_RSF <- FALSE
#   Lauf 5:  N_THIS_RUN <- 10000   INCLUDE_RSF <- FALSE
#
# Danach: n_vergleich_merge_v5.R ausführen.
# ===================================================================

# ===================== KONFIGURATION ===============================
N_THIS_RUN   <- 10000   # 100 / 300 / 1000 / 3000 / 10000
INCLUDE_RSF  <- FALSE    # FALSE ab n = 3000
RSF_TREES    <- 500
RSF_MIN_NODE <- 3
n_iterations <- 50
N_TEST       <- 5000
TAU          <- 5       # Auswertungshorizont für Uno's C
# ===================================================================

library(survival)
library(ranger)
library(survivalmodels)
library(reticulate)
library(dplyr)
library(splines)

if (!requireNamespace("survC1", quietly = TRUE)) {
  stop("Paket 'survC1' fehlt. Bitte installieren: install.packages('survC1')")
}

use_condaenv("r-reticulate", required = TRUE)

# -------------------------------------------------------------------
# 1. Simulation (identisch zu v4)
# -------------------------------------------------------------------
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
           error = function(e) { warning("Kalibrierung fehlgeschlagen"); 1e-4 })
}

# -------------------------------------------------------------------
# Metriken: Harrell UND Uno
# -------------------------------------------------------------------
harrell_c <- function(time, status, risk) {
  if (all(is.na(risk)) || length(unique(risk)) < 2) return(NA_real_)
  concordance(Surv(time, status) ~ risk, reverse = TRUE)$concordance
}

uno_c <- function(time, status, risk, tau = 5) {
  if (all(is.na(risk)) || length(unique(risk)) < 2) return(NA_real_)
  df <- data.frame(time = time, status = status, marker = risk)
  tryCatch(
    survC1::Est.Cval(mydata = df, tau = tau, nofit = TRUE)$Dhat,
    error = function(e) NA_real_
  )
}

# -------------------------------------------------------------------
# 2. Setup
# -------------------------------------------------------------------
target_censoring <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.70, 0.90)

form_naive  <- Surv(time, status) ~ x1 + x2 + x3 + x4
form_oracle <- Surv(time, status) ~ ns(x1, df = 6)*x3 + I(x2^2)

outfile <- paste0("nVergleich_v5_n", N_THIS_RUN,
                  if (!INCLUDE_RSF) "_ohneRSF" else "", ".csv")

cat("========================================\n")
cat(" n-Vergleich v5 | n =", N_THIS_RUN, "\n")
cat(" RSF:", if (INCLUDE_RSF)
  paste0("an (", RSF_TREES, " Baeume, min.node.size = ", RSF_MIN_NODE, ")")
  else "AUS", "\n")
cat(" Iterationen:", n_iterations, "| Testset:", N_TEST, "| tau =", TAU, "\n")
cat(" Masse: Harrell (_H) UND Uno (_U) | Oracle df=6 | NN-Seed gesetzt\n")
cat(" Ausgabe:", outfile, "\n")
cat("========================================\n\n")

cat("Kalibriere Zensierungs-Raten...\n")
rates <- sapply(target_censoring, calibrate_rate)
names(rates) <- paste0(target_censoring * 100, "%")
print(round(rates, 4))

results <- data.frame()
t_start <- Sys.time()

# -------------------------------------------------------------------
# 3. Hauptschleife
# -------------------------------------------------------------------
for (k in seq_along(target_censoring)) {
  target <- target_censoring[k]
  rate   <- rates[k]

  cat("\n-------------------------------------------\n")
  cat("n =", N_THIS_RUN, "| Zensierung:", target * 100, "% |",
      format(Sys.time(), "%H:%M:%S"), "\n")
  cat("-------------------------------------------\n")

  for (i in 1:n_iterations) {
    cat("  Iter", i, "/", n_iterations, "-", format(Sys.time(), "%H:%M:%S"), "\n")

    iter_seed <- N_THIS_RUN * 10000 + round(target * 1000) + i
    set.seed(iter_seed)
    set_seed(iter_seed)
    train <- simulate_censoring(n = N_THIS_RUN, cens_rate = rate)
    test  <- simulate_censoring(n = N_TEST,     cens_rate = rate)

    n_events <- sum(train$status)

    # Risk-Scores je Modell zwischenspeichern (NULL, falls Modell scheitert)
    risk <- list(CoxPH = NULL, OracleCox = NULL, RSF = NULL,
                 DeepSurv = NULL, DeepHit = NULL)

    capture.output({
      tryCatch({
        m <- coxph(form_naive, data = train)
        risk$CoxPH <- predict(m, newdata = test, type = "risk")
      }, error = function(e) invisible(NULL))

      tryCatch({
        m <- coxph(form_oracle, data = train)
        risk$OracleCox <- predict(m, newdata = test, type = "risk")
      }, error = function(e) invisible(NULL))

      if (INCLUDE_RSF) {
        tryCatch({
          m <- ranger(form_naive, data = train,
                      num.trees = RSF_TREES,
                      min.node.size = RSF_MIN_NODE,
                      importance = "none")
          p <- predict(m, data = test)
          risk$RSF <- rowSums(p$chf)
          rm(m, p); gc(verbose = FALSE)
        }, error = function(e) invisible(NULL))
      }

      tryCatch({
        m <- deepsurv(form_naive, data = train, epochs = 50L,
                      num_nodes = c(32, 32), verbose = FALSE)
        risk$DeepSurv <- predict(m, newdata = test, type = "risk")
      }, error = function(e) invisible(NULL))

      tryCatch({
        m <- deephit(form_naive, data = train, epochs = 50L,
                     num_nodes = c(32, 32), verbose = FALSE)
        risk$DeepHit <- predict(m, newdata = test, type = "risk")
      }, error = function(e) invisible(NULL))
    })

    # Beide Masse aus den gespeicherten Risk-Scores
    get_H <- function(r) if (is.null(r)) NA_real_ else harrell_c(test$time, test$status, r)
    get_U <- function(r) if (is.null(r)) NA_real_ else uno_c(test$time, test$status, r, tau = TAU)

    results <- rbind(results, data.frame(
      n_train     = N_THIS_RUN,
      Target_Cens = target,
      Iteration   = i,
      n_events    = n_events,
      CoxPH_H     = get_H(risk$CoxPH),     CoxPH_U     = get_U(risk$CoxPH),
      OracleCox_H = get_H(risk$OracleCox), OracleCox_U = get_U(risk$OracleCox),
      RSF_H       = get_H(risk$RSF),       RSF_U       = get_U(risk$RSF),
      DeepSurv_H  = get_H(risk$DeepSurv),  DeepSurv_U  = get_U(risk$DeepSurv),
      DeepHit_H   = get_H(risk$DeepHit),   DeepHit_U   = get_U(risk$DeepHit)
    ))
  }

  write.csv(results, outfile, row.names = FALSE)
  cat("  -> Zwischenstand gespeichert.\n")
}

write.csv(results, outfile, row.names = FALSE)

cat("\n========================================\n")
cat(" FERTIG fuer n =", N_THIS_RUN, "\n")
cat(" Dauer:", round(difftime(Sys.time(), t_start, units = "mins"), 1), "Minuten\n")
cat(" Datei:", outfile, "\n")
cat("========================================\n")
cat("\nNaechster Schritt: Konfiguration oben anpassen und erneut laufen lassen.\n")
