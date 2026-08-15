# Survival Analysis: Cox vs. Machine Learning

Begleitcode zur Bachelorarbeit **„Prädiktive Performanz vs. Interpretierbarkeit: Vergleich klassischer und Deep-Learning-Ansätze in der Überlebenszeitanalyse"**.

Dieses Repository enthält den vollständigen Quellcode einer Simulationsstudie, die klassische Cox-Modelle mit maschinellen Lernverfahren (Random Survival Forest, DeepSurv, DeepHit) hinsichtlich Vorhersagegüte und Interpretierbarkeit vergleicht. Sämtliche in der Arbeit berichteten Ergebnisse und Abbildungen lassen sich mit diesem Code reproduzieren.

---

## Überblick

Die Studie basiert auf einem kontrollierten Datengenerierungsprozess, dessen Effektstruktur bewusst die Annahmen des Standard-Cox-Modells verletzt (nichtlineare Effekte, Wechselwirkung, Verletzung der Proportional-Hazards-Annahme). Verglichen werden fünf Modelle:

- **naiv-Cox** – ein linear spezifiziertes Cox-Modell ohne Kenntnis der wahren Struktur
- **Oracle-Cox** – ein Cox-Modell mit korrekter Spezifikation (Splines, quadratischer Term, Interaktionsterm)
- **Random Survival Forest (RSF)**
- **DeepSurv** – neuronales Netz mit Cox-partieller Likelihood
- **DeepHit** – neuronales Netz mit diskreter Zeitbehandlung

Bewertet wird über den C-Index (Diskriminierung) und den Integrierten Brier Score (Kalibrierung), über ein breites Spektrum von Zensierungsgraden und Trainingsgrößen.

---

## Voraussetzungen

### R
Der Code wurde mit **R Version 4.4.3** entwickelt. Folgende Pakete werden benötigt:

```r
install.packages(c(
  "survival",        # Cox-Modell, Schoenfeld-Residuen
  "ranger",          # Random Survival Forest
  "survivalmodels",  # DeepSurv, DeepHit (Schnittstelle zu pycox)
  "survC1",          # Uno's C-Index
  "pec",             # Integrierter Brier Score  [falls verwendet – ggf. anpassen]
  "ggplot2",         # Abbildungen
  "dplyr",           # Datenaufbereitung
  "reshape2"       # Umformen der Ergebnistabellen (melt/dcast) für die Auswertung
))
```

### Python / PyTorch (für die neuronalen Verfahren)
DeepSurv und DeepHit werden über das R-Paket `survivalmodels` angesprochen, das im Hintergrund **Python mit PyTorch und dem Paket `pycox`** benötigt. Die Einrichtung erfolgt über `reticulate`:

```r
install.packages("reticulate")
reticulate::install_miniconda()
```
Anschließend werden PyTorch und pycox eingerichtet:

```r
library(survivalmodels)
install_pycox(pip = TRUE, install_torch = TRUE)
```

Ohne eine funktionierende Python-Umgebung lassen sich nur die klassischen Modelle und der RSF ausführen.

---

## Reproduzierbarkeit

Vor jeder Wiederholung wird der Startwert des Zufallszahlengenerators deterministisch gesetzt, sowohl für R als auch für PyTorch. Dadurch sind alle Datensätze und Modellanpassungen exakt nachbildbar. Ein wiederholter Durchlauf mit identischen Einstellungen liefert dieselben Ergebnisse.

---

## Struktur und Ausführreihenfolge

```
.
├── README.md
├── hauptlauf/
│   ├── hauptlauf_v4_repro.R        # zentraler Performance-Vergleich (alle 5 Modelle, 7 Zensierungsstufen)
│   └── plots_hauptlauf_v4.R        # Abbildungen: C-Index und IBS nach Zensierung
├── skalierung/
│   ├── n_vergleich_single_v5.R     # Skalierungsvergleich über Trainingsgrößen
│   ├── n_vergleich_merge_v5.R      # Zusammenführung + Ereigniszahl-Auswertung (Uno's C) + Plot
├── interpretierbarkeit/
│   ├── forest_plot_hr_v5_untereinander.R # Forest-Plot der Hazard Ratios mit Wahrheit (naiv-Cox)
│   ├── schoenfeld_ph_v3.R          # PH-Diagnostik über Schoenfeld-Residuen
│   ├── forestplotest.R             # HR's rechnen und Plot ohne Wahrheit
│   └── rsf_importance_v3.R         # Permutation Variable Importance (RSF)
├── sensitivität/
│   └── uno_c_sensitivitaet.R       # Sensitivitätsanalyse Harrell vs. Uno
└── synthese/
    └── tradeoff_plot_v5.R          # Trade-off Vorhersagegüte vs. Interpretierbarkeit
```


**Empfohlene Reihenfolge:**

1. **Hauptlauf** ausführen (`hauptlauf_v4_repro.R`) → erzeugt die zentralen Ergebnis-CSVs
2. Dann (`plots_hauptlauf_v4.R`) → erzeugt die Plots
3. **Skalierungsvergleich** ausführen (`n_vergleich_single_v5.R`, dann `n_vergleich_merge_v5.R`)
4. **Interpretierbarkeits-Analysen** ausführen (`forestplotest.R`) → (`forest_plot_hr_v5_untereinander.R`) → (`schoenfeld_ph_v3.R`) → (`rsf_importance_v3.R`) 
5. **Sensitivitätsanalyse** ausführen (`uno_c_sensitivitaet.R`)
6. **Synthese** ausführen (`tradeoff_plot_v5.R`)

---

## Hinweise

- **Pfade:** Alle Skripte verwenden relative Pfade. Das Arbeitsverzeichnis sollte auf den Repository-Ordner gesetzt werden.
- **Rechenzeit:** Der Random Survival Forest wird im Skalierungsvergleich aus Rechenzeitgründen nur bis zu einer Trainingsgröße von n = 1000 einbezogen.
- **DeepHit unter hoher Zensierung:** Bei 90 % Zensierung kann das Training von DeepHit fehlschlagen; solche Ausfälle werden protokolliert und in der Auswertung als fehlende Werte behandelt.

---

## Zitation

Bei Verwendung dieses Codes bitte die zugehörige Bachelorarbeit angeben:

> Fleischer, Kilian (2026). *Prädiktive Performanz vs. Interpretierbarkeit: Vergleich klassischer und Deep-Learning-Ansätze in der Überlebenszeitanalyse.* Bachelorarbeit, Technische Hochschule Augsburg.
