# ===================================================================
# n-VERGLEICH MERGE v5 - Harrell + Uno
# -------------------------------------------------------------------
# Liest nVergleich_v5_n*.csv (Spalten *_H und *_U).
# Erzeugt:
#   Plot_v5_nVergleich_facet_Harrell.png   (Facetten nach Zensierung, Harrell)
#   Plot_v5_nVergleich_Events_Uno.png      (Ereigniszahl-Plot, UNO'S C)
#   nVergleich_v5_Zusammenfassung.csv       (beide Masse)
# ===================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# -------------------------------------------------------------------
# 1. Dateien einlesen
# -------------------------------------------------------------------
dateien <- list.files(pattern = "^nVergleich_v5_n[0-9]+(_ohneRSF)?\\.csv$")
if (length(dateien) == 0) stop("Keine nVergleich_v5_n*.csv gefunden.")

n_aus_name <- as.numeric(sub("^nVergleich_v5_n([0-9]+).*$", "\\1", dateien))
dateien    <- dateien[order(n_aus_name)]
cat("Gefundene Dateien:\n"); print(dateien)

roh <- do.call(rbind, lapply(dateien, read.csv))

modell_levels <- c("CoxPH", "OracleCox", "RSF", "DeepSurv", "DeepHit")
farben <- c("CoxPH"="#E41A1C", "OracleCox"="#377EB8", "RSF"="#4DAF4A",
            "DeepSurv"="#984EA3", "DeepHit"="#FF7F00")

# -------------------------------------------------------------------
# 2. In Langformat bringen - getrennt nach Mass
# -------------------------------------------------------------------
# Harrell-Spalten
long_H <- roh %>%
  select(n_train, Target_Cens, Iteration, n_events, ends_with("_H")) %>%
  pivot_longer(ends_with("_H"), names_to = "Modell", values_to = "C") %>%
  mutate(Modell = sub("_H$", "", Modell), Mass = "Harrell")

# Uno-Spalten
long_U <- roh %>%
  select(n_train, Target_Cens, Iteration, n_events, ends_with("_U")) %>%
  pivot_longer(ends_with("_U"), names_to = "Modell", values_to = "C") %>%
  mutate(Modell = sub("_U$", "", Modell), Mass = "Uno")

long <- bind_rows(long_H, long_U)
long$Modell <- factor(long$Modell, levels = modell_levels)

# -------------------------------------------------------------------
# 3. Zusammenfassung (Mittel je n / Zensierung / Modell / Mass)
# -------------------------------------------------------------------
summ <- long %>%
  group_by(n_train, Target_Cens, Modell, Mass) %>%
  summarise(C_Mean = mean(C, na.rm = TRUE),
            C_SD   = sd(C, na.rm = TRUE),
            Events = mean(n_events, na.rm = TRUE),
            n_valid = sum(!is.na(C)),
            .groups = "drop")

write.csv(summ, "nVergleich_v5_Zusammenfassung.csv", row.names = FALSE)
cat("Gespeichert: nVergleich_v5_Zusammenfassung.csv\n")

# -------------------------------------------------------------------
# 4. FACETTEN-PLOT (nach Zensierung) - mit HARRELL
#    (Vergleich innerhalb einer Stufe -> Harrell korrekt)
# -------------------------------------------------------------------
facet_dat <- summ %>% filter(Mass == "Harrell", n_valid > 0)

p_facet <- ggplot(facet_dat,
                  aes(x = n_train, y = C_Mean, color = Modell, group = Modell)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  facet_wrap(~ paste0(Target_Cens*100, "% Zensierung")) +
  scale_x_log10() +
  scale_color_manual(values = farben, drop = FALSE) +
  theme_minimal() +
  labs(title = "C-Index nach Trainingsgroesse (Harrell)",
       subtitle = "Vergleich innerhalb jeder Zensierungsstufe",
       x = "Trainingsgroesse n (log)", y = "Harrell C-Index",
       color = "Modell",
       caption = "RSF nur bis n = 1000 berechnet.") +
  theme(legend.position = "bottom", text = element_text(size = 11))

ggsave("Plot_v5_nVergleich_facet_Harrell.png", p_facet,
       width = 12, height = 7, dpi = 300)
cat("Gespeichert: Plot_v5_nVergleich_facet_Harrell.png\n")

# -------------------------------------------------------------------
# 5. EREIGNISZAHL-PLOT - mit UNO'S C
#    (Vergleich ueber Zensierungsstufen -> Uno noetig)
# -------------------------------------------------------------------
ev_dat <- summ %>% filter(Mass == "Uno", n_valid > 0)

p_ev <- ggplot(ev_dat,
               aes(x = Events, y = C_Mean, color = Modell)) +
  geom_point(aes(shape = factor(Target_Cens)), size = 2.4, alpha = 0.8) +
  geom_smooth(aes(group = Modell), method = "loess", se = FALSE,
              linewidth = 1, span = 1) +
  scale_x_log10() +
  scale_color_manual(values = farben, drop = FALSE) +
  scale_shape_discrete(name = "Zensierung") +
  theme_minimal() +
  labs(title = "C-Index nach mittlerer Ereigniszahl im Training (Uno's C)",
       subtitle = "Zensierungsrobustes Mass, tau = 5 - über Stufen vergleichbar",
       x = "Mittlere Ereigniszahl im Training (log)",
       y = "Uno's C-Index",
       color = "Modell") +
  theme(legend.position = "bottom", text = element_text(size = 11))

ggsave("Plot_v5_nVergleich_Events_Uno.png", p_ev,
       width = 11, height = 7, dpi = 300)
cat("Gespeichert: Plot_v5_nVergleich_Events_Uno.png\n")

# -------------------------------------------------------------------
# 6. NA-Quote pro Zelle berichten (fuer DeepHit-Ausfaelle)
# -------------------------------------------------------------------
na_report <- long %>%
  filter(Mass == "Harrell") %>%
  group_by(n_train, Target_Cens, Modell) %>%
  summarise(n_total = n(), n_na = sum(is.na(C)),
            na_quote = round(n_na / n_total, 2), .groups = "drop") %>%
  filter(n_na > 0) %>%
  arrange(desc(na_quote))

if (nrow(na_report) > 0) {
  cat("\n=== Zellen mit Ausfaellen (NA-Quote) ===\n")
  print(as.data.frame(na_report))
  write.csv(na_report, "nVergleich_v5_NA_Report.csv", row.names = FALSE)
} else {
  cat("\nKeine Ausfaelle.\n")
}

cat("\nFertig.\n")

# ===================================================================
# FACETTEN-PLOT (Harrell) MIT SD-BAENDERN - 2 SPALTEN, GROSS
# x-Achse (100/1000/10000) unter JEDEM Panel
# -------------------------------------------------------------------
# Erzeugt: Plot_v5_facet_Baender_2spalten.png
# ===================================================================

library(ggplot2)
library(dplyr)

zus <- read.csv("nVergleich_v5_Zusammenfassung.csv")

dat <- zus %>%
  filter(Mass == "Harrell", n_valid > 0, !is.na(C_Mean))

modell_levels <- c("CoxPH", "OracleCox", "RSF", "DeepSurv", "DeepHit")
dat$Modell <- factor(dat$Modell, levels = modell_levels)

dat$Facette <- factor(paste0(dat$Target_Cens*100, "% Zensierung"),
                      levels = paste0(c(10,20,30,40,50,70,90), "% Zensierung"))

farben <- c("CoxPH"="#E41A1C", "OracleCox"="#377EB8", "RSF"="#4DAF4A",
            "DeepSurv"="#984EA3", "DeepHit"="#FF7F00")

dat <- dat %>%
  mutate(ymin = C_Mean - C_SD, ymax = C_Mean + C_SD)

# feste Achsen-Ticks, damit in jedem Panel dieselben Werte stehen
x_breaks <- c(100, 300, 1000, 3000, 10000)

p <- ggplot(dat, aes(x = n_train, y = C_Mean,
                     color = Modell, fill = Modell, group = Modell)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  facet_wrap(~ Facette, ncol = 2, scales = "free_x") +   # <-- free_x: Achse unter jedem Panel
  scale_x_log10(breaks = x_breaks, labels = x_breaks) +  # <-- feste, identische Ticks
  scale_color_manual(values = farben, drop = FALSE) +
  scale_fill_manual(values = farben, drop = FALSE) +
  theme_minimal(base_size = 15) +
  labs(
    title = "C-Index nach Trainingsgröße (Harrell)",
    subtitle = "Vergleich innerhalb jeder Zensierungsstufe; Band = Standardabweichung",
    x = "Trainingsgröße n (log)", y = "Harrell C-Index",
    color = "Modell", fill = "Modell"
  ) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 11),
        axis.text.x = element_text(size = 10),
        panel.spacing = unit(1.4, "lines"))

ggsave("Plot_v5_facet_Bänder_2spalten_DE_100_1000_10000.png", p,
       width = 10, height = 15, dpi = 300)
cat("Gespeichert: Plot_v5_facet_Baender_2spalten.png\n")
cat("x-Achse (100/300/1000/3000/10000) steht jetzt unter JEDEM Panel.\n")