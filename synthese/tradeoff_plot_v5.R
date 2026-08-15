# ===================================================================
# TRADE-OFF-PLOT: Vorhersageguete vs. Interpretierbarkeit (v5)
# -------------------------------------------------------------------
# y-Achse: C-Index, gemittelt über alle 7 Zensierungsstufen (v4-Daten)
# x-Achse: Interpretierbarkeit (dreistufig, aus Abschnitt 3.4.4)
#
# Liest Ergebnisse_v4_CIndex.csv, mittelt pro Modell über die Stufen.
# Caption entfernt (kommt als Word-Bildunterschrift).
#
# Erzeugt: Plot_v5_Tradeoff.png
# ===================================================================

library(ggplot2)
library(dplyr)

ci <- read.csv("Ergebnisse_v4_CIndex.csv")

# Pro Modell über alle gültigen Zensierungsstufen mitteln
agg <- ci %>%
  filter(!is.na(Mean), n_valid > 0) %>%
  group_by(Modell) %>%
  summarise(C_Mean = mean(Mean), .groups = "drop")

# Interpretierbarkeits-Einstufung (aus 3.4.4): hoch=3, mittel=2, niedrig=1
interp <- data.frame(
  Modell = c("CoxPH", "OracleCox", "RSF", "DeepSurv", "DeepHit"),
  Interp = c(3, 2, 1, 1, 1),
  Interp_Label = c("hoch", "mittel", "niedrig", "niedrig", "niedrig")
)

dat <- merge(agg, interp, by = "Modell")

# schöne Anzeigenamen
disp <- c(CoxPH="naiv-Cox", OracleCox="Oracle-Cox", RSF="RSF",
          DeepSurv="DeepSurv", DeepHit="DeepHit")
dat$Label <- disp[dat$Modell]

farben <- c("CoxPH"="#E41A1C", "OracleCox"="#377EB8", "RSF"="#4DAF4A",
            "DeepSurv"="#984EA3", "DeepHit"="#FF7F00")

# leichtes Jitter auf x, damit die drei "niedrig"-Punkte nicht übereinander liegen
set.seed(1)
dat$Interp_jit <- dat$Interp + ifelse(dat$Interp==1, c(-0.12,0,0.12)[rank(dat$C_Mean[dat$Interp==1])], 0)

p <- ggplot(dat, aes(x = Interp_jit, y = C_Mean, color = Modell)) +
  geom_point(size = 6) +
  geom_text(aes(label = Label), vjust = -1.4, size = 5, fontface = "bold",
            show.legend = FALSE) +
  scale_x_continuous(breaks = c(1,2,3),
                     labels = c("niedrig", "mittel", "hoch"),
                     limits = c(0.5, 3.5)) +
  scale_y_continuous(limits = c(0.55, 0.90)) +
  scale_color_manual(values = farben, guide = "none") +
  theme_minimal(base_size = 15) +
  labs(
    title = "Trade-off zwischen Vorhersageguete und Interpretierbarkeit",
    subtitle = "C-Index gemittelt über alle sieben Zensierungsstufen",
    x = "Interpretierbarkeit", y = "C-Index (Vorhersagegüte)"
  ) +
  theme(text = element_text(size = 14),
        panel.grid.minor = element_blank())

ggsave("Plot_v5_Tradeoff.png", p, width = 10, height = 7, dpi = 300)
cat("Gespeichert: Plot_v5_Tradeoff.png\n\n")
cat("Werte:\n")
print(dat[order(-dat$C_Mean), c("Label","C_Mean","Interp_Label")])
