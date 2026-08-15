# ===================================================================
# HINWEIS: Zuerst das Skript "forestplotest.R" ausführen
# ===================================================================
# FOREST-PLOT der Hazard Ratios (naives Cox) - UNTEREINANDER
# mit echten mathematischen Formelzeichen (pi, Indizes, Hochzahlen)
# -------------------------------------------------------------------
# Erzeugt: Plot_ForestHR_v5_untereinander.png
# ===================================================================

library(ggplot2)
library(dplyr)

df <- read.csv("Ergebnisse_v3_ForestHR.csv")
df <- df %>% filter(grepl("Zensierung", Zensierung), Zensierung != "Zensierung")

var_labels <- c("x1" = "x[1]", "x2" = "x[2]", "x3" = "x[3]", "x4" = "x[4]")
df$Variable <- factor(df$Variable, levels = c("x4","x3","x2","x1"))
stufen <- unique(df$Zensierung)

wahrheit <- do.call(rbind, lapply(stufen, function(s) {
  data.frame(
    Zensierung = s,
    Variable   = factor(c("x1","x2","x4"), levels = levels(df$Variable)),
    HR_wahr    = c(exp(1.5), exp(1.2), 1),
    Skala = c("Skala: ~ sin(pi * x[1])",
              "Skala: ~ x[2]^2",
              "Skala: ~ x[4]~(Original)"),
    stringsAsFactors = FALSE
  )
}))

farbe_hr   <- "#2C3E50"
farbe_wahr <- "#C0392B"

p <- ggplot(df, aes(x = HR, y = Variable)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbarh(aes(xmin = HR_low, xmax = HR_high), height = 0.20,
                 color = farbe_hr, linewidth = 0.9) +
  # Punkt: Typ über aes gemappt -> erzeugt Legendeneintrag
  geom_point(aes(color = "Geschätztes HR", shape = "Geschätztes HR"),
             size = 3.4) +
  # Raute: Typ über aes gemappt -> erzeugt Legendeneintrag
  geom_point(data = wahrheit, aes(x = HR_wahr, y = Variable,
                                  color = "Wahre Effektstärke", shape = "Wahre Effektstärke"),
             size = 5.5, inherit.aes = FALSE) +
  geom_text(data = wahrheit, aes(x = HR_wahr, y = Variable, label = Skala),
            color = farbe_wahr, size = 3.8, vjust = -1.2,
            parse = TRUE, inherit.aes = FALSE) +
  facet_wrap(~ Zensierung, ncol = 1) +
  scale_x_log10() +
  scale_y_discrete(labels = function(v) parse(text = var_labels[v])) +
  # Manuelle Zuweisung von Farbe und Form je Typ
  scale_color_manual(
    name = NULL,
    values = c("Geschätztes HR" = farbe_hr, "Wahre Effektstärke" = farbe_wahr)
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Geschätztes HR" = 16, "Wahre Effektstärke" = 18)
  ) +
  theme_minimal(base_size = 15) +
  labs(
    title = "Hazard Ratios des naiven Cox-Modells nach Zensierungsstufe",
    subtitle = "Punkt = geschätztes HR (pro Einheit Originalvariable), Balken = 95%-KI",
    x = "Hazard Ratio (log-Skala)", y = NULL
  ) +
  theme(text = element_text(size = 14),
        strip.text = element_text(face = "bold", size = 15),
        axis.text = element_text(size = 13),
        panel.spacing = unit(1.5, "lines"),
        legend.position = "bottom")

print(p)
ggsave("Plot_ForestHR_v5_untereinander_Legende.png", p,
       width = 11, height = 10, dpi = 300)
cat("Gespeichert: Plot_ForestHR_v5_untereinander.png\n")
