# ===================================================================
# HAUPTLAUF v4 - FINALE PLOTS (C-Index & IBS)
# -------------------------------------------------------------------
# Nutzt die v4-CSVs (reproduzierbar, df=6). Kein Neurechnen.
# ===================================================================

library(ggplot2)
library(dplyr)

# -------------------------------------------------------------------
# 1. Daten laden  (v4!)
# -------------------------------------------------------------------
summary_c   <- read.csv("Ergebnisse_v4_CIndex.csv")
summary_ibs <- read.csv("Ergebnisse_v4_IBS.csv")

# Modell-Reihenfolge und feste Farben
modell_levels <- c("CoxPH", "OracleCox", "RSF", "DeepSurv", "DeepHit")
farben <- c("CoxPH"="#E41A1C", "OracleCox"="#377EB8", "RSF"="#4DAF4A",
            "DeepSurv"="#984EA3", "DeepHit"="#FF7F00")

summary_c$Modell   <- factor(summary_c$Modell,   levels = modell_levels)
summary_ibs$Modell <- factor(summary_ibs$Modell, levels = modell_levels)

# -------------------------------------------------------------------
# 2. Ausfallstellen identifizieren (n_valid == 0)
# -------------------------------------------------------------------
ausfall_c   <- summary_c   %>% filter(n_valid == 0)
ausfall_ibs <- summary_ibs %>% filter(n_valid == 0)

plot_c_data   <- summary_c   %>% filter(n_valid > 0)
plot_ibs_data <- summary_ibs %>% filter(n_valid > 0)

y_marker_c   <- min(plot_c_data$Mean   - plot_c_data$SD,   na.rm = TRUE)
y_marker_ibs <- max(plot_ibs_data$Mean + plot_ibs_data$SD, na.rm = TRUE)

# -------------------------------------------------------------------
# 3. C-INDEX PLOT
# -------------------------------------------------------------------
plot_c <- ggplot(plot_c_data,
                 aes(x = Mean_Cens, y = Mean, color = Modell, group = Modell)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.02, alpha = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey50")


plot_c <- plot_c +
  scale_color_manual(values = farben, drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = c(0.1,0.2,0.3,0.4,0.5,0.7,0.9)) +
  theme_minimal() +
  labs(title = "Diskriminierung: Harrells C-Index",
       subtitle = "n = 1000, 100 Iterationen, Hazard-Faktor 0.3",
       x = "Empirischer Anteil zensierter Patienten",
       y = "Mittlerer C-Index",
       color = "Modell") +
  theme(legend.position = "bottom", text = element_text(size = 13))

ggsave("Plot_v4_CIndex_final.png", plot = plot_c, width = 11, height = 7, dpi = 300)
print(plot_c)

# -------------------------------------------------------------------
# 4. IBS PLOT
# -------------------------------------------------------------------
plot_ibs <- ggplot(plot_ibs_data,
                   aes(x = Mean_Cens, y = Mean, color = Modell, group = Modell)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.02, alpha = 0.5)


plot_ibs <- plot_ibs +
  scale_color_manual(values = farben, drop = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = c(0.1,0.2,0.3,0.4,0.5,0.7,0.9)) +
  theme_minimal() +
  labs(title = "Kalibrierung + Diskriminierung: Integrierter Brier Score",
       subtitle = "n = 1000, 100 Iterationen, t in [0.1, 4.5]",
       x = "Empirischer Anteil zensierter Patienten",
       y = "Mittlerer IBS",
       color = "Modell") +
  theme(legend.position = "bottom", text = element_text(size = 13))

ggsave("Plot_v4_IBS_final.png", plot = plot_ibs, width = 11, height = 7, dpi = 300)
print(plot_ibs)

cat("\nFertig!\n")
cat("  Plot_v4_CIndex_final.png\n")
cat("  Plot_v4_IBS_final.png\n")
if (nrow(ausfall_c) > 0 || nrow(ausfall_ibs) > 0) {
  cat("\nGekennzeichnete Ausfälle:\n")
  if (nrow(ausfall_c) > 0)
    cat("  C-Index:", paste(ausfall_c$Modell, "bei", ausfall_c$Mean_Cens*100, "%"), "\n")
  if (nrow(ausfall_ibs) > 0)
    cat("  IBS:    ", paste(ausfall_ibs$Modell, "bei", ausfall_ibs$Mean_Cens*100, "%"), "\n")
}
