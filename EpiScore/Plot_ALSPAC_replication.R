library(tidyverse)

setwd("C:/Users/s1654019/Desktop/_old")
ALSPAC = read_csv("Replication_in_ALSPAC.csv")

# Custom colors
colours = c("#f3b473", "#c94905", "#218f60", "#155da2", "#4297d8")

# Facet labels
new_labels <- c(
  "current_former" = "Current vs Former",
  "current_never" = "Current vs Never",
  "former_never" = "Former vs Never"
)

# Factor levels for sorting
ALSPAC$Comparison <- factor(ALSPAC$Comparison, levels = c("current_never", "current_former", "former_never"))
ALSPAC$`Time code` <- factor(ALSPAC$`Time code`, levels = c("antenatal", "FOM_FOF", "F17", "F24"))

# Plotting
ggplot(ALSPAC, aes(fill = Score, x = `Time code`, y = Estimate)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(
    ~ Comparison, 
    scales = "free_y", 
    labeller = labeller(Comparison = new_labels)
  ) +
  scale_fill_manual(
    values = colours, 
    labels = c(
      "EpiSmokEr" = "EpiSmokEr",
      "McCartney.et.al." = expression("McCartney " * italic("et al.")), 
      "mCigarette" = "mCigarette",
      "Trejo.et.al." = expression("Trejo " * italic("et al."))
    )
  ) +
  labs(
    x = "ALSPAC Collection Name",
    y = "AUC",
    fill = "Score"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",  # Place legend on top
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_continuous(limits = c(0.5, 1), expand = c(0, 0), oob = scales::squish)
