library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)

mean_coverages = read.csv('/Users/shirin/Documents/Edinburgh/Lab/TWIST/Descriptives/Coverage.csv')
new_vs_old = read.csv('/Users/shirin/Documents/Edinburgh/Lab/TWIST/Descriptives/16new_vs_old.csv')
new_vs_old = new_vs_old %>% rename("Approach" = "X")
new_vs_old_single_run = read.csv('/Users/shirin/Documents/Edinburgh/Lab/TWIST/Descriptives/16new_vs_old_single_run.csv')

# Plot
new_vs_old %>%
  ggplot( aes(x=Id, y=Coverage, group=Approach, color=Approach)) +
  geom_line() +
  geom_point() +
  ggtitle("Coverage across samples") +
  ylab("Coverage (X)") +
  xlab("Sample")

new_vs_old$Coverage = ifelse(new_vs_old$Approach == "Old", new_vs_old$Coverage / 2, new_vs_old$Coverage)

new_vs_old %>%
  ggplot( aes(x=Id, y=Coverage, group=Approach, color=Approach)) +
  geom_line() +
  geom_point() +
  ggtitle("Coverage across samples") +
  ylab("Coverage (X)") +
  xlab("Sample")


# Plot
new_vs_old_single_run %>%
  ggplot( aes(x=Id, y=Coverage, group=Method, color=Method)) +
  geom_line() +
  geom_point() +
  ggtitle("Coverage across samples") +
  ylab("Coverage (X)") +
  xlab("Sample")


