library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
library(ggpubr)
library(ggsignif)


main_theme <- theme(
  plot.title = element_text(size = 30, hjust = 0.5),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 16),
  strip.text = element_text(size = 20),
  legend.title = element_text(size = 16, hjust = 0.5),
  legend.text = element_text(size = 12),
  panel.grid = element_blank()
)

# Summary Figures ---------------------------------------------------------


data <- read.csv("data/calculations.csv", check.names = FALSE) %>%
  mutate(pen_number = as.factor(pen_number))

data_sum <- data %>%
  summarize(
    across(c(MPR, MS, AUC), list(mean = mean, stdev = sd)), 
    count = n(),
    .by = c(pen_number, mpi)
  ) %>%
  pivot_longer(
    cols = matches("[A-Z]{2,3}_[a-z]+", ignore.case = FALSE),
    names_to = c("metric", ".value"),
    names_sep = "_"
  )

results <- read.csv("data/results.csv", check.names = FALSE) %>%
  separate(contrast, c("pen_number", "pen_ctrl"), " - ", remove=FALSE) %>%
  mutate(across(c("pen_number", "pen_ctrl"), ~as.factor(str_remove(., "pen_number"))))

data_joined <- right_join(data_sum, results)


# Metric trends over time ------------------------------------------------


sig_points <- data_joined %>% filter(p.adj < 0.05)

data_norm <- data_sum %>%
  filter(!is.na(pen_number)) %>%
  pivot_wider(names_from=pen_number, values_from=c(mean, stdev)) %>%
  mutate(across(matches("mean_[2-4]"), ~ . - mean_1)) %>%
  pivot_longer(matches("^[a-z]+_\\d$"), names_to = c(".value", "pen_number"), names_sep = "_") %>%
  mutate(pen_number = factor(pen_number, c("1", "2", "3", "4"), c("Pen 1", "Pen 2", "Pen 3", "Pen 4")))

data_norm %>%
  filter(!is.na(pen_number) & pen_number != "Pen 1", !is.na(mean)) %>%
  arrange(mpi) %>%
  group_by(metric) %>%
  mutate(across(c(mean, stdev), ~ scale(., FALSE)), metric = as.factor(metric)) %>%
  ggplot(aes(mpi, mean, color = metric, group = metric)) +
  geom_hline(
    yintercept = 0, linetype="dashed", linewidth = 1, color="black", alpha = 0.5
  ) +
  geom_path(linewidth=1) +
  geom_point(size = 2, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 13, 1)) +
  facet_grid(rows = vars(pen_number)) +
  scale_color_brewer(palette = "Set1", name = "Metric") +
  labs(
    x = "Months Post-Inoculation",
    y = "Normalized Mean Value",
    title = "RT-QuIC Metrics Over Time by Pen",
  ) +
  main_theme

ggsave("figures/metrics_over_time.png", width = 16, height = 12, dpi = 300)

# pen_comparisons <- t(expand.grid(paste("Pen", seq(2, 4)), "Pen 1", stringsAsFactors = FALSE))
pen_comparisons <- list(c("Pen 1", "Pen 2"), c("Pen 1", "Pen 3"), c("Pen 1", "Pen 4"))

data %>%
  filter(!is.na(pen_number), !is.na(mpi)) %>%
  pivot_longer(c(MPR, MS, AUC), names_to = "metric") %>%
  mutate(
    pen_number = paste("Pen", pen_number),
    mpi = factor(mpi, levels = seq(0, 13), labels = paste(seq(0, 13), "MPI"))
  ) %>%
  ggplot(aes(pen_number, value, color=pen_number)) +
  geom_boxplot() +
  # geom_point() +
  facet_grid(
    cols=vars(mpi), 
    rows = vars(metric), 
    scales="free"
  ) +
  stat_compare_means(
    size = 6,
    comparisons = pen_comparisons, 
    label = "p.signif",
    method = "wilcox.test",
    method.args = list(alternative = "less")
  ) +
  main_theme +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

ggsave("figures/comparisons_over_time.png", width = 16, height = 12, dpi = 300)

# Figure 2: Significance heatmap ------------------------------------------


results %>%
  filter(!is.na(p.adj)) %>%
  mutate(metric = factor(metric, levels = c("MPR", "MS", "AUC"))) %>%
  ggplot(aes(mpi, metric, fill = -log10(p.adj), label = paste("p =", signif(p.adj, 3)))) +
  geom_tile(color="black") +
  geom_text(
    size = 5, color = "white", fontface = "bold"
  ) +
  facet_grid(vars(paste0("Pen #", pen_number))) +
  scale_fill_gradient(
    low = "darkcyan", high = "#c0392b",
    limits = c(-log10(1), -log10(0.01)),
    name = expression(-log[10](p[adj]))
  ) +
  scale_x_continuous(breaks = seq(0, 13, 1), expand = FALSE) +
  scale_y_discrete(expand=FALSE) +
  labs(
    x = "Months Post-Inoculation",
    title = "Significance Over Time vs. Negative Control (Pen 1)"
  ) +
  main_theme +
  theme(
    axis.title.y = element_blank()
  )

ggsave("figures/significance_heatmap.png", width = 16, height = 12, dpi = 300)
