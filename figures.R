library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)


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
  group_by(pen_number, mpi) %>%
  summarize(
    across(
      c(MPR, MS, AUC), 
      mean
    ), 
    .groups = "drop_last"
  ) %>%
  ungroup() %>%
  pivot_longer(c(MPR, MS, AUC), names_to = "metric")

results <- read.csv("data/results.csv", check.names = FALSE) %>%
  separate(contrast, c("pen_number", "pen_ctrl"), " - ", remove=FALSE) %>%
  mutate(across(c("pen_number", "pen_ctrl"), ~as.factor(str_remove(., "pen_number"))))

data_joined <- right_join(data_sum, results)


# Figure 1: Metric trends over time ---------------------------------------


sig_points <- data_joined %>% filter(p.adj < 0.05)

data_norm <- data_sum %>%
  filter(!is.na(pen_number)) %>%
  pivot_wider(names_from=pen_number, values_from=value) %>%
  mutate(across(c("2", "3", "4"), ~ . - `1`)) %>%
  pivot_longer(c("1", "2", "3", "4"), names_to = "pen_number")

data_norm %>%
  filter(!is.na(pen_number) & pen_number != "1") %>%
  group_by(metric) %>%
  mutate(value = scale(value)) %>%
  ggplot(aes(mpi, value, color = metric, fill = metric, group = metric)) +
  geom_line(
    data = NULL,
    aes(y=0), linetype="dashed", linewidth = 1.5, color="black", alpha = 0.5
  ) +
  geom_line(linewidth = 1.5, alpha = 0.5) +
  geom_point(size = 3, alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 13, 1)) +
  facet_wrap(vars(pen_number), nrow=1) +
  scale_color_brewer(palette = "Set1", name = "Metric") +
  labs(
    x = "Months Post-Inoculation",
    y = "Normalized Mean Value",
    title = "RT-QuIC Metrics Over Time by Pen",
  ) +
  main_theme

ggsave("figures/metrics_over_time.png", width = 16, height = 12, dpi = 300)


# Figure 2: Significance heatmap ------------------------------------------


results %>%
  filter(!is.na(p.adj)) %>%
  mutate(metric = factor(metric, levels = c("MPR", "MS", "AUC"))) %>%
  ggplot(aes(mpi, pen_number, fill = -log10(p.adj), label = paste("p =", signif(p.adj, 3)))) +
  geom_tile(color="black") +
  geom_text(
    size = 5, color = "white", fontface = "bold"
  ) +
  facet_grid(vars(metric)) +
  scale_fill_gradient(
    low = "darkcyan", high = "#c0392b",
    limits = c(-log10(1), -log10(0.01)),
    name = expression(-log[10](p[adj]))
  ) +
  scale_x_continuous(breaks = seq(0, 13, 1), expand = FALSE) +
  scale_y_discrete(expand=FALSE) +
  labs(
    x = "Months Post-Inoculation",
    y = "Pen",
    title = "Significance vs. Negative Control (Pen 1)"
  ) +
  main_theme

ggsave("figures/significance_heatmap.png", width = 16, height = 12, dpi = 300)
