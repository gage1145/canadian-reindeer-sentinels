library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)
library(ggpubr)
library(ggsignif)
library(broom)


main_theme <- theme(
	plot.title = element_text(size=30),
	axis.title = element_text(size=24),
	axis.text = element_text(size=20),
	legend.title = element_text(size=24),
	legend.text = element_text(size=20),
	strip.text = element_text(size=24)
)


data <- read.csv("data/calculations.csv", check.names = FALSE) %>%
  mutate(pen_number = as.factor(pen_number)) %>%
  filter(!is.na(pen_number))


# Principal Component Analysis -------------------------------------------


cor(select(data, MPR, MS, AUC))

pca_res <- data %>%
  select(MPR, MS, AUC) %>%
  prcomp(scale. = TRUE)

explained_var_ratio <- summary(pca_res)[["importance"]]["Proportion of Variance",] * 100

pca_df <- as.data.frame(pca_res$x) %>%
  mutate(
    across(c(PC1), ~ -.),
    pen_number = data$pen_number,
    mpi = data$mpi
  )

pca_long <- pca_df %>%
  select(starts_with("PC"), pen_number, mpi) %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(starts_with("PC"), names_to = "PC", values_to = "value") 

pca_comp_df <- inner_join(
  pca_long %>% rename(x_pc = PC, x_val = value),
  pca_long %>% rename(y_pc = PC, y_val = value),
  relationship = "many-to-many"
) %>%
  filter(x_pc != y_pc)

# PCA Figure
pca_comp_df %>%
  ggplot(aes(x_val, y_val, color = pen_number)) +
  geom_point() +
  facet_grid(rows = vars(y_pc), cols = vars(x_pc), scales = "free") +
  labs(
    x = NULL, y = NULL,
    title = "Principal Component Pairs",
    color = "Pen"
  ) +
  main_theme
ggsave("pca.png", path="figures", width=12, height=12)

# PC1 with MPI and Pen #
scalar <- abs(min(-pca_long$value[which(pca_long$PC == "PC1")])) + 1

pca_long_sum <- pca_long %>%
  filter(PC == "PC1") %>%
  mutate(
    value = log(-value + scalar),
    inoculated = pen_number != "1"
  ) %>%
  summarize(
    mean_value = mean(value),
    stdev = sd(value),
    upper = mean_value + stdev,
    lower = mean_value - stdev,
    lower = ifelse(lower < 0, 0, lower),
    .by = c(mpi, inoculated)
  )

pca_long %>%
  filter(PC == "PC1") %>%
  mutate(
    value = log(-value + scalar),
    inoculated = pen_number != "1"
  ) %>%
  ggplot(aes(mpi, value, color = inoculated)) +
  geom_ribbon(aes(y=mean_value, ymin=lower, ymax=upper, fill = inoculated), data = pca_long_sum, alpha = 0.1, linetype=NA, show.legend = FALSE) +
  geom_line(aes(y=mean_value), data=pca_long_sum, alpha = 0.5, show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.25, position = position_jitter(0.2)) +
  scale_color_manual(values = c("navy", "red")) +
  scale_fill_manual(values = c("navy", "red")) +
  scale_x_continuous(breaks=seq(0, 60, 1)) +
  guides(color = guide_legend(override.aes = list(alpha = 1))) +
  # facet_grid(cols = vars(cutoff), rows = vars(Tissue)) +
  labs(
    x = "Months",
    y = sprintf("Log(pc1 + %s)", round(scalar, 2)),
    color = "Inoculated"
  ) +
  main_theme +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position = "top",
    legend.text=element_text(size=20),
    axis.text = element_text(size = 16)
  )
ggsave("pca_mpi.png", path="figures", width=12, height=9)



# Summary Figures ---------------------------------------------------------


data_sum <- data %>%
  mutate(PCA = pca_df$PC1) %>%
  summarize(
    across(c(MPR, MS, AUC, PCA), list(mean = mean, stdev = sd)), 
    count = n(),
    .by = c(pen_number, mpi, week)
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
  mutate(across(matches("mean_[2-5]"), ~ . - mean_1)) %>%
  pivot_longer(matches("^[a-z]+_\\d$"), names_to = c(".value", "pen_number"), names_sep = "_") %>%
  mutate(pen_number = factor(pen_number, c("1", "2", "3", "4"), c("Pen 1", "Pen 2", "Pen 3", "Pen 4")))

data_norm %>%
  filter(
    !is.na(pen_number) & pen_number != "Pen 1", 
    !is.na(mean),
    metric == "PCA"
  ) %>%
  arrange(week) %>%
  ggplot(aes(week, mean, group = pen_number)) +
  geom_hline(
    yintercept = 0, linetype="dashed", linewidth = 1.5, color="blue"
  ) +
  geom_point(size = 4, color="red") +
  geom_segment(aes(xend = week, y = 0, yend = mean), color = "red", linewidth=2) +
  # stat_smooth(
  #   method = "glm",
  #   linewidth=1.5, color="red") +
  # stat_regline_equation(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), size = 6, label.x=0, label.y=5) +
  stat_cor(aes(label = paste(..p.label.., ..rr.label.., sep = "~','~")), size = 9, label.x=0, label.y = 5) +
  scale_x_continuous(breaks = seq(0, max(data_norm$week), 4)) +
  scale_y_continuous(limits = c(-7.5, 7.5), breaks=seq(-8, 8, 2), expand=expansion()) +
  facet_grid(rows = vars(pen_number), scales="free_y") +
  labs(
    x = "Weeks Post-Inoculation",
    y = "Normalized Mean PC1",
  ) +
  main_theme

ggsave("figures/metrics_over_time.png", width = 16, height = 12, dpi = 300)

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
