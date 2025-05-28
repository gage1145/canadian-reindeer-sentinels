library(dplyr)
library(tidyr)
library(quicR)
library(stringr)
library(openxlsx)
library(readxl)
library(ggplot2)


# Curation ----------------------------------------------------------------



files <- list.files("raw/salt_block", ".xlsx", full.names = TRUE)

curate <- function(x) {
  meta <- convert_tables(organize_tables(x))
  get_real(x, ordered = FALSE)[[1]] %>%
    transpose_real() %>%
    mutate(
      `Sample IDs` = meta$`Sample IDs`,
      Dilutions = meta$Dilutions
    )
}

data <- lapply(files, curate) %>%
  bind_rows()
write.csv(data, "data/salt_data.csv", row.names = FALSE)



# Calculations ------------------------------------------------------------


samples <- data$`Sample IDs` %>%
  factor(levels = c("P", "BP", "SP", "N", "B", "S"))
norm <- normalize_RFU(select(data, -Dilutions), transposed = TRUE)
calcs <- data.frame(
  `Sample IDs` = samples,
  Dilutions = factor(
    -log10(as.integer(data$Dilutions)),
    levels = seq(-1, -8)
  ),
  check.names = FALSE
) %>%
  mutate(
    MPR = calculate_MPR(norm),
    MS  = calculate_MS(norm),
    TtT = calculate_TtT(norm, threshold = 2),
    RAF = 1/TtT,
    crossed = MPR >= 2,
    pos_dilution = ifelse(
      `Sample IDs` == "P", 
      as.character(Dilutions),
      ifelse(str_detect(`Sample IDs`, "P"), -3, NA)
    ),
    block_dilution = ifelse(
      str_detect(`Sample IDs`, "B"),
      as.character(Dilutions),
      NA
    )
  ) 
  # arrange(`Sample IDs`)
write.csv(calcs, "data/salt_calculations.csv", row.names = FALSE)



# Summarize ---------------------------------------------------------------



df_sum <- calcs %>%
  group_by(`Sample IDs`, Dilutions, pos_dilution, block_dilution) %>%
  summarize(
    median_TtT = median(TtT)
  )

pos_ttt <- df_sum[[which(df_sum$pos_dilution == -3 & df_sum$`Sample IDs` == "P"), "median_TtT"]]
bp2_ttt <- df_sum[[which(df_sum$block_dilution == -2 & df_sum$`Sample IDs` == "BP"), "median_TtT"]]
bp3_ttt <- df_sum[[which(df_sum$block_dilution == -3 & df_sum$`Sample IDs` == "BP"), "median_TtT"]]
bp4_ttt <- df_sum[[which(df_sum$block_dilution == -4 & df_sum$`Sample IDs` == "BP"), "median_TtT"]]
sp_ttt  <- df_sum[[which(df_sum$`Sample IDs` == "SP"), "median_TtT"]]
s_ttt   <- df_sum[[which(df_sum$`Sample IDs` == "S"), "median_TtT"]]



# Plot Real-time ----------------------------------------------------------



# Positive curves
samp_order <- c(
  "P", "BP_-2", "BP_-3", "BP_-4", "SP", 
  "N", "B_-2", "B_-3", "B_-4", "S"
)

norm %>%
  # arrange(`Sample IDs`) %>%
  mutate(
    # Dilutions = calcs$Dilutions,
    pos_dilution = calcs$pos_dilution,
    block_dilution = calcs$block_dilution
  ) %>%
  unite("Sample IDs", `Sample IDs`, block_dilution, na.rm = TRUE) %>%
  group_by(`Sample IDs`, pos_dilution) %>%
  summarize_all(median) %>% 
  filter(pos_dilution == -3 | is.na(pos_dilution)) %>% 
  gather("Time (h)", "Normalized RFU", -c("Sample IDs", "pos_dilution")) %>%
  mutate_at(c("Time (h)", "Normalized RFU"), as.numeric) %>%
  mutate(
    `Sample IDs` = factor(`Sample IDs`, levels = samp_order),
    pos_dilution = ifelse(is.na(pos_dilution), "Not Spiked", "Spiked")
  ) %>% 
  arrange(`Sample IDs`) %>%
  ggplot(aes(`Time (h)`, `Normalized RFU`, color = `Sample IDs`)) +
  
  geom_vline(xintercept = pos_ttt, color = "#5b0085", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = bp2_ttt, color = "#c2003a", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = bp3_ttt, color = "#e13000", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = bp4_ttt, color = "#ffb800", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = sp_ttt,  color = "#a40084", linewidth = 1, linetype = "dashed") +
  # geom_vline(xintercept = s_ttt,   color = "#00c3ff", linewidth = 1, linetype = "dashed") +
  
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 2, linetype = "dashed", linewidth = 1.2) +
  facet_grid(rows = vars(pos_dilution)) +
  scale_color_manual(
    values = c(
      "#5b0085", "#c2003a", "#e13000", "#ffb800", "#a40084", 
      "#4f9e00", "#00af5d", "#00b6a0", "#00bed2", "#00c3ff"
    )
  ) +
  scale_x_log10(
    limits = c(5, 72), 
    breaks = seq(0, 75, 10),
    minor_breaks = seq(0, 75, 2),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, 14), 
    breaks = seq(0, 20, 2)
  ) +
  ggtitle("Real-time Positive Comparisons") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_blank(),
    legend.direction = "horizontal",
    legend.position = c(.8, .95),
    strip.text = element_text(size = 14)
  ) +
  guides(color=guide_legend(nrow=2, byrow = TRUE))
ggsave("figures/salt_real-time.png", width = 12, height = 10)



# # Plot Metrics ------------------------------------------------------------
# 
# 
# 
# calcs %>%
#   ggplot(aes(`Sample IDs`, RAF, fill = pos_dilution)) +
#   geom_boxplot(varwidth = TRUE)
# 
# 
# 
# # Titration Curve ---------------------------------------------------------
# 
# 
# 
# df_sum <- calcs %>%
#   # filter(Dilutions != -3) %>%
#   group_by(`Sample IDs`, pos_dilution) %>%
#   summarize(
#     mean_RAF = mean(RAF),
#     median_RAF = median(RAF),
#     sd_RAF = sd(RAF)
#   )
# 
# model <- glm(
#   RAF ~ pos_dilution,
#   data = calcs
#   # data = filter(calcs, `Sample IDs` == "P")
# ) %>%
#   summary()
# 
# r_square <- with(model, 1 - deviance / null.deviance)
# 
# # Plot the titration curve
# sample_key <- list(
#   "Positive Control" = "P",
#   "Block + Positive" = "BP",
#   "Swab + Positive"  = "SP",
#   "Negative Control" = "N",
#   "Block Only"       = "B",
#   "Swab Only"        = "S"
# ) %>% unlist()
# 
# calcs %>%
#   filter(`Sample IDs` %in% c("P", "BP", "SP")) %>%
#   mutate (
#     `Sample IDs` = names(sample_key)[match(`Sample IDs`, sample_key)] %>%
#       factor(levels = names(sample_key))
#   ) %>%
#   # na.omit() %>%
#   # mutate(pos_dilution = as.integer(pos_dilution)) %>%
#   # filter(`Sample IDs` == "P") %>%
#   ggplot(aes(
#     pos_dilution, RAF, fill = `Sample IDs`, shape = `Sample IDs`
#     # ymin = median_RAF - sd_RAF, ymax = median_RAF + sd_RAF
#   )) +
#   geom_boxplot(width = 0.5, outlier.shape = NA) +
#   geom_point(aes(fill = `Sample IDs`), size = 2, stroke = 0.5) +
#   # scale_fill_manual(values = c("#d10000", "#e60079", "#ba57d3", "#578bf8", "#00a8e0", "#00b8b8")) +
#   facet_grid(~ block_dilution, scales = "free")
#   # geom_text(aes(label = `Sample IDs`)) +
#   # geom_point(size = 3, position = "jitter") +
#   # geom_errorbar(width = 0.05, position = "jitter") +
#   # geom_smooth(
#   #   method = glm,
#   #   se = FALSE,
#   #   linetype = "dashed",
#   #   colour = "black"
#   # ) 
#   # geom_text(label = expression(R^2 == 0.616), x = -6.5, y = 3e-5) +
#   # ggtitle("Positive Titration Curve") +
#   # xlab("Log Dilution Factor") +
#   # ylab("Rate of Amyloid Formation (1/s)") 
#   # scale_y_continuous(breaks = seq(0, 3e-5, 5e-6)) +
#   # coord_cartesian(xlim = c(-7.2, -2.8), ylim = c(0, 3.3e-5), expand = FALSE) +
#   # theme_classic() +
#   # theme(
#   #   legend.position = "none",
#   #   axis.line = element_line(colour = "black"),
#   #   axis.text = element_text(size = 12),
#   #   axis.title = element_text(size = 14),
#   #   plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#   #   panel.background = element_blank(),
#   #   panel.grid = element_line(colour = "lightgrey"),
#   #   panel.border = element_rect(colour = "black", fill = NA, size = 1),
#   #   strip.background = element_blank(),
#   #   strip.placement = "outside"
#   # )






