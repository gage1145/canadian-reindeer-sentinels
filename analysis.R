library(tidyverse)
library(janitor)
library(quicR)
library(stringr)
library(openxlsx)
library(readxl)


# Curation ----------------------------------------------------------------


files <- list.files("raw", ".xlsx", full.names = TRUE)

metadata <- read_xlsx("data/metadata.xlsx") %>%
  clean_names() %>%
  mutate(
    across(c("sample_id", "random_id", "pen_number"), as.factor),
    across(c("week", "mpi"), as.integer)
  ) %>%
  select(sample_id, random_id, pen_number, week, mpi, cell)

curate <- function(file) {
  filename <- file %>%
    str_split_i("/", str_count(., "/") + 1) %>%
    str_remove(".xlsx")
  cli::cli_inform(filename)
  get_quic(file) %>%
    mutate(
      `Sample IDs` = trimws(as.character(`Sample IDs`)),
      `Rxn IDs` = filename
    ) %>%
    relocate(`Rxn IDs`, .after = `Sample IDs`) %>%
    suppressMessages()
}

df_ <- map_dfr(files, curate)

df_clean <- df_ %>%
  clean_names() %>%
  rename("sample_id" = 1, "rxn" = 2)

write.csv(df_clean, "data/data.csv", row.names = FALSE)


# Calculations ------------------------------------------------------------


groupings <- c("sample_id", "dilutions", "wells", "rxn")

calcs <- df_clean %>%
  calculate_metrics(
    groupings,
    threshold = 3,
    time_col = "time",
    ttt_values = "norm",
    auc_values = "norm",
    norm_col = "norm",
    deriv_col = "deriv"
  ) %>%
  left_join(metadata) %>%
  arrange(random_id)
write.csv(calcs, "data/calculations.csv", row.names = FALSE)


# Analysis ----------------------------------------------------------------


library(lme4)
library(emmeans)

df_sum <- calcs %>%
  filter(!is.na(pen_number))

auc_formula <- AUC ~ mpi * pen_number + (1 | random_id) + (1 | rxn)
mpr_formula <- MPR ~ mpi * pen_number + (1 | random_id) + (1 | rxn)
ms_formula  <- MS  ~ mpi * pen_number + (1 | random_id) + (1 | rxn)

auc_mod <- lmer(auc_formula, df_sum)
mpr_mod <- lmer(mpr_formula, df_sum)
ms_mod  <- lmer(ms_formula,  df_sum)

mpis <- sort(unique(df_sum$mpi))

auc_emm <- emmeans(auc_mod, ~ pen_number | mpi, at = list(mpi = mpis))
mpr_emm <- emmeans(mpr_mod, ~ pen_number | mpi, at = list(mpi = mpis))
ms_emm  <- emmeans(ms_mod,  ~ pen_number | mpi, at = list(mpi = mpis))

auc_con <- contrast(auc_emm, method = "trt.vs.ctrl", ref = 1)
mpr_con <- contrast(mpr_emm, method = "trt.vs.ctrl", ref = 1)
ms_con  <- contrast(ms_emm,  method = "trt.vs.ctrl", ref = 1)

auc_sum <- summary(auc_con, adjust = "none", side = ">") %>% mutate(metric = "AUC")
mpr_sum <- summary(mpr_con, adjust = "none", side = ">") %>% mutate(metric = "MPR")
ms_sum  <- summary(ms_con,  adjust = "none", side = ">") %>% mutate(metric = "MS")

results <- bind_rows(auc_sum, mpr_sum, ms_sum) %>%
  arrange(mpi, contrast, metric) %>%
  group_by(mpi) %>%
  mutate(p.adj = p.adjust(p.value, method = "holm")) %>%
  ungroup()

write.csv(results, "data/results.csv", row.names = FALSE)
