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
    str_split_i("/", 2) %>%
    str_split_i("\\.", 1)
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
    "sample_id", "dilutions", "wells", "rxn", 
    threshold = 3, time_col = "time", ttt_values = "norm", 
    auc_values = "norm", norm_col = "norm", deriv_col = "deriv"
  ) %>%
  left_join(metadata) %>%
  arrange(random_id)
write.csv(calcs, "data/calculations.csv", row.names = FALSE)



# Analysis ----------------------------------------------------------------


  
# samples <- unique(calcs$sample_id)

# get_results <- function(samples, metric) {
#   t_test <- function(sample) {
#     data <- filter(calcs, `Sample IDs` %in% c(sample, "N"))
#     t.test(
#       data[[metric]], 
#       data[which(data$`Sample IDs` == "N"),][[metric]],
#       alternative = "greater"
#     )
#   }
  
#   sapply(samples, t_test) %>%
#     # bind_rows() %>%
#     t() %>%
#     as.data.frame() %>%
#     mutate(
#       `Sample IDs` = ordered(
#         rownames(.),
#         levels = samples
#       ),
#       p.value = signif(as.numeric(p.value), 3)
#     ) %>%
#     relocate(`Sample IDs`, 1) %>%
#     select(`Sample IDs`, p.value) %>%
#     arrange(`Sample IDs`)
# }

# results_MPR <- get_results(samples, "MPR") %>%
#   rename(MPR_pvalue = p.value) %>%
#   left_join(metadata)
# results_MS  <- get_results(samples, "MS") %>%
#   rename(MS_pvalue = p.value) %>%
#   left_join(metadata)
# results_AUC  <- get_results(samples, "AUC") %>%
#   rename(AUC_pvalue = p.value) %>%
#   left_join(metadata)



# Summary -----------------------------------------------------------------


# counts <- calcs %>%
#   group_by(`Sample IDs`) %>%
#   summarize(
#     counts = n(),
#     crossed = sum(crossed) / counts,
#     thres_pos = crossed >= 0.5
#   )

# results <- calcs %>%
#   group_by(`Sample IDs`, `Pen #`, Week) %>%
#   summarise_at(c("MPR", "MS", "AUC", "TtT", "RAF"), list(mean=mean, sd=sd)) %>%
#   left_join(counts) %>%
#   left_join(results_MPR) %>%
#   left_join(results_MS) %>%
#   left_join(results_AUC) %>%
#   mutate(
#     # `Sample IDs` = ordered(
#     #   `Sample IDs`,
#     #   levels = samples
#     # ),
#     Result = ifelse(
#       MPR_pvalue < 0.05 & MS_pvalue < 0.05 & AUC_pvalue < 0.05 & thres_pos,
#       "Detected",
#       ifelse(
#         thres_pos & (MPR_pvalue > 0.05 | MS_pvalue > 0.05 | AUC_pvalue > 0.05),
#         "Suspect",
#         "Not Detected"
#       )
#     )
#   ) %>%
#   arrange(`Sample IDs`)

# write.csv(results, "data/results.csv", row.names = FALSE)










