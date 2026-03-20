library(dplyr)
library(tidyr)
library(quicR)
library(stringr)
library(openxlsx)
library(readxl)



# Curation ----------------------------------------------------------------



files <- list.files("raw", ".xlsx", full.names = TRUE)

metadata <- read_xlsx("data/metadata.xlsx") %>%
  mutate_at("Sample IDs", as.character) %>%
  mutate_at("Random IDs", as.integer) %>%
  select(`Sample IDs`, `Random IDs`, `Pen #`, Week, MPI)

curate <- function(file) {
  filename <- file %>%
    str_split_i("/", 2) %>%
    str_split_i("\\.", 1)
  cli::cli_inform(filename)
  get_quic(file) %>%
    mutate(
      `Sample IDs` = as.character(`Sample IDs`),
      `Rxn IDs` = filename
    ) %>%
    relocate(`Rxn IDs`, .after = `Sample IDs`) %>%
    suppressMessages()
}

data <- lapply(files, curate) %>%
  bind_rows()

write.csv(data, "data/data.csv", row.names = FALSE)



# Calculations ------------------------------------------------------------



# norm <- normalize_RFU(data[-2], transposed = TRUE)
calcs <- data %>%
  calculate_metrics(
    "Sample IDs", "Dilutions", "Wells", "Rxn IDs", threshold = 3
  ) %>%
  left_join(metadata) %>%
  # filter(`Sample IDs` %in% c("N", "P") | `Random IDs` %in% 117:228) %>%
  arrange(`Random IDs`)
write.csv(calcs, "data/calculations.csv", row.names = FALSE)



# Analysis ----------------------------------------------------------------



# non_numerics <- calcs$`Sample IDs`[
#   which(is.na(as.numeric(calcs$`Sample IDs`)))
# ] %>%
#   unique()
# 
# numerics <- calcs$`Sample IDs` %>%
#   as.numeric() %>%
#   sort() %>%
#   unique() %>%
#   as.character()
  
samples <- unique(calcs$`Sample IDs`)

get_results <- function(samples, metric) {
  t_test <- function(sample) {
    data <- filter(calcs, `Sample IDs` %in% c(sample, "N"))
    t.test(
      data[[metric]], 
      data[which(data$`Sample IDs` == "N"),][[metric]],
      alternative = "greater"
    )
  }
  
  sapply(samples, t_test) %>%
    # bind_rows() %>%
    t() %>%
    as.data.frame() %>%
    mutate(
      `Sample IDs` = ordered(
        rownames(.),
        levels = samples
      ),
      p.value = signif(as.numeric(p.value), 3)
    ) %>%
    relocate(`Sample IDs`, 1) %>%
    select(`Sample IDs`, p.value) %>%
    arrange(`Sample IDs`)
}

results_MPR <- get_results(samples, "MPR") %>%
  rename(MPR_pvalue = p.value) %>%
  left_join(metadata)
results_MS  <- get_results(samples, "MS") %>%
  rename(MS_pvalue = p.value) %>%
  left_join(metadata)
results_AUC  <- get_results(samples, "AUC") %>%
  rename(AUC_pvalue = p.value) %>%
  left_join(metadata)



# Summary -----------------------------------------------------------------


counts <- calcs %>%
  group_by(`Sample IDs`) %>%
  summarize(
    counts = n(),
    crossed = sum(crossed) / counts,
    thres_pos = crossed >= 0.5
  )

results <- calcs %>%
  group_by(`Sample IDs`, `Pen #`, Week) %>%
  summarise_at(c("MPR", "MS", "AUC", "TtT", "RAF"), list(mean=mean, sd=sd)) %>%
  left_join(counts) %>%
  left_join(results_MPR) %>%
  left_join(results_MS) %>%
  left_join(results_AUC) %>%
  mutate(
    # `Sample IDs` = ordered(
    #   `Sample IDs`,
    #   levels = samples
    # ),
    Result = ifelse(
      MPR_pvalue < 0.05 & MS_pvalue < 0.05 & AUC_pvalue < 0.05 & thres_pos,
      "Detected",
      ifelse(
        thres_pos & (MPR_pvalue > 0.05 | MS_pvalue > 0.05 | AUC_pvalue > 0.05),
        "Suspect",
        "Not Detected"
      )
    )
  ) %>%
  arrange(`Sample IDs`)

write.csv(results, "data/results.csv", row.names = FALSE)










