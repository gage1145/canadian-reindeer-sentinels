library(dplyr)
library(tidyr)
library(quicR)
library(stringr)
library(openxlsx)
library(readxl)



# Curation ----------------------------------------------------------------



files <- list.files("raw", ".xlsx", full.names = TRUE)

metadata <- read_xlsx("data/metadata.xlsx") %>%
  rename(`Sample IDs` = `Sample #`) %>%
  mutate_at("Sample IDs", as.character)

curate <- function(file) {
  filename <- file %>%
    str_split_i("/", 2) %>%
    str_split_i("\\.", 1)
  get_real(file, ordered = FALSE)[[1]] %>%
    transpose_real() %>%
    mutate(
      `Sample IDs` = convert_tables(organize_tables(file))$`Sample IDs`,
      `Rxn IDs` = filename
    ) %>%
    relocate(`Rxn IDs`, .after = `Sample IDs`)
}

data <- lapply(files, curate) %>%
  bind_rows()

write.csv(data, "data/data.csv", row.names = FALSE)



# Calculations ------------------------------------------------------------



norm <- normalize_RFU(data[-2], transposed = TRUE)
calcs <- data.frame(
  `Sample IDs` = trimws(norm$`Sample IDs`),
  `Rxn IDs` = data$`Rxn IDs`,
  check.names = FALSE
) %>%
  mutate(
    MPR = calculate_MPR(norm),
    MS  = calculate_MS(norm),
    TtT = calculate_TtT(norm, threshold = 2),
    RAF = 1/TtT,
    crossed = MPR >= 2
  ) %>%
  arrange(`Sample IDs`) %>%
  left_join(metadata, "Sample IDs")
write.csv(calcs, "data/calculations.csv", row.names = FALSE)



# Analysis ----------------------------------------------------------------



non_numerics <- calcs$`Sample IDs`[
  which(is.na(as.numeric(calcs$`Sample IDs`)))
] %>%
  unique()

numerics <- calcs$`Sample IDs` %>%
  as.numeric() %>%
  sort() %>%
  unique() %>%
  as.character()
  
samples <- c(non_numerics, numerics)

get_results <- function(samples, metric) {
  t_test <- function(sample) {
    data <- subset(calcs, `Sample IDs` %in% c(sample, "N"))
    t.test(
      data[[metric]], 
      data[which(data$`Sample IDs` == "N"),][[metric]],
      alternative = "greater"
    )
  }
  
  sapply(samples, t_test) %>%
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
  rename(MPR_pvalue = p.value)
results_MS  <- get_results(samples, "MS") %>%
  rename(MS_pvalue = p.value)



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
  summarise_at(vars(MPR:RAF), list(mean=mean, sd=sd)) %>%
  left_join(counts) %>%
  left_join(results_MPR) %>%
  left_join(results_MS) %>%
  mutate(
    `Sample IDs` = ordered(
      `Sample IDs`,
      levels = samples
    ),
    Result = ifelse(
      MPR_pvalue < 0.05 & MS_pvalue < 0.05 & thres_pos,
      "Detected",
      ifelse(
        thres_pos & (MPR_pvalue > 0.05 | MS_pvalue > 0.05),
        "Suspect",
        "Not Detected"
      )
    )
  ) %>%
  arrange(`Sample IDs`)

write.csv(results, "data/results.csv", row.names = FALSE)










