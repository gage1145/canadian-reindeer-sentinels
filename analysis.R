library(dplyr)
library(tidyr)
library(quicR)
library(stringr)
library(openxlsx)
library(readxl)
library(lme4)
library(emmeans)



# Curation ----------------------------------------------------------------



files <- list.files("raw", ".xlsx", full.names = TRUE)

curate <- function(x) {
  get_real(x, ordered = FALSE)[[1]] %>%
    transpose_real() %>%
    mutate(`Sample IDs` = convert_tables(organize_tables(x))$`Sample IDs`)
}

data <- lapply(files, curate) %>%
  bind_rows()
write.csv(data, "data/data.csv", row.names = FALSE)



# Calculations ------------------------------------------------------------



norm <- normalize_RFU(data, transposed = TRUE)
calcs <- data.frame(
  `Sample IDs` = trimws(norm$`Sample IDs`),
  check.names = FALSE
) %>%
  mutate(
    MPR = calculate_MPR(norm),
    MS  = calculate_MS(norm),
    TtT = calculate_TtT(norm, threshold = 2),
    RAF = 1/TtT,
    crossed = MPR >= 2
  ) %>%
  arrange(`Sample IDs`)
write.csv(calcs, "data/calculations.csv", row.names = FALSE)



# Analysis ----------------------------------------------------------------



samples <- unique(calcs$`Sample IDs`)



get_results <- function(samples, metric) {
  
  t_test <- function(name, metric) {
    data <- subset(calcs, `Sample IDs` %in% c(name, "N"))
    t.test(
      data[[metric]], 
      data[which(data$`Sample IDs` == "N"),][[metric]],
      alternative = "greater"
    )
  }
  
  sapply(samples, t_test, metric = metric) %>%
    t() %>%
    as.data.frame() %>%
    mutate(
      `Sample IDs` = ordered(
        rownames(.),
        levels = c("P", "N", as.character(seq(1, 88)))
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
    thres_pos = sum(crossed) / counts > 0.5
  )

results <- calcs %>%
  group_by(`Sample IDs`) %>%
  summarise_at(vars(MPR:RAF), list(mean=mean, sd=sd, min=min, max=max)) %>%
  left_join(counts) %>%
  left_join(results_MPR) %>%
  left_join(results_MS) %>%
  mutate(
    `Sample IDs` = ordered(
      `Sample IDs`,
      levels = c("P", "N", as.character(seq(1, 88)))
    ),
    Result = ifelse(
      MPR_pvalue < 0.05 & MS_pvalue < 0.05 & thres_pos,
      "Detected",
      ifelse(
        (MPR_pvalue < 0.05) + (MS_pvalue < 0.05) + thres_pos == 2,
        "Strong Suspect",
        ifelse(
          (MPR_pvalue < 0.05) + (MS_pvalue < 0.05) + thres_pos == 1,
          "Weak Suspect",     
          "Not Detected"
        )
      )
    )
  ) %>%
  arrange(`Sample IDs`)

write.csv(results, "data/results.csv", row.names = FALSE)










