library(quicR)
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
library(ggh4x)





# Plate Views -------------------------------------------------------------



files <- list.files("raw", ".xlsx", full.names = TRUE, recursive = TRUE)

save_plates <- function(file) {
  
  filename <- file %>%
    str_remove(".xlsx") %>%
    str_split_i("/", str_count(., "/") + 1)
  
  full_name <- paste0("figures/plateviews/", filename, "_plateview.png")
  
  show_dilution <- str_detect(filename, "block")
  
  if (!file.exists(full_name)) {
    plate_view(
      get_real(file), 
      get_sample_locations(
        file, 
        dilution_bool = show_dilution, 
        dilution_fun = function(x) -log10(x)
      )
    ) +
      ggtitle(filename) +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    ggsave(
      full_name,
      width = 12, height = 8
    )
  } else{
    cat(paste("\nPlate view for", file, "already exists. "))
    return()
  }
}

lapply(files, save_plates)



# Summary Figures ---------------------------------------------------------


data <- read.csv("data/calculations.csv", check.names = FALSE)
results <- read.csv("data/results.csv", check.names = FALSE) %>%
  select(`Sample IDs`, Result)

# For factor levels
non_numerics <- data$`Sample IDs`[
  which(is.na(as.numeric(data$`Sample IDs`)))
] %>%
  unique()
numerics <- data$`Sample IDs` %>%
  as.numeric() %>%
  sort() %>%
  unique() %>%
  as.character()
samples <- c(non_numerics, numerics)

data <- left_join(data, results) %>%
  mutate(
    `Sample IDs` = ordered(
      `Sample IDs`, levels = samples
    )
  ) %>%
  arrange(`Sample IDs`)

y_scales = list(
  scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, 2)),
  scale_y_continuous(limits = c(0, 4)),
  scale_y_continuous(limits = c(0, 0.11), breaks = seq(0, 0.11, 0.02)),
  scale_y_continuous(limits = c(0, 72), breaks = seq(0, 72, 12))
)

data %>%
  select(-crossed) %>%
  mutate(
    Result = factor(
      Result, 
      levels = c("Not Detected", "Suspect", "Detected")
    ),
    Week = as.factor(Week)
  ) %>%
  filter(!is.na(`Pen #`)) %>%
  # gather("key", "value", -c("Sample IDs", "Result", "Rxn IDs")) %>%
  ggplot(aes(
    # fct_inorder(`Sample IDs`), 
    Week,
    RAF, 
    fill = Result
  )) +
  geom_boxplot() +
  scale_fill_manual(values = c("#08D2FF", "#5F5589", "#99023B")) +
  facet_grid(
    rows=vars(`Pen #`), 
    # cols=vars(Week), 
    scales = "free_x"
  ) +
  # facetted_pos_scales(y = y_scales) +
  labs(
    # x = "Sample IDs",
    x = "Week",
    y = "Rate of Amyloid Formation (1/s)",
    title = "PrioScout Swab Results"
  ) +
  theme(
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 9, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# ggsave("summary.png", path = "figures", width = 10.43, height = 15)
ggsave("summary.png", path = "figures", width = 12, height = 8)



# No pen IDs --------------------------------------------------------------


  
data %>%
  select(-crossed) %>%
  mutate(
    Result = factor(
      Result, 
      levels = c("Not Detected", "Suspect", "Detected")
    ),
    Week = as.factor(Week)
  ) %>%
  filter(is.na(`Pen #`)) %>%
  # gather("key", "value", -c("Sample IDs", "Result", "Rxn IDs")) %>%
  ggplot(aes(
    fct_inorder(`Sample IDs`),
    # Week,
    RAF, 
    fill = Result
  )) +
  geom_boxplot() +
  scale_fill_manual(values = c("#08D2FF", "#5F5589", "#99023B")) +
  labs(
    x = "Sample IDs",
    # x = "Week",
    y = "Rate of Amyloid Formation (1/s)",
    title = "PrioScout Swab Results"
  ) +
  theme(
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 9, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
# ggsave("summary.png", path = "figures", width = 10.43, height = 15)
ggsave("summary_no_pens.png", path = "figures", width = 12, height = 8)  

  
  
  
  
  
  