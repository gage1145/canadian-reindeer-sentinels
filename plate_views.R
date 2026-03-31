library(quicR)
library(tidyverse)


files <- list.files("raw", ".xlsx", full.names = TRUE, recursive = TRUE)

save_plates <- function(file) {
  filename <- file %>%
    str_remove(".xlsx") %>%
    str_split_i("/", str_count(., "/") + 1)

  full_name <- paste0("figures/plateviews/", filename, "_plateview.png")

  show_dilution <- str_detect(filename, "block")

  if (!file.exists(full_name)) {
    plate_view(get_quic(file), sep = "\n", plot_deriv = FALSE) +
      ggtitle(filename) +
      theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    ggsave(
      full_name,
      width = 12, height = 8
    )
  } else {
    cat(paste("\nPlate view for", file, "already exists. "))
    return()
  }
}

lapply(files, save_plates)