library(dplyr)
library(tidyr)
library(stringr)

setwd("~/Documents/GitHub/lab/romero/KRLA/samples/")

krla <- read.csv("KRLA_SW_202210_samples.csv") %>%
  filter(sample_ID != "EMPTY") %>%
  mutate(nums = str_sub(cell_location, start = 2),
         lets = str_sub(cell_location, end = 1),
         plt = str_sub(well_plate, start = 6),
         rep = seq(1, 497, by = 1),
         trueID = paste("KL", str_sub(sample_ID, end = 2), as.character(rep), sep = "_")) %>%
  pivot_wider(id_cols = c("plt", "lets"), names_from = nums, values_from = trueID)

write.csv(krla, "KRLA_extraction_plates.csv", row.names = FALSE)
