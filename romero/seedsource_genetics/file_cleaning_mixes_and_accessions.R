library (tidyverse)
library (readxl)
library (dplyr)
library (stringr)
library (lubridate)

setwd("C:/Users/Seth/Documents/seed_mix_genetics")

culumber_LECI_accessions <- read_xlsx("csc2cropsci2012060396-sup-0002.xlsx", skip = 2) %>%
  rename(accession_id = "Sample ID",
         bayesian_group = "Bayesian K=4 group",
         ploidy_num = "Ploidy (2n)",
         LAT = Latitude,
         LONG = Longitude) %>%
  filter(str_starts(.$accession_id, "Lcin"),
         LAT != "-",
         LONG != "-",
         ploidy_num != "-",
         bayesian_group != "-") %>%
  mutate(race = ifelse(bayesian_group == "1", "Rocky Mtn.",
                       ifelse(bayesian_group == "2", "Columbia", "Great Basin")),
         ploidy_chr = ifelse(ploidy_num == "28", "Tetraploid", "Octoploid"),
         LAT = as.numeric(LAT),
         LONG = as.numeric(LONG))

gb_psz <- data.frame(PROV_SEED_ZONE = c("10 - 15 Deg. F. / 6 - 12",
                                        "15 - 20 Deg. F. / 3 - 6",
                                        "15 - 20 Deg. F. / 6 - 12",
                                        "20 - 25 Deg. F. / 6 - 12"),
                     gb_psz = c("Zone 8", "Zone 11", "Zone 12", "Zone 16"))

# NOTE: two different collection dates for Elko, NV pop. (2001 & 2004) - only keeping 2001 record for simplicity
pops_LECI_seed_mix <- read_xlsx("Benson LECI4 Location and Collection records.xlsx", sheet = 2) %>%
  filter(str_starts(.$Lot_Number, "LECI")) %>%
  mutate(COLLECT_DATE = as_date(as.numeric(COLLECT_DATE)) - years(70)) %>%
  distinct(Location_Name, .keep_all = TRUE)

ploidy_LECI <- read_xlsx("Basin Wildrye Ploidy.xlsx", skip = 2) %>%
  rename(Location_Name = "Location Name",
         notes = ...8) %>%
  filter(is.na(Location_Name) == FALSE)

combined_LECI_seed_mix <- pops_LECI_seed_mix %>%
  full_join(ploidy_LECI, by = "Location_Name") %>%
  select(-'Provisional Seed Zone', -'Collection Lot Number') %>%
  rename(tz_2014_pct = '2014 TZ (%)',
         seed_weight = 'Weight (.lbs)',
         ploidy_chr = Ploidy) %>%
  full_join(gb_psz, by = "PROV_SEED_ZONE") %>%
  arrange(as.numeric(str_sub(gb_psz, 6))) %>%
  group_by(gb_psz) %>%
  mutate(prop_weight = seed_weight/sum(seed_weight))

mapping_data <- bind_rows(combined_LECI_seed_mix, culumber_LECI_accessions) %>%
  mutate(color_group = ifelse(is.na(gb_psz) == FALSE, gb_psz, race),
         LONG = ifelse(LONG > 0, LONG*-1, LONG)) %>%
  rename(X = LONG,
         Y = LAT) %>%
  arrange(X)

write.csv(mapping_data, "LECI_mapping_data.csv")
