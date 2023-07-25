# Plotting steelhead returns at Bonneville

library(tidyverse)
library(here)

BON_passage_counts <- read.csv(here::here("figures", "exploratory", "adultannual_1689958375_841.csv"))

# replace 0s with NAs so that we start the time series at the right point
# Drop 2023 because the run year is not yet complete
BON_passage_counts %>% 
  dplyr::select(Year, Steelhead, Wild.Steelhead) %>% 
  pivot_longer(cols = c("Steelhead", "Wild.Steelhead")) %>% 
  mutate(value = ifelse(value == 0, NA, value)) %>% 
  filter(Year <= 2022) -> BON_steelhead_counts


ggplot(BON_steelhead_counts, aes(x = Year, y = value, color = name)) +
  geom_line()


# Conclusions: 
# Steelhead are not doing well. Right around the full time series low


# How does that compare to other salmonids?
BON_passage_counts %>% 
  dplyr::select(Year, Chinook, Steelhead, Sockeye, Coho) %>% 
  pivot_longer(cols = c("Chinook", "Steelhead", "Sockeye", "Coho")) %>% 
  mutate(value = ifelse(value == 0, NA, value)) %>% 
  filter(Year <= 2022) -> BON_salmonid_counts

ggplot(BON_salmonid_counts, aes(x = Year, y = value, color = name)) +
  geom_line()

# Sockeye are doing really well - 2022 was the all time high.
# Chinook and Coho are also doing well - above the long term average.
# Steelhead? Worst of the bunch. There was a time when Steelhead and Cnihook were 
# jointly the most abundant salmonids, but now Steelhead are the least 
# (and Chinook are the still the most)


# So why are Steelhead doing so badly? What makes them different from Chinook and Coho and Sockeye?
# Well... one part of it might be their adult life history.