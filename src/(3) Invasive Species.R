# Connectivity - Nutrients Example

# Author(s): Reilly O'Connor
# Version: 2025-07-20

# Load Pkgs
library(tidyverse)
library(cowplot)

# load data





df_fert <- read.csv("../Connectivity/Data/Nutrients/Global_Fertilizer_Consumption.csv", header = T)
df_nitr <- read.csv("../Connectivity/Data/Nutrients/Nitrogen to Coastal Zone.csv", header = T)
df_wetland <- read.csv("../Connectivity/Data/Nutrients/Wetland_Loss.csv", header = T)









gg_eid <- ggplot(df_eid, aes(x = x, y = y)) +
  geom_smooth(se = F, color = "black") +
  labs(x = "Year",
       y = "Number of Local Disease Events",
       title = "Flux") +
  xlim(1950, 2010) +
  theme_bw(base_size = 14)

gg_eid

gg_airline <- ggplot(df_airplane_sum, aes(x = Year, y = sum_passengers/1e10)) +
  geom_smooth(se = F, color = "black") +
  labs(x = "Year",
       y = "Number of Airline Passengers (Billions)",
       title = "Connectivity") +
  xlim(1950, 2010) +
  theme_bw(base_size = 14)

gg_airline

gg_outbreaks <- ggplot(df_outbreaks, aes(x = x, y = y)) +
  geom_smooth(se = F, color = "black") +
  labs(x = "Year",
       y = "Number of Large Scale Disease Outbreaks",
       title = "Imbalance") +
  xlim(1950, 2010) +
  theme_bw(base_size = 14)

gg_outbreaks

gg_disease_grid <- plot_grid(gg_eid, gg_airline, gg_outbreaks, nrow = 1, align = "hv")

gg_disease_grid






