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
df_dz <- read.csv("../Connectivity/Data/Nutrients/Num_Deadzones.csv", header = T)


gg_fert <- ggplot(df_fert, aes(x = Year, y = World_Historic)) +
  geom_smooth(se = F, color = "black") +
  labs(x = "Year",
       y = "Global Fertilizer Consumption (Tonnes/Yr)",
       title = "Flux") +
  xlim(1850, 2010) +
  theme_bw(base_size = 14)

gg_fert


gg_nitr <- ggplot(df_nitr, aes(x = Year, y = Nitrogen.flux.Mtons.yr.1)) +
  geom_smooth(se = F, color = "black") +
  labs(x = "Year",
       y = "Nitrogen Flux (Mtons per year)",
       title = "Flux") +
  xlim(1850, 2010) +
  theme_bw(base_size = 14)

gg_nitr

gg_wetland <- ggplot(df_wetland, aes(x = x, y = y)) +
  geom_smooth(se = F, color = "black") +
  labs(x = "Year",
       y = "Global Percent Wetland Loss",
       title = "Connectivity") +
  xlim(1850, 2010) +
  theme_bw(base_size = 14)

gg_wetland

gg_dz <- ggplot(df_dz, aes(x = x, y = y)) +
  geom_smooth(se = F, color = "black") +
  labs(x = "Year",
       y = "Reported # of Deadzones",
       title = "Imbalance") +
  xlim(1850, 2010) +
  theme_bw(base_size = 14)

gg_dz

gg_nutrient_grid <- plot_grid(gg_nitr, gg_wetland, gg_dz, nrow = 1, align = "hv", labels = c("a", "b", "c"))

gg_nutrient_grid






