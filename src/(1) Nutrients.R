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
  # geom_point(size = 2, color = "grey20", color = "black", stroke = 0.5, alpha = 0.2) +
  geom_smooth(se = F, color = "#F79D1E", linewidth = 3) +
  labs(x = "Year") +
  xlim(1850, 2010) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())

gg_fert


gg_nitr <- ggplot(df_nitr, aes(x = Year, y = Nitrogen.flux.Mtons.yr.1)) +
  # geom_point(size = 2, color = "grey20", color = "black", stroke = 0.5, alpha = 0.2) +
  geom_smooth(se = F, color = "#F79D1E", linewidth = 3) +
  labs(x = "Year") +
  xlim(1850, 2010) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())

gg_nitr

gg_wetland <- ggplot(df_wetland, aes(x = x, y = y)) +
  # geom_point(size = 2, color = "grey20", color = "black", stroke = 0.5, alpha = 0.2) +
  geom_smooth(se = F, color = "#075149", linewidth = 3) +
  labs(x = "Year") +
  xlim(1850, 2010) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())

gg_wetland

gg_dz <- ggplot(df_dz, aes(x = x, y = y)) +
  # geom_point(size = 2, color = "grey20", color = "black", stroke = 0.5, alpha = 0.2) +
  geom_smooth(se = F, color = "#800000", linewidth = 3) +
  labs(x = "Year") +
  xlim(1850, 2010) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())

gg_dz

gg_nutrient_grid <- plot_grid(gg_fert, gg_wetland, gg_dz, nrow = 1, align = "hv", labels = c("a", "b", "c"))

gg_nutrient_grid

ggsave("../Connectivity/Figures/Figure 3a.jpeg", plot = gg_fert, dpi = 300, width = 3, height = 3)
ggsave("../Connectivity/Figures/Figure 3b.jpeg", plot = gg_wetland, dpi = 300, width = 3, height = 3)
ggsave("../Connectivity/Figures/Figure 3c.jpeg", plot = gg_dz, dpi = 300, width = 3, height = 3)





