# Connectivity - Nutrients Example

# Author(s): Reilly O'Connor
# Version: 2025-07-20

# Load Pkgs
library(tidyverse)
library(cowplot)

# load data
df_alien <- read.csv("../Connectivity/Data/Invasives/Figure 2.04_Data_TimeSeries_Aliens.csv", header = TRUE)

df_ships <- read.table("../Connectivity/Data/Invasives/Figure 2.01_Data_DriverTimeSeries.csv", 
                   header = TRUE, sep = "", stringsAsFactors = FALSE)

df_invasive <- read.csv("../Connectivity/Data/Invasives/Negative_Repoted Impacts.csv", header = TRUE)

df_alien_count <- df_alien %>%
  group_by(FirstRecord) %>%
  summarise(number_aliens = n()) %>%
  rename(Year = FirstRecord) %>%
  filter(Year > 1499) %>%
  arrange(Year)

gg_alien <- ggplot(df_alien_count, aes(x = Year, y = number_aliens)) +
  geom_smooth(se = F, color = "#FFC107", linewidth = 3) +
  labs(x = "Year") +
  xlim(1700, 2010) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())

gg_alien

gg_ships <- ggplot(df_ships, aes(x = year, y = shipcalls)) +
  geom_smooth(se = F, color = "#1E88E5", linewidth = 3) +
  labs(x = "Year") +
  xlim(1700, 2010) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())

gg_ships

gg_invasive <- ggplot(df_invasive, aes(x = x, y = y)) +
  geom_smooth(se = F, color = "#D81B60", linewidth = 3) +
  labs(x = "Year") +
  xlim(1700, 2010) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())

gg_invasive

gg_invasive_grid <- plot_grid(gg_alien, gg_ships, gg_invasive, nrow = 1, align = "hv", labels = c("g", "h", "i"))

gg_invasive_grid

ggsave("../Connectivity/Figures/Figure 3g.jpeg", plot = gg_alien, dpi = 300, width = 3, height = 3)
ggsave("../Connectivity/Figures/Figure 3h.jpeg", plot = gg_ships, dpi = 300, width = 3, height = 3)
ggsave("../Connectivity/Figures/Figure 3i.jpeg", plot = gg_invasive, dpi = 300, width = 3, height = 3)








