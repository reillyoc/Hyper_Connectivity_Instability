# Connectivity - Nutrients Example

# Author(s): Reilly O'Connor
# Version: 2025-07-20

# Load Pkgs
library(tidyverse)
library(cowplot)
library(reshape2)

# load data
df_eid <- read.csv("../Connectivity/Data/Disease/Local_EID_Outbreaks.csv", header = T)
df_airplane <- read.csv("../Connectivity/Data/Disease/global_airline_passengers.csv", header = T)
df_outbreaks <- read.csv("../Connectivity/Data/Disease/Global_Outbreaks.csv", header = T)

#Transpose Airline Data
df_airplane_melt <- melt(df_airplane %>% select(- Country.Code, - Indicator.Name, - Indicator.Code),
                         id.vars = c("Country.Name"),
                         variable.name = "Year",
                         value.name = "airline_passengers")

df_airplane_melt$Year <- sub("^X", "", df_airplane_melt$Year)

df_airplane_sum <- df_airplane_melt %>%
  group_by(Year) %>%
  summarize(sum_passengers = sum(airline_passengers, na.rm = T)) %>%
  mutate(Year = as.numeric(Year))

gg_eid <- ggplot(df_eid, aes(x = x, y = y)) +
  geom_smooth(se = F, color = "#FFC107", linewidth = 3) +
  labs(x = "Year") +
  xlim(1950, 2010) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())

gg_eid

gg_airline <- ggplot(df_airplane_sum, aes(x = Year, y = sum_passengers/1e10)) +
  geom_smooth(se = F, color = "#1E88E5", linewidth = 3) +
  labs(x = "Year") +
  xlim(1950, 2010) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())

gg_airline

gg_outbreaks <- ggplot(df_outbreaks, aes(x = x, y = y)) +
  geom_smooth(se = F, color = "#D81B60", linewidth = 3) +
  labs(x = "Year") +
  xlim(1950, 2010) +
  theme_bw(base_size = 14) +
  theme(axis.title.y = element_blank())

gg_outbreaks

gg_disease_grid <- plot_grid(gg_eid, gg_airline, gg_outbreaks, nrow = 1, align = "hv", labels = c("d", "e", "f"))

gg_disease_grid

ggsave("../Connectivity/Figures/Figure 3d.jpeg", plot = gg_eid, dpi = 300, width = 3, height = 3)
ggsave("../Connectivity/Figures/Figure 3e.jpeg", plot = gg_airline, dpi = 300, width = 3, height = 3)
ggsave("../Connectivity/Figures/Figure 3f.jpeg", plot = gg_outbreaks, dpi = 300, width = 3, height = 3)






