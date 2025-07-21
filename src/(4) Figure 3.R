# Description

# Author(s): Reilly O'Connor
# Version: 2025-07-20

# Load Pkgs
library(tidyverse)

# load data
source("../Connectivity/src/(1) Nutrients.R")
source("../Connectivity/src/(2) Disease.R")
source("../Connectivity/src/(3) Invasive Species.R")


gg_fig_3 <- plot_grid(gg_nutrient_grid, gg_disease_grid, gg_invasive_grid,
                      nrow = 3, align = "hv")

gg_fig_3

ggsave("../Connectivity/Figures/Figure 3.jpeg", plot = gg_fig_3, dpi = 300, width = 10, height = 12)






