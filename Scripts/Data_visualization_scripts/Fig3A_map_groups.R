setwd("~/Box Sync/Lampropeltis project/Writing_gent:tri/SuppMaterials/4_Data_visualization/Scripts")

library(tidyverse)
library(cowplot)

theme_set(theme_cowplot())

##  The following code generates the map for samples in KS-MO with data points colorized by group assignment  (Fig. 3A).

##    FILES REQUIRED:
##    Sample IDs, coordinates, and group assignment:
##          ../data_files_input_into_scripts/admixture_data.txt

##    STRUCTURE OF CODE:
##              (1) Import data
##              (2) Build background map
##              (3) Add samples
##              (4) Export plot


# (1) Import data ---------------------------------------------------------

dat <- read_tsv("../data_files_input_into_scripts/admixture_data.txt", col_names = TRUE)



# (2) Build background map ------------------------------------------------

## Build background map
us <- map_data("state")

colors <- c("ks_1"="#fcae60", "ks_2"="#89d062", "ks_3"="#35b779", "ks_4"="#22908c", "ks_5"="#443a83",
            "pure_gent"="#b67431", "pure_sys"="#a8a2ca")



# (3) Add samples ---------------------------------------------------------

p_base <-
  us %>% 
  ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group),
               color = 'white', fill = '#e6e6e6', size = 2) +
  coord_fixed(xlim = c(-102.5,-92.75), ylim = c(36.5,40.5), ratio = 1.3) +
  theme_map()

p_full <-
  p_base +
  geom_point(data = dat, aes(long, lat, group=SW_onedeglong, color=SW_onedeglong), size = 4, alpha=.8) +
  scale_color_manual(values=colors)



# (4) Export plot ---------------------------------------------------------

p_full
ggsave("Fig3A_sampling_map.pdf", width=11, height=8.5)
