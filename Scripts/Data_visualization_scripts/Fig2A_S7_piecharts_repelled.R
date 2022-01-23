setwd("~/Box Sync/Lampropeltis project/Writing_gent:tri/SuppMaterials/4_Data_visualization/Scripts")

library(cowplot)
library(tidyverse)
# remotes::install_github("hms-dbmi/repel")
library(repel)

theme_set(theme_cowplot())

##  The following code generates piecharts from the sNMF and conStruct analyses (Figs. 2a and S7).
##  A large portion of this code was developed with Spencer J Fox (https://spncrfx.wordpress.com/)

##    FILES REQUIRED:
##    Results from snmf analysis (ancestry proportions (from Q matrices), lat, long, sample ID):
##          ../data_files_input_into_scripts/lampro_ksmo_k2_run5_data.txt
##    Results from conStruct analysis:
##          ../data_files_input_into_scripts/conStruct_k2_results.txt

##    STRUCTURE OF CODE:
##              (1) Imports data
##              (2) Write functions to build pies
##              (3) Input data into functions
##              (4) Build background plot with samples and lines
##              (5) Build pies and export plots



# (1) Import data --------------------------------------------------------

k2_data_ksmo <- read_tsv("../data_files_input_into_scripts/lampro_ksmo_k2_run5_data.txt", col_names = TRUE)
k2_construct_ksmo <- read_tsv("../data_files_input_into_scripts/conStruct_k2_results.txt", col_names = TRUE)

## Build background map
us <- map_data("state")



# (2) Function to make pies ---------------------------------------------------

# Make pies with data
make_pies <- function(bar_df, alpha = 0.5, legend.position = 'none'){
  bar_df %>% 
    gather(key, value, k1:k2) %>% 
    ggplot() +
      geom_bar(aes(x = "", y = value, fill = key),
             stat = "identity", width = 1, alpha = alpha) +
    coord_polar("y") +
    theme_void() +
    theme(legend.position = legend.position) +
    scale_fill_manual(values = c("#b67431", "#a8a2ca"))
}

# Draw pies
draw_pies <- function(plot, lat, long, height = 1, width = 1) {
  draw_plot(plot = plot, 
            x = long, 
            y = lat,
            height = height,
            width = width, 
            hjust = 0.5, 
            vjust = 0.5)
}


# (3) Input data into functions -----------------------------------------------

# sNMF data (do the same for conStruct data)
df <- k2_data_ksmo %>%
  bind_cols(k2_data_ksmo %>% 
              dplyr::select(x = long, y = lat) %>% 
              mutate(label = 'XX') %>% 
              repel_text(box.padding = .5, force = 4) %>% # there are many params you could set here
              dplyr::select(new_lat = y, new_long = x))

## If there are samples you don't want repelled:
df %>% 
  mutate(to_repel = long>-100) %>% 
  mutate(new_lat = ifelse(to_repel, new_lat, lat),
         new_long = ifelse(to_repel, new_long, long)) -> df


pies_to_add <-  df %>% 
  nest(data = c('k1','k2')) %>% 
  mutate(pies = map(data, make_pies, alpha = 1)) %>% 
  mutate(drawn_pies = pmap(tibble(plot=pies, lat=new_lat, long=new_long), draw_pies, height = .4, width = .4))


# (4) Build background map with samples and lines -----------------------------

us %>% 
  ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group),
               color = 'white', fill = '#e6e6e6', size = 2) +
  coord_fixed(xlim = c(-102.5,-92.75), ylim = c(36.5,40.5), ratio = 1.3) +
  theme_map() +
  theme(legend.position="none") +
  geom_segment(data = df, aes(x = long, y = lat, xend =new_long, yend = new_lat), size = 0.75) +
  geom_point(data = df, aes(long, lat), size = 3) -> kans_map

# For plot with just samples mapped replace with:
# And comment out geom_segment() line
# geom_point(data = df, aes(long, lat), size = 7, stroke=0.1, color="white", fill="black") 


# (5) Build and export plots --------------------------------------------------------------

# Finally, plot the pies you've drawn on the map you've drawn
list(kans_map, pies_to_add$drawn_pies) %>% 
  reduce(.f = `+`) 
ggsave("sNMF_piechart_repelled.pdf", width=11, height=8.5)

# list(kans_map, pies_to_add$drawn_pies) %>% 
#   reduce(.f = `+`)
# ggsave("conStruct_piechart_repelled.pdf", width=11, height=8.5)
