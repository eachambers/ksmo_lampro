library(quantreg)
library(MCMCpack)
library(hzar)
library(geodist)
library(geosphere)
library(tidyverse)

## The following code generates the HZAR input files.

##    FILES REQUIRED:
##          ../3_Admixture_index/admix_index_input_files/seq_data.txt
##          ../../4_Data_visualization/data_files_input_into_scripts/admixture_data.txt
##          pure_allele_dict object from Admixture_index_analysis.R script


##    STRUCTURE OF CODE:
##              (1) Import data
##              (2) Generate average admix index data
##              (3) Calculate distances between samples
##              (4) Generate HZAR input file



# (1) Import data ---------------------------------------------------------

# Below requires both lat/long data; contact the author if needed
dist_data <- read_tsv("admixture_data.txt", col_names = TRUE) %>% 
  arrange(ordered_no_long)

seq <- read_tsv("../3_Admixture_index/admix_index_input_files/seq_data.txt", col_names = TRUE)



# (2) Generate average admix index data -----------------------------------

seq %>% 
  gather(locus, value, -population, -sample_ID, -read_data) %>% 
  filter(value != 'N') %>% 
  inner_join(pure_allele_dict, by = 'locus') %>% 
  group_by(sample_ID, population) %>% 
  summarize(frac_gent_ind = sum(value==pure_gent)/n(),
            frac_sys_ind = sum(value==pure_sys)/n()) %>% 
  filter(population != 'alterna', population !='elapsoides') %>% 
  left_join(dist_data %>% 
              group_by(SW_onedeglong) %>% 
              summarize(avg_long = mean(long),
                        avg_lat = mean(lat)) %>% 
              rename(population = SW_onedeglong), by = 'population') %>% 
  group_by(population, avg_long, avg_lat) %>% 
  summarize(avg_pop_sys = mean(frac_sys_ind),
            avg_pop_gent = mean(frac_gent_ind),
            no_samples = n()) -> freq_data_pops



# (3) Calculate distances -------------------------------------------------

# Calculate distances between populations
coords_pops <-
  freq_data_pops %>% 
  dplyr::select(avg_lat, avg_long) %>% 
  rename(lat = avg_lat,
         long = avg_long)

distance_pops <- geodist(coords_pops, measure="geodesic")
distance_pops <- as.data.frame(distance_pops[,1]) %>% 
  summarize(distance = `distance_pops[, 1]`/1000)



# (4) Create input file -------------------------------------------------------

# Combine group, average admixture index values, distances, and no. samples together:
input_hzar_avg <-
  cbind(as.data.frame(freq_data_pops), distance_pops) %>% 
  dplyr::select(population, avg_pop_sys, avg_pop_gent, distance, no_samples)
