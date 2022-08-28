library(quantreg)
library(MCMCpack)
# install_github("GrahamDB/hzar")
library(hzar)
library(geodist)
library(geosphere)
library(tidyverse)
library(cowplot)

theme_set(theme_cowplot())

## The following code generates the HZAR input files.

##    FILES REQUIRED:
##          3_Analyses/3_Admixture_index/admix_index_input_files/seq_data.txt
##          3_Analyses/metadata_ksmo.txt
##          pure_allele_dict object from Admixture_index_analysis.R script


##    STRUCTURE OF CODE:
##              (1) Import data
##              (2) Generate average admix index data
##              (3) Calculate distances between samples
##              (4) Generate HZAR input file



# (1) Import data ---------------------------------------------------------

load("3_Analyses/3_Admixture_index/pure_allele_dict.rda")

# Below requires both lat/long data; contact the author if needed
dat <- read_tsv("3_Analyses/metadata_ksmo.txt", col_names = TRUE) %>% 
  arrange(ordered_no_long) %>% 
  dplyr::select(sample_ID, lat, long, SW_onedeglong)

seq <- read_tsv("3_Analyses/3_Admixture_index/admix_index_input_files/seq_data.txt", col_names = TRUE)



# (2) Generate average admix index data -----------------------------------

seq %>% 
  gather(locus, value, -population, -sample_ID, -read_data) %>% 
  filter(value != 'N') %>% 
  inner_join(pure_allele_dict, by = 'locus') %>% 
  group_by(sample_ID, population) %>% 
  # generate proportions according to sample and locus value 
  # compared to pure dictionary value across all loci
  summarize(frac_gent_ind = sum(value==pure_gent)/n(),
            frac_sys_ind = sum(value==pure_sys)/n()) %>% 
  # append on locality info, calculating avg lat and long per group/population
  left_join(dat %>% 
              group_by(SW_onedeglong) %>% 
              summarize(avg_long = mean(long),
                        avg_lat = mean(lat)) %>% 
              rename(population = SW_onedeglong)) %>% 
  # remove outgroup species
  filter(population != 'alterna', population !='elapsoides') %>% 
  group_by(population, avg_long, avg_lat) %>% 
  # create final table with one row per group and a col with no. inds per group
  summarize(avg_pop_sys = mean(frac_sys_ind),
            avg_pop_gent = mean(frac_gent_ind),
            no_samples = n()) -> freq_data_pops

# Save object to use in generating Fig. S8
save(freq_data_pops, file = "3_Analyses/4_HZAR/freq_data_pops.rda")



# (3) Calculate distances -------------------------------------------------

# Extract lat, long, and population values, making sure that ordering
# is according to longitudinal grouping (this is how geodist will set 0)
coords_pops <-
  freq_data_pops %>% 
  dplyr::select(avg_lat, avg_long) %>% 
  rename(lat = avg_lat,
         long = avg_long) %>% 
  arrange(long)

# Calculate distances
distance_pops <- geodist::geodist(coords_pops, measure="geodesic")

# Convert into df, convert from m to km as required by HZAR
distance_pops <- as.data.frame(distance_pops[,1]) %>% 
  summarize(distance = `distance_pops[, 1]`/1000)



# (4) Create input file -------------------------------------------------------

# Combine group, average admixture index values, distances, and no. samples together:
input_hzar_avg <-
  cbind(as.data.frame(freq_data_pops %>% 
                        arrange(avg_long)), distance_pops) %>% 
  dplyr::select(population, avg_pop_sys, avg_pop_gent, distance, no_samples)

save(input_hzar_avg, file = "3_Analyses/4_HZAR/input_hzar_avg.rda")
