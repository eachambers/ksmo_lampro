library(tidyverse)
library(cowplot)
library(hzar)

theme_set(theme_cowplot())

## The following code generates Fig. S8, the HZAR results.

##    FILES REQUIRED:
##          pure_allele_dict object from Admixture_index_analysis.R script
##          freq output from HZAR_analysis.R script
##          3_Analyses/metadata_ksmo.txt
##          3_Analyses/3_Admixture_index/admix_index_input_files/seq_data.txt

##    STRUCTURE OF CODE:
##              (1) Import data files and tidy them
##              (2) Calculate distances for individual samples (for plotting)
##              (3) Build plot (Fig. S8)


# (1) Load and import data ----------------------------------------------------

load("3_Analyses/3_Admixture_index/pure_allele_dict.rda")
load("3_Analyses/4_HZAR/hzar_best_model.rda")

# metadata file requires lat and long values; please contact author if needed
dat <- read_tsv("3_Analyses/metadata_ksmo.txt", col_names = TRUE) %>% 
  arrange(ordered_no_long)

seq <- read_tsv("3_Analyses/3_Admixture_index/admix_index_input_files/seq_data.txt", col_names = TRUE, col_types = cols(.default = "c"))

# Make ind freq df
seq %>% 
  gather(locus, value, -population, -sample_ID, -read_data) %>% 
  filter(value != 'N') %>% 
  inner_join(pure_allele_dict, by = 'locus') %>% 
  group_by(sample_ID, population) %>% 
  # generate proportions according to sample and locus value 
  # compared to pure dictionary value across all loci
  summarize(frac_gent_ind = sum(value==pure_gent)/n(),
            frac_sys_ind = sum(value==pure_sys)/n()) %>% 
  # remove outgroup species
  filter(population != 'alterna', population !='elapsoides') %>% 
  # append on locality info, calculating avg lat and long per group/population
  left_join(dat %>% 
              group_by(SW_onedeglong) %>% 
              rename(population = SW_onedeglong)) %>% 
  group_by(sample_ID, population, frac_gent_ind, frac_sys_ind, lat, long) %>% 
  # create final table with one row per group and a col with no. inds per group
  summarize(no_samples = n()) %>% 
  arrange(long) -> freq_data_inds



# (2) Calculate distances for individual samples -------------------------

# Calculate distances
distance_inds <- geodist::geodist(freq_data_inds, measure="geodesic")

# Convert into df, convert from m to km
distance_inds <- as.data.frame(distance_inds[,1]) %>% 
  summarize(distance = `distance_inds[, 1]`/1000)

plot_data <-
  cbind(as.data.frame(freq_data_inds %>% 
                        arrange(long)), distance_inds) %>% 
  dplyr::select(sample_ID, frac_gent_ind, frac_sys_ind, distance, no_samples, population)

plot_data <-
plot_data %>% 
  mutate(color = case_when(population == "ks_1" ~ "#fcae60",
                           population == "ks_2" ~ "#89d062",
                           population == "ks_3" ~ "#35b779",
                           population == "ks_4" ~ "#22908c",
                           population == "ks_5" ~ "#443a83",
                           population == "pure_gent" ~ "#b67431",
                           population == "pure_sys" ~ "#a8a2ca"))

# Because HZAR has 0 distance set to average of westernmost group, the ind
# values need to extend below that so they're plotted correctly
# How far is the farthest west sample from the pure gent group average?
plot_data %>% 
  filter(distance==0) %>% 
  summarize(sample_ID) -> low_ind # F11400_KS

# Get lat/long values for the farthest west individual
freq_data_inds %>% 
  filter(sample_ID==low_ind)

# Get lat/long values from HZAR farthest west pop
freq_data_pops %>% 
  filter(population=="pure_gent")

# Request latitude values from author
farwest_coords <- c(-102.01200, XXX)
farwest_avg_coords <- c(-101.12071, XXX)

geodist::geodist(farwest_coords, farwest_avg_coords) # 80.23658 km

# Now, correct the plot_data to go lower than the HZAR 0 value
plot_data %>% 
  mutate(corr_dist = distance-80.23658) -> plot_data



# (3) Build and export plot (Fig. S8) ------------------------------------

# Build base plot with HZAR curve and averages per group as black points
# Xlim needs to go low enough to encompass individual samples' distances
hzar.plot.fzCline(freq$avg$analysis$oDG$data.groups$modelI, pch = 19, xlim = c(-82,675),
                  xlab = "Distance (km)", ylab = "Admixture index")

# Add each individual's data point based on distance
points(plot_data$corr_dist, plot_data$frac_sys_ind, pch = 19, col = alpha(plot_data$color, 0.7))

# Export the plot 6x4
ggsave("FigS8_HZAR.pdf", width = 6, height = 4)

