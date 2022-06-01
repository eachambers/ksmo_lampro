library(devtools)
# install_github("gbradburd/conStruct", ref="covariance_fix")
library(conStruct)
library(tidyverse)
library(fields)

## The following code runs the conStruct analysis for KS-MO transect samples and is modified
## from Marshall et al. (2021) Mol. Phylogenet. Evol. 162:107194.

## The following file:
##    (1) Subsamples ustr file for KS-MO subset and exports file
##    (2) Converts a Structure file (.ustr) to conStruct file (allele freq data)
##    (3) Calculates geographic distances between samples
##    (4) Runs conStruct

##    FILES REQUIRED:
##            ../../../1_Bioinformatics/iPyrad_output_files/all_lampro_pruned.ustr
##      Geographic sampling coordinates (long, lat) in tsv format and in order of VCF file:
##            ../conStruct_input_files/ksmo_coords_MD50.txt



# (1) Create subsampled datasets ----------------------------------------------

# Load ustr file; 210 obs. (105 samples) and 3254 variables
all_data <- read_tsv("../../../1_Bioinformatics/iPyrad_output_files/all_lampro_pruned.ustr", col_names=FALSE)

# remove >50% MD samples (n=12)
MD_50_remove_samples <- c("F10982_KS", "F10983_KS", "F8366_KS", "F9743_KS",
                   "F10984_KS", "F11005_KS", "F10996_KS", "F8504_KS",
                   "F13643_KS", "F14200_KS", "F10981_KS", "F10987_KS")

# Now, retain only KS and MO samples
# n=170 (85 samples)
MD50_ksmo <-
  all_data %>% 
  filter(!`X1` %in% MD_50_remove_samples) %>% 
  filter(grepl("_KS", `X1`) |
         grepl("_MO", `X1`))

# Save the modified file
write_delim(MD50_ksmo, "../conStruct_input_files/n85_ksmo_MD50.ustr", col_names=FALSE)



# (2) Convert ustr to conStruct file ------------------------------------------

# Saves an RData file in working directory
conStruct.data <- structure2conStruct(infile="../conStruct_input_files/n85_ksmo_MD50.ustr",
                                      start.loci=2, # first col that contains data
                                      start.samples=1, # first row that contains samples
                                      onerowperind=FALSE, # two rows per individual
                                      missing.datum=-9, # missing data
                                      outfile="../conStruct_input_files/n85_ksmo_MD50")



# (3) Calculate geographic distances ------------------------------------------

# Make sure input file is a) long then lat and b) in the same order as the ustr file!
# The following file requires both lat and long values; contact author if needed
ksmo_coords <- read_tsv("../conStruct_input_files/ksmo_coords_MD50.txt", col_names=TRUE) # 85 obs.
ksmo_coords <- as.matrix(ksmo_coords)
ksmo_geoDist <- rdist.earth(x1=ksmo_coords, miles=FALSE, R=NULL)



# (4) Run conStruct -----------------------------------------------------------

k2_ksmo <- conStruct(spatial = TRUE,
                        K = 2,
                        freqs = freqs, # be sure to load the data
                        geoDist = md50_geoDist,
                        coords = md50_coords,
                        prefix = "md50",
                        control = setNames(list(0.9),"adapt_delta"),
                        n.chains = 2,
                        n.iter = 12000,
                        make.figs = TRUE,
                        save.files = TRUE)
