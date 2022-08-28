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
##    (5) Runs cross-validation analysis comparing spatial and non-spatial models

##    FILES REQUIRED:
##            1_Bioinformatics/iPyrad_output_files/all_lampro_pruned.ustr
##      Geographic sampling coordinates (long, lat) in tsv format and in order of VCF file:
##            3_Analyses/metadata_ksmo.txt



# (1) Create subsampled datasets ----------------------------------------------

# Load ustr file; 210 obs. (105 samples) and 3254 variables
all_data <- read_tsv("1_Bioinformatics/iPyrad_output_files/n105_all_samples/all_lampro_pruned.ustr", col_names=FALSE)

# Remove >50% MD samples (n=12)
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
write_delim(MD50_ksmo, "3_Analyses/2_Population_genetic_structure/conStruct_input_files/n85_ksmo_MD50.ustr", col_names=FALSE)



# (2) Convert ustr to conStruct file ------------------------------------------

# Saves an RData file in working directory
conStruct.data <- conStruct::structure2conStruct(infile="3_Analyses/2_Population_genetic_structure/conStruct_input_files/n85_ksmo_MD50.ustr",
                                      start.loci=2, # first col that contains data
                                      start.samples=1, # first row that contains samples
                                      onerowperind=FALSE, # two rows per individual
                                      missing.datum=-9, # missing data
                                      outfile="3_Analyses/2_Population_genetic_structure/conStruct_input_files/n85_ksmo_MD50")


# (3) Calculate geographic distances ------------------------------------------

# Make sure input file is a) long then lat and b) in the same order as the ustr file!
# The following requires lat/long values; please contact author if needed
ksmo_coords <- read_tsv("3_Analyses/metadata_ksmo.txt", col_names = TRUE) %>% 
  dplyr::select(long, lat)

ksmo_coords <- as.matrix(ksmo_coords)
ksmo_geoDist <- fields::rdist.earth(x1=ksmo_coords, miles=FALSE, R=NULL)


# (4) Run conStruct -----------------------------------------------------------

# Load frequency data if haven't run above code:
load("3_Analyses/2_Population_genetic_structure/conStruct_input_files/n85_ksmo_MD50.RData")

k2_ksmo <- conStruct(spatial = TRUE,
                     K = 2,
                     freqs = freqs, # be sure to load the data
                     geoDist = ksmo_geoDist,
                     coords = ksmo_coords,
                     prefix = "md50",
                     control = setNames(list(0.9),"adapt_delta"),
                     n.chains = 2,
                     n.iter = 12000,
                     make.figs = TRUE,
                     save.files = TRUE)

# Save results
# save(k2_ksmo, file = "3_Analyses/2_Population_genetic_structure/conStruct_input_files/md50_conStruct.results.Robj")
# load("3_Analyses/2_Population_genetic_structure/conStruct_input_files/md50_conStruct.results.Robj")

### Export results for plotting purposes
ksmo_coords <- read_tsv("3_Analyses/metadata_ksmo.txt", col_names = TRUE) %>% 
  dplyr::select(-SW_onedeglong)

results <- conStruct.results$chain_1$MAP$admix.proportions %>% 
  as.data.frame() %>% 
  rename(k1 = V1,
         k2 = V2)
  
dat <- cbind(ksmo_coords, results)

write_tsv(dat, "4_Data_visualization/data_files_input_into_scripts/conStruct_k2_results.txt")


# (5) Cross-validation analysis -------------------------------------------

# Run x-validation analysis
my.xvals <- x.validation(train.prop = 0.8,
                         n.reps = 8,
                         K = 2,
                         freqs = freqs,
                         data.partitions = NULL,
                         geoDist = md50_ksmo_geoDist,
                         coords = md50_ksmo_coords,
                         prefix = "md50_ksmo",
                         control = setNames(list(0.9),"adapt_delta"),
                         n.iter = 5000,
                         make.figs = TRUE,
                         save.files = FALSE,
                         parallel = TRUE,
                         n.nodes = 4)

# Read in results of xvalid analysis from text files
sp.results <- as.matrix(
  read.table("md50_ksmo_sp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

nsp.results <- as.matrix(
  read.table("md50_ksmo_nsp_xval_results.txt",
             header = TRUE,
             stringsAsFactors = FALSE)
)

# First, get the 95% confidence intervals for the spatial and nonspatial
#	models over values of K (mean +/- 1.96 the standard error)
sp.CIs <- apply(sp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})
nsp.CIs <- apply(nsp.results,1,function(x){mean(x) + c(-1.96,1.96) * sd(x)/length(x)})

# Then, plot cross-validation results with 8 replicates
par(mfrow = c(1,3))

plot(rowMeans(sp.results),
     pch=19,col="dodgerblue3",
     ylab = "",xlab = "",
     ylim=range(sp.results,nsp.results),
     main = "")

points(rowMeans(nsp.results),col="darkorange2",pch = 19)

