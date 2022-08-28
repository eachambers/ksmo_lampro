library(tidyverse)
library(SNPRelate)
library(dartR)

theme_set(theme_cowplot())

## The following code takes VCF files with a single SNP per locus and converts to genlight object, adds strata for
## population designations, filters and analyzes data using dartR, removes admixed individuals, and performs a fixed 
## difference analysis. We performed two analyses for this: calculating fixed differences between reference groups
## and outgroups, excluding the L. alterna hybrid (n=21), and a full fixed difference analysis (collapsing, significance
## testing, etc.) with all samples, but once again excluding the L. alterna hybrid (n=92).

##    FILES REQUIRED:
##          3_Analyses/5_Fixed_difference_analysis//FDA_input_files/lampro_all_usnps_noMD_NOALTALLELES.vcf
##          3_Analyses/5_Fixed_difference_analysis//FDA_input_files/ind_assignments_recode_data.xlsx # use this to delete select inds from genlight object
##          3_Analyses/metadata_n93.txt

##    STRUCTURE OF CODE:
##              (1) Read in file and convert vcf to genlight
##              (2) Add strata to gl object for populations
##              (3) Filter data using dartR
##              (4) Remove irrelevant samples
##              (5) Calculate fixed differences and COMPARE
##              (6) AMALGAMATE & REITERATE
##              (7) TEST


# (1) Read in file and convert ------------------------------------------------------------

# File should have 93 individuals and 3173 locus labels
gl_usnps_all <- dartR::gl.read.vcf("3_Analyses/5_Fixed_difference_analysis/FDA_input_files/lampro_all_usnps_noMD_NOALTALLELES.vcf")

# Import data for inds
sampling <- read_tsv("3_Analyses/metadata_n93.txt", col_names = TRUE)



# (2) Add strata to gl object -----------------------------------------------------------

# Add strata to the genlight object defining the different windows
strata(gl_usnps_all) <- sampling %>% 
  select(sample_ID, SW_onedeglong)

setPop(gl_usnps_all) <- ~SW_onedeglong
popNames(gl_usnps_all) # should be 9

# List no. inds in each group by:
# table(pop(gl_usnps_all))



# (3) Filter out monomorphs ----------------------------------------

# There are monomorphic loci; verify this by:
dartR::gl.report.monomorphs(gl_usnps_all) # 682 monomorphic loci in uSNPs

# Remove monomorphic loci:
gl_usnps_all_new <- dartR::gl.filter.monomorphs(gl_usnps_all, verbose=5) # removed 682 monomorphic loci

gl_usnps_all_new_2 <- dartR::gl.recalc.metrics(gl_usnps_all_new) # calc metrics for each locus
gl_usnps_all_new_2$other$loc.metrics <- as.data.frame(gl_usnps_all_new_2$other$loc.metrics) # changes "loc.metrics" slot from a list to a df

dartR::gl.compliance.check(gl_usnps_all_new_2, verbose=5)

## Take a look at the structure of the data
# pcoa <- dartR::gl.pcoa(gl_usnps_all_new_2)
# dartR::gl.pcoa.plot(pcoa, gl_usnps_all_new_2)

# Save the gl object
# saveRDS(gl_usnps_all_new_2, file="3_Analyses/5_Fixed_difference_analysis/all_genlight_unfiltered_usnps.RData")

# Load gl object
gl_usnps_all_new_2 <- readRDS("3_Analyses/5_Fixed_difference_analysis/all_genlight_unfiltered_usnps.RData")



# (4) Remove irrelevant samples for FDAs -------------------------------------------------

# Only reference groups + outgroups (-L. alterna hybrid)
# Easiest way to do this is to alter the recode file:

# Generate recode files (these are already included in supp mats)
# dartR::gl.make.recode.ind(gl_usnps_all_new_2, out.recode.file = "3_Analyses/5_Fixed_difference_analysis/FDA_input_files/ind_assignments_refs.csv", outpath = getwd())
# dartR::gl.make.recode.ind(gl_usnps_all_new_2, out.recode.file = "3_Analyses/5_Fixed_difference_analysis/FDA_input_files/ind_assignments_all.csv", outpath = getwd())

# Go into excel and manually add "Delete" to the second column for samples you want removed;
# You can just copy the second column from the "FDA_input_files/ind_assignments_recode_data.xlsx" file
gl_usnps_refs <- dartR::gl.recode.ind(gl_usnps_all_new_2, ind.recode="3_Analyses/5_Fixed_difference_analysis/FDA_input_files/ind_assignments_refs.csv") # new: 4 pops, 21 individuals
gl_usnps_all <- dartR::gl.recode.ind(gl_usnps_all_new_2, ind.recode="3_Analyses/5_Fixed_difference_analysis/FDA_input_files/ind_assignments_all.csv") # new: 9 pops, 93/92??? individuals



# (5) Calculate fixed diffs and COMPARE -----------------------------------------

##### Analysis (1) is for only reference samples + outgroups
# Set 'populations' to be groups
setPop(gl_usnps_refs) <- ~SW_onedeglong
popNames(gl_usnps_refs) # should be four pop names

# Calculate fixed diffs for ref samples + OGs
D_refs <- dartR::gl.fixed.diff(gl_usnps_refs)
# D_refs$fd # look at fixed diffs for reference groups

##### Analysis (2) all samples + OGs (-L. alterna hybrid)
## Set 'populations' to be groups
setPop(gl_usnps_all) <- ~SW_onedeglong
D_all <- dartR::gl.fixed.diff(gl_usnps_all)
# D_all$fd # results in Table S4

save(D_all, file = "4_Data_visualization/data_files_input_into_scripts/D_all.rda")

##### Analysis (3) all samples + OGs with no groups assigned (not used in manuscript)
## Set 'populations' to be samples, even though they're admixed
setPop(gl_usnps_all) <- ~sample_ID
popNames(gl_usnps_all) # should be 92 (no alterna hybrid)
D_all_inds <- dartR::gl.fixed.diff(gl_usnps_all)
# D_all_inds$fd

save(D_all_inds, file = "4_Data_visualization/data_files_input_into_scripts/D_all_inds.rda")


# (6) AMALGAMATE & REITERATE --------------------------------------------------------------

# Amalgamate populations based on fixed diffs between groups
D_all_2 <- dartR::gl.collapse(D_all, verbose=5)
D_all_3 <- dartR::gl.collapse(D_all_2, verbose=5)
# D_all_4 <- dartR::gl.collapse(D_all_3, verbose=5) # "no further amalgamation of populations"
# D_all_3$fd # results in Table S5

# (7) TEST --------------------------------------------------------------------

D_all_3 <- dartR::gl.fixed.diff(D_all_3$gl, test=TRUE, reps=10000) # all are 0 (significant)
D_all_3$sdfpos
D_all_3$fd

save(D_all_3, file = "4_Data_visualization/data_files_input_into_scripts/D_all_3.rda")
