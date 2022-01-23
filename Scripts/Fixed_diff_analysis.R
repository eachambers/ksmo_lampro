setwd("~/Box Sync/Lampropeltis project/Writing_gent:tri/SuppMaterials/3_Analyses/5_Fixed_difference_analysis/")

library(adegenet)
library(vcfR)
library(tidyverse)
library(dartR)
library(plotly)
library(directlabels)
library(ggthemes)
library(pheatmap)
library(cowplot)
library(viridis)

theme_set(theme_cowplot())

## The following code takes VCF files with a single SNP per locus and converts to genlight object, adds strata for
## population designations, filters and analyzes data using dartR, removes admixed individuals, and performs a fixed 
## difference analysis. We performed two analyses for this: calculating fixed differences between reference groups
## and outgroups, excluding the L. alterna hybrid (n=21), and a full fixed difference analysis (collapsing, significance
## testing, etc.) with all samples, but once again excluding the L. alterna hybrid (n=92).

##    FILES REQUIRED:
##          FDA_input_files/lampro_all_usnps_noMD_NOALTALLELES.vcf
##          FDA_input_files/FDA_data.txt
##          FDA_input_files/ind_assignments_recode_data.xlsx # use this to delete select inds from genlight object

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
gl_usnps_all <- gl.read.vcf("./FDA_input_files/lampro_all_usnps_noMD_NOALTALLELES.vcf")

# Load data for inds
sampling <- read_tsv("FDA_input_files/FDA_data.txt", col_names = TRUE)



# (2) Add strata to gl object -----------------------------------------------------------

# Add strata to the genlight object defining the different windows
strata(gl_usnps_all) <- sampling %>% 
  select(-ordered_no_long, -ordered_no_vcf)

setPop(gl_usnps_all) <- ~SW_onedeglong
popNames(gl_usnps_all) # should be 9
# table(pop(gl_usnps_all))



# (3) Filter out monomorphs ----------------------------------------

# There are monomorphic loci; verify this by:
gl.report.monomorphs(gl_usnps_all) # 682 monomorphic loci in uSNPs

# Remove monomorphic loci:
gl_usnps_all_new <- gl.filter.monomorphs(gl_usnps_all, verbose=5) # removed 682 monomorphic loci

gl_usnps_all_new_2 <- gl.recalc.metrics(gl_usnps_all_new) # calc metrics for each locus
gl_usnps_all_new_2$other$loc.metrics <- as.data.frame(gl_usnps_all_new_2$other$loc.metrics) # changes "loc.metrics" slot from a list to a df

gl.compliance.check(gl_usnps_all_new_2, verbose=5)

## Take a look at the structure of the data
# pcoa <- gl.pcoa(gl_usnps_all_new_2)
# gl.pcoa.plot(pcoa, gl_usnps_all_new_2)

# Save the gl object
# saveRDS(gl_usnps_all_new_2, file="all_genlight_unfiltered_usnps.RData")

# Load gl object
# gl_usnps_all_new_2 <- readRDS("all_genlight_unfiltered_usnps.RData")



# (4) Remove irrelevant samples for FDAs -------------------------------------------------

# Only reference groups + outgroups (-L. alterna hybrid)
# Easiest way to do this is to alter the recode file:

# Generate recode file
gl.make.recode.ind(gl_usnps_all_new_2, out.recode.file = "FDA_input_files/ind_assignments_refs.csv", outpath = getwd())
gl.make.recode.ind(gl_usnps_all_new_2, out.recode.file = "FDA_input_files/ind_assignments_all.csv", outpath = getwd())

# Go into excel and manually add "Delete" to the second column for samples you want removed;
# You can just copy the second column from the "FDA_input_files/ind_assignments_recode_data.xlsx" file
gl_usnps_refs <- gl.recode.ind(gl_usnps_all_new_2, ind.recode="FDA_input_files/ind_assignments_refs.csv") # 4 pops, 21 individuals
gl_usnps_all <- gl.recode.ind(gl_usnps_all_new_2, ind.recode="FDA_input_files/ind_assignments_all.csv") # 9 pops, 92 individuals



# (5) Calculate fixed diffs and COMPARE -----------------------------------------

##### Analysis (1) is for only reference samples + outgroups
# Set 'populations' to be groups
setPop(gl_usnps_refs) <- ~SW_onedeglong
popNames(gl_usnps_refs) # should be four pop names

# Calculate fixed diffs for ref samples + OGs
D_refs <- gl.fixed.diff(gl_usnps_refs)
# D_refs$fd # look at fixed diffs for reference groups

##### Analysis (2) all samples + OGs
## Set 'populations' to be samples, even though they're admixed
setPop(gl_usnps_all) <- ~sample_ID
popNames(gl_usnps_all)
D_all_inds <- gl.fixed.diff(gl_usnps_all)
# D_all_inds$fd

## Set 'populations' to be groups
setPop(gl_usnps_all) <- ~SW_onedeglong
D_all <- gl.fixed.diff(gl_usnps_all)
# D_all$fd


# (6) AMALGAMATE & REITERATE --------------------------------------------------------------

# Amalgamate populations based on fixed diffs between groups
D_all_2 <- gl.collapse(D_all, verbose=5)
D_all_3 <- gl.collapse(D_all_2, verbose=5)
# D_all_4 <- gl.collapse(D_all_3, verbose=5) # "no further amalgamation of populations"



# (7) TEST --------------------------------------------------------------------

D_all_3 <- gl.fixed.diff(D_all_3$gl, test=TRUE, reps=10000) # all are 0 (significant)
D_all_3$sdfpos
D_all_3$fd
