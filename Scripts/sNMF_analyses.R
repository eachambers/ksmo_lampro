setwd("~/Box Sync/Lampropeltis project/Writing_gent:tri/SuppMaterials/3_Analyses/2_Population_genetic_structure/Scripts")

library(LEA)
library(cowplot)

##    Written by E. Anne Chambers
##    The following code runs the sNMF analyses. There are two sets of analyses run in sNMF:
##          (1) Only samples from KS-MO contact zone, including alterna (lampro_ksmo_alt) - 50% MD samples n=90
##          (2) Only samples from KS-MO contact zone (lampro_ksmo) - 50% MD samples n=85
##          (3) This code also generates cross-entropy score comparisons (Fig. S5)

##    FILES REQUIRED:
##          vcf files (still have MD samples and are missing commented lines): 
##              ../../1_Bioinformatics/iPyrad_output_files/n90_ksmo_alterna/ksmo_alterna_noMD.vcf # n=90
##              ../../1_Bioinformatics/iPyrad_output_files/n85_ksmo_transect/lampro_ksmo.vcf # n=85

##    STEPS IN CODE:
##          (1) Convert vcf to lfmm file
##          (2) Run sNMF analyses
##          (3) Load sNMF projects (can bypass steps 1 and 2)
##          (4) Examine cross-entropy values (generate Fig. S5)



# (1) Convert vcf to lfmm file ------------------------------------------------

### 1435 invariant sites removed
ksmo_alt_lfmm <- 
  vcf2lfmm(input.file = "../../1_Bioinformatics/iPyrad_output_files/n90_ksmo_alterna/ksmo_alterna_noMD.vcf", 
         output.file = "../sNMF_input_files/lampro_ksmo_alterna.lfmm", 
         force = TRUE)

ksmo_lfmm <- 
  vcf2lfmm(input.file = "../../1_Bioinformatics/iPyrad_output_files/n85_ksmo_transect/lampro_ksmo.vcf", 
           output.file = "../sNMF_input_files/lampro_ksmo.lfmm",
           force = TRUE)



# (2) Run sNMF analyses -------------------------------------------------------

### Output from sNMF analysis:
### For each run, there is a .G file (), a .Q file (Q matrix), and a .snmfClass file ()

## Analysis (1): ks, mo + alterna (n=90)
lampro_snmf_ksmo_alt_5runs = snmf("../sNMF_input_files/lampro_ksmo_alterna.lfmm",
                                iterations=10000, K=1:10, rep=5,
                                entropy=T, CPU=8, ploidy=2) # rep is no. of runs for each value of K

# lampro_snmf_ksmo_alt_5runs <- load.snmfProject("lampro_ksmo_alterna.snmfProject")

## Analysis (2): ks and mo only (n=85)
lampro_snmf_ksmo_5runs = snmf("../sNMF_input_files/lampro_ksmo.lfmm",
                            iterations=10000, K=1:10, rep=5,
                            entropy=T, CPU=8, ploidy=2) # rep is no. of runs for each value of K



# (3) Load sNMF projects --------------------------------------------------

## Bypass steps 1-2 and load projects:
lampro_snmf_ksmo_5runs <- load.snmfProject("../sNMF_input_files/lampro_ksmo.snmfProject")
lampro_snmf_ksmoalt_5runs <- load.snmfProject("../sNMF_input_files/lampro_ksmo_alterna.snmfProject")



# (4) Look at cross-entropy values (Fig. S5) --------------------------------------------------

# figure out which of the runs for each K is the best (minimized)
which.min(cross.entropy(lampro_snmf_ksmoalt_5runs, K = 1)) # do for K 1 through 10

## Take a look at the entropy values for each K
p_alt <- plot(lampro_snmf_ksmoalt_5runs, col = "dodgerblue4", pch = 19, cex = 1.2, ylim=c(0.28, 0.33))
p_ksmo <- plot(lampro_snmf_ksmo_5runs, col = "dodgerblue2", pch = 19, cex = 1.2, ylim=c(.28, 0.33))
