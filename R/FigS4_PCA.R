library(adegenet)
library(tidyverse)
library(cowplot)
library(vcfR)
library(fields)

theme_set(theme_cowplot())

## The following code generates the PCA figure (Fig. S4), color-coded based on locality

##  Steps in this script:
##    Pull out hybrid individuals from conStruct analysis for color-coding in subsequent PCA plots
##    (1) Import data and convert to genind
##    (2) Add strata defining populations to genind objects
##    (3a) Run PCA on full dataset
##    (3b) Run PCA on dataset excluding elapsoides
##    (3c) Run PCA only on CZ samples from KS and MO
##    (4) Visualize results for all PCAs using ggplot
##    (5) Calculate distances
##    (6) Run correlation test between PC1 loading and geographic distance

##    FILES REQUIRED:
##          1_Bioinformatics/iPyrad_output_files/n93_no_missing_data/all_lampro_pruned.str
##          1_Bioinformatics/iPyrad_output_files/n90_ksmo_alterna/ksmo_alterna_noMD.str (3 elapsoides samples removed)
##          1_Bioinformatics/iPyrad_output_files/n85_ksmo_transect/lampro_ksmo.str (only samples from KS and MO CZ included; n=85)
##          3_Analyses/metadata_n93.txt
##          ../data_files_input_into_scripts/corr_test_data_F11400.txt


# (1) Import data -------------------------------------------------------------

# This will automatically detect a .str file and convert to genind
# Answers to output questions:
# no. genotypes: 93 (all_lampro) or 90 (ksmo_alterna) or 85 (lampro_ksmo)
# no. markers: 16383
# labels: 1
# pop factor: 0
# marker names: 0
# single row: n (there are two rows per individual)
obj <- import2genind("1_Bioinformatics/iPyrad_output_files/n93_no_missing_data/all_lampro_pruned.str")
obj_noelap <- import2genind("1_Bioinformatics/iPyrad_output_files/n90_ksmo_alterna/ksmo_alterna_noMD.str")
obj_ksmo <- import2genind("1_Bioinformatics/iPyrad_output_files/n85_ksmo_transect/lampro_ksmo.str")



# (2) Add strata to genind ----------------------------------------------------

# Add strata for values; make sure this is in the ind order (i.e., how the str file was ordered)
# 3 inds were removed by above conversion; these had high amounts of missing data
# The three samples were: F10981_KS, F10987_KS, and F14200_KS
# The population assignment in this strata file is simply based on locality, where
# KS+OK=gentilis, TX=annulata, FL=elapsoides, and alterna
strata <- read_tsv("3_Analyses/metadata_n93.txt", col_names = TRUE) %>% 
  dplyr::select(sample_ID, strata)

strata_noelap <-
  strata %>% 
  filter(strata!="elapsoides")

strata_ksmo <-
  strata_noelap %>% 
  filter(strata!="alterna")

# Assign strata to strata within the obj genind object
strata(obj) <- strata
strata(obj_noelap) <- strata_noelap
strata(obj_ksmo) <- strata_ksmo

# Can take a look at separate parts of the strata by doing:
# head(strata(obj, ~locality))

# Now, let's take one of these strata columns (population) and assign it as the pop
setPop(obj) <- ~strata
setPop(obj_noelap) <- ~strata
setPop(obj_ksmo) <- ~strata

# Check that it worked
# popNames(obj)
# popNames(obj_noelap)
# popNames(obj_ksmo)


# (3a) Run PCA on full dataset -------------------

sum(is.na(obj$tab)) # 181924
X_full <- scaleGen(obj, NA.method="mean")

# Verify that it worked:
# X_full[1:20, 1:20]

pca1_full <- dudi.pca(X_full, cent=FALSE, scale=FALSE, scannf=FALSE, nf=3)



# (3b) PCA on no elap -----------------------------------------------------

sum(is.na(obj_noelap$tab)) # 159231
X_noelap <- scaleGen(obj_noelap, NA.method="mean")
pca1_noelap <- dudi.pca(X_noelap, cent=FALSE, scale=FALSE, scannf=FALSE, nf=3)



# (3c) Run PCA on only KS-MO samples --------------------------

sum(is.na(obj_ksmo$tab)) # 148383
X_ksmo <- scaleGen(obj_ksmo, NA.method="mean")
pca1_ksmo <- dudi.pca(X_ksmo, cent=FALSE, scale=FALSE, scannf=FALSE, nf=3)



# (4) Build PCA plots ---------------------------------------------------------

# Define color palettes
palette_ksmo <- c("KS"="#b67431", "MO"="#a8a2ca")
palette_noelap <- c("KS"="#b67431", "MO"="#a8a2ca", "alterna"="gray74")
palette <- c("KS"="#b67431", "MO"="#a8a2ca", "alterna"="gray74", "elapsoides"="gray30")

draw_PCA_plot_PC2 <- function(sample_names, strata_dat, other_axis, palette_vals){
  cbind(sample_ID=rownames(sample_names), sample_names, 
        locality=strata_dat, row.names=NULL) %>% 
    ggplot(aes(x=Axis1, y=Axis2, group=locality, color=locality)) +
    geom_vline(xintercept=0, size=1, color="lightgrey") +
    geom_hline(yintercept = 0, size=1, color="lightgrey") +
    geom_point(alpha=0.5, size=5) +
    theme(legend.position="none",
          axis.text = element_text(size=16),
          axis.title = element_text(size=20)) +
    scale_color_manual(values=palette_vals) +
    xlab("PC 1") +
    ylab("PC 2")
}

draw_PCA_plot_PC3 <- function(sample_names, strata_dat, other_axis, palette_vals){
  cbind(sample_ID=rownames(sample_names), sample_names, 
        locality=strata_dat, row.names=NULL) %>% 
    ggplot(aes(x=Axis1, y=Axis3, group=locality, color=locality)) +
    geom_vline(xintercept=0, size=1, color="lightgrey") +
    geom_hline(yintercept = 0, size=1, color="lightgrey") +
    geom_point(alpha=0.5, size=5) +
    theme(legend.position="none",
          axis.text = element_text(size=16),
          axis.title = element_text(size=20)) +
    scale_color_manual(values=palette_vals) +
    xlab("PC 1") +
    ylab("PC 2")
}
    
p_full_PC12 <- draw_PCA_plot_PC2(sample_names=pca1_full$li, strata_dat=strata$strata, palette_vals=palette)
p_full_PC13 <- draw_PCA_plot_PC3(sample_names=pca1_full$li, strata_dat=strata$strata, palette_vals=palette)

p_noelap_PC12 <- draw_PCA_plot_PC2(sample_names=pca1_noelap$li, strata_dat=strata_noelap$strata, palette_vals=palette_noelap)
p_noelap_PC13 <- draw_PCA_plot_PC3(sample_names=pca1_noelap$li, strata_dat=strata_noelap$strata, palette_vals=palette_noelap)

p_ksmo_PC12 <- draw_PCA_plot_PC2(sample_names=pca1_ksmo$li, strata_dat=strata_ksmo$strata, palette_vals=palette_ksmo)
p_ksmo_PC13 <- draw_PCA_plot_PC3(sample_names=pca1_ksmo$li, strata_dat=strata_ksmo$strata, palette_vals=palette_ksmo)

plot_grid(p_ksmo_PC12, p_ksmo_PC13)
ggsave("FigS4EF_KSMO_PCA.pdf", width=11, height=6)

plot_grid(p_noelap_PC12, p_noelap_PC13)
ggsave("FigS4CD_noelap_PCA.pdf", width=11, height=6)

plot_grid(p_full_PC12, p_full_PC13)
ggsave("FigS4AB_all_PCA.pdf", width=11, height=6)



# (5) Calculate distances ----------------------------------------------------

# The following code is also contained within the conStruct analysis script
# Make sure input file is a) long then lat and b) in the same order as the ustr file (ordered_no_vcf)!
ksmo_coords <- read_tsv("3_Analyses/metadata_n93.txt", col_names = TRUE) %>% 
  dplyr::filter(strata != "elapsoides",
                strata != "alterna") %>% 
  dplyr::select(long, lat) %>% 
  as.matrix()

# Calculate geo distances
ksmo_geoDist <- fields::rdist.earth(x1=ksmo_coords, miles=FALSE, R=NULL)

# We just need the column for the farthest west sample (F11400_KS) for the PCA correlation test
# This sample is no. 22 if samples are ordered by vcf
geodists <- as.data.frame(ksmo_geoDist[,22])
write_tsv(geodists, "dist_matrix.txt", col_names = FALSE)

# Now, combine distance data with PCA results
corr_data <-
  pca1_ksmo$li %>% 
  rownames_to_column(var = "sample_ID") %>% 
  cbind(geodists) %>% 
  rename(distance = `ksmo_geoDist[, 22]`) %>% 
  dplyr::select(-Axis2, -Axis3)


# (6) Run correlation test ------------------------------------------------

res <- cor.test(corr_data$Axis1, corr_data$distance,
                method = "pearson")
res

corr_data %>% 
  ggplot(aes(x=distance, y=Axis1)) +
  geom_point() +
  ylab("PCA axis 1") +
  xlab("Distance (km)") +
  geom_smooth(method="lm")
