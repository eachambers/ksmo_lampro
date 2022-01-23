# Contact zones in the American Milksnake
This repository contains scripts for analyses and data visualization from Chambers et al. (2022) on resolving a species boundary in the *Lampropeltis triangulum* complex. It is organized with scripts for [data processing and analysis](https://github.com/eachambers/ksmo_lampro/tree/main/Scripts/Data_analysis_scripts) and the scripts required for [data visualization](https://github.com/eachambers/ksmo_lampro/tree/main/Scripts/Data_visualization_scripts). All raw data and input files can be found on Dryad [here](LINK_HERE).

*Please note: for scripts requiring sampling coordinates, latitude values are not provided because of concerns related to the commercial exploitation of these animals; I can provide latitude values to qualified investigators upon request for the purposes of re-analysis and reproducibility.*

**Subset coding:**

Most analyses were performed on subsets of the data (in all cases, samples with >50% missing data were removed). These are as follows:
* **Full dataset:** 93 samples (KS-MO transect samples (n=85), *L. alterna* (n=5), *L. elapsoides* (n=3))
* **KS-MO + alterna dataset:** (KS-MO transect samples (n=85), *L. alterna* (n=5))
* **KS-MO transect dataset:** (KS-MO transect samples (n=85))

## Data processing
* [Missing data calculations script](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_analysis_scripts/Missing_data_calcs.R)
    1. Takes PAUP* output and calculates missing data per sample and on average
* [Missing data and loci in assembly visualization](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_visualization_scripts/FigS1_Missing_data.R)
    1. Plots missing data and loci in assembly per sample (**Fig. S1**)

## Phylogenetic tree
* [Data visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_visualization_scripts/FigS2_phylotree.R)
    1. Builds tree figure with node support values colorized according to bootstrap support
    2. Colorizes terminal edges according to group assignment (**Fig. S2**)

## Population Genetics Analyses

### PCA and correlation test
* [Analysis & visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_visualization_scripts/FigS4_PCA.R)
    1. Runs a PCA on three subsets of the data: all samples, KS-MO + alterna, KS-MO transect only
    2. Visualizes PCA for PC1, PC2, and PC3 for each subset **Fig. S4**

### sNMF analysis
* [Analysis script](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_analysis_scripts/sNMF_analyses.R)
    1. Runs sNMF on all three data subsets

### conStruct analysis
* [Analysis script](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_analysis_scripts/conStruct_analysis.R)
    1. Runs conStruct on KS-MO transect samples
* [Data visualization script (**Fig SX**)](LINK_HERE)

### Structure-based data visualization scripts
* [Piecharts, repelled (**Figs. 2a, S7**)](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_visualization_scripts/Fig2A_S7_piecharts_repelled.R)
* [STRUCTURE-style plots (**Figs. 2b, S6**)](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_visualization_scripts/Fig2B_S6_popgen.R)

## Other Analyses

### Admixture index analysis
* [Analysis script](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_analysis_scripts/Admixture_index_analysis.R)
* [Data visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_visualization_scripts/Fig3_admixture_index.R)
    1. Build histograms of admixture index values per group (**Fig. 3b**)
    2. Plots per-locus admixture index values (**Fig. 3c**)
    3. Plots average (across all loci) admixture index values over longitude (**Fig. 3d**)

### Fixed difference analysis
* [Analysis script](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_analysis_scripts/Fixed_diff_analysis.R)
* [Data visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_visualization_scripts/Fig4_S9_S10_Fixed_diff_analysis.R)
    1. Builds PCoA plots for pre- and post-collapsed groups (**Fig. 4**)
    2. Builds heat map with raw fixed difference values among individuals (**Fig. S9**)
    3. Builds heat map with raw fixed difference values among groups (**Fig. S10**)

### HZAR cline fitting
* [Generate HZAR input file from admixture index values](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_analysis_scripts/Generate_HZAR_input_file.R)
* [Run HZAR](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_analysis_scripts/HZAR_analysis.R)

### Morphological re-assessment of Armstrong et al. (2001)
* [Analysis and data visualization script](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_visualization_scripts/Fig5_Armstrong_analysis.R)
    1. Takes in canonical axis 1 data from Armstrong et al. (2001)
    2. Builds histograms for reference and transect samples (**Fig. 5**)

### Map figure
* [Samples into group assignments on map (**Fig. 3a**)](https://github.com/eachambers/ksmo_lampro/blob/main/Scripts/Data_visualization_scripts/Fig3A_map_groups.R)
