library(cowplot)
library(tidyverse)
library(ggmap)
library(rgdal)
library(precrec)

theme_set(theme_cowplot())

##  The following code generates STRUCTURE-style plots from the sNMF and conStruct analyses (Figs. 2b and S6)

##    FILES REQUIRED:
##    Results from snmf analysis (ancestry proportions (from Q matrices), lat, long, sample ID):
##          4_Data_visualization/data_files_input_into_scripts/lampro_ksmo_alt_k3_run1_data.txt
##          4_Data_visualization/data_files_input_into_scripts/lampro_ksmo_k2_run5_data.txt
##    Results from conStruct analysis:
##          4_Data_visualization/data_files_input_into_scripts/conStruct_k2_results.txt

##    STRUCTURE OF CODE:
##              (1) Imports data
##              (2) Tidies data
##              (3) Builds STRUCTURE-style plots
##              (4) Exports plots



# (1) Import data -------------------------------------------------------------

# Upload data and define variables

### Visualization (1): sNMF KS-MO and alterna
k3_data_ksmoalt <- read_tsv("4_Data_visualization/data_files_input_into_scripts/lampro_ksmo_alt_k3_run1_data.txt", col_names = TRUE)

### Visualization (2): sNMF KS-MO contact zone samples
k2_data_ksmo <- read_tsv("4_Data_visualization/data_files_input_into_scripts/lampro_ksmo_k2_run5_data.txt", col_names = TRUE)

### Visualization (3): conStruct KS-MO contact zone samples
k2_construct_ksmo <- read_tsv("4_Data_visualization/data_files_input_into_scripts/conStruct_k2_results.txt", col_names = TRUE)



# (2) Tidy data ---------------------------------------------------------------

tidy_structure_data <- function(data){
  data %>% 
    gather(key="cluster", value="proportion", -ordered_no_long, -ordered_no_vcf, -sample_ID, -long, -lat)
}

### Visualization (1): KS-MO and alterna
k3_data_ksmoalt <- tidy_structure_data(k3_data_ksmoalt)

### Visualization (2): only KS-MO contact zone samples
k2_data_ksmo <- tidy_structure_data(k2_data_ksmo)

### conStruct visualization: ks-mo
k2_construct_ksmo <- tidy_structure_data(k2_construct_ksmo)



# (3) Build the STRUCTURE-style plots ----------------------------------------------------------

structure_plot <- function(data, colors, order) {
  data %>% 
    ggplot(aes(x=order, y=proportion, fill=cluster)) +
    geom_bar(stat="identity") + 
    panel_border() + 
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    scale_fill_manual(values=colors) +
    theme_cowplot() %+replace% theme(axis.line=element_line(colour="black"),
                                     axis.text.x=element_blank(),
                                     axis.text.y=element_blank(),
                                     axis.ticks=element_blank(),
                                     axis.title.x=element_blank(),
                                     axis.title.y=element_blank(),
                                     legend.position="none",
                                     panel.border = element_rect(fill=NA, colour="black", linetype="solid",size=1.5),
                                     strip.text.y = element_text(size=30, face="bold"),
                                     strip.background = element_rect(colour="white", fill="white"),
                                     panel.spacing=unit(-0.1, "lines"))
}


### Define colors
#318c84 # green (alterna)
#a8a2ca # purple (syspila)
#b67431 # brown (gentilis)

k2_cols = c("#b67431", "#a8a2ca")
k2_cols_construct = c("#a8a2ca", "#b67431")
k3_cols = c("#b67431", "#a8a2ca", 
            "#318c84")

## Build plots using the above function

### Visualization (1) ks-mo, alterna
p_k3_ksmoalt_long <- structure_plot(k3_data_ksmoalt, k3_cols, k3_data_ksmoalt$ordered_no_long)

### Visualization (2): ks-mo CZ only
p_k2_ksmo_long <- structure_plot(k2_data_ksmo, k2_cols, k2_data_ksmo$ordered_no_long)

### conStruct visualization: ks-mo
p_k2_construct_ksmo <- structure_plot(k2_construct_ksmo, k2_cols_construct, k2_construct_ksmo$ordered_no_long)


# (4) Export plots --------------------------------------------------------

p_k3_ksmoalt_long
ggsave("FigS6_sNMF_K3_ksmoalt_orderedlong.pdf", width=8.6, height=7.856)

p_k2_ksmo_long
ggsave("Fig2b_sNMF_K2_ksmo_orderedlong.pdf", width=8.6, height=7.856)

p_k2_construct_ksmo
ggsave("Fig2b_conStruct_K2_ksmo_orderedlong.pdf", width=8.6, height=7.856)
