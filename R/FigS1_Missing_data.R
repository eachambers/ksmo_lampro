library(devtools)
library(tidyverse)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())

## The following code generates Fig. S1, which illustrates the amount of missing data
## and number of loci in assembly for each sample. Amounts of missing data were calculated
## using the PAUP missdata function.

##    FILES REQUIRED:
##          ../../2_Data_processing/data_files_input_into_scripts/missing_data_snps_n105.txt

# (1) Import data -------------------------------------------------------------

miss_all <- read_tsv("../../2_Data_processing/data_files_input_into_scripts/missing_data_snps_n105.txt")



# (2) Build base plot with MD -------------------------------------------------------------

# This will build the base (missing data) plot with samples (x-axis) ordered
# by increasing % missing data
p_miss <-
  miss_all %>% 
  filter(!percent_missing>=50) %>% 
  arrange(-desc(percent_missing)) %>% 
  mutate(sample_ID = factor(sample_ID, levels=sample_ID)) %>% 
  ggplot(aes(x=sample_ID, y=percent_missing, fill=percent_missing)) +
    geom_bar(stat="identity") +
    scale_fill_distiller(palette="Spectral", breaks=c(0,25,50,75,100)) +
    theme_cowplot() %+replace% theme(strip.text.x = element_text(size=15, face="bold"),
                                     axis.text.x = element_text(angle=90, size=6),
                                     strip.background = element_rect(colour="white", fill="white"),
                                     #panel.border = element_rect(fill=NA, colour="black", linetype="solid",size=1.5),
                                     axis.title.x=element_text(size=15, face="bold"),
                                     axis.title.y=element_text(size=15, face="bold", angle=90)) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits=c(0,50), breaks=seq(0,50, by=10)) +
    labs(y="Missing data (%)") +
    labs(x="Sample")


# (3) Add loci in assembly as data points with righthand y-axis -------------------------------------

### Multiply second y-axis so it scales up to 50; there are 3230 loci total so 3230*0.015=~100
p_miss <-
  p_miss +
  geom_point(aes(y=loci_in_assembly*0.015), size=1, color="black") +
  scale_y_continuous(sec.axis = sec_axis(~./0.015), expand = c(0,0), limits=c(0,50), breaks=seq(0,50, by=10))



# (4) Export plot ------------------------------------------------------------

p_miss
ggsave("FigS1_missingdata_lociassembly.pdf", width=12, height = 8)
