setwd("~/Box Sync/Lampropeltis project/Writing_gent:tri/SuppMaterials/3_Analyses/3_Admixture_index")

library(tidyverse)

## The following code calculates the admixture index values for samples.

##    FILES REQUIRED:
##    Sequence data (each site is a column) including sample IDs and population assignment:
##          admix_index_input_files/seq_data.txt # made using PGDSpider and the lampro_n93_noMD_usnps.vcf file
##    Sample information (mostly important to have longitude data):
##          ../5_Fixed_difference_analysis/FDA_input_files/FDA_data.txt

##    STRUCTURE OF CODE:
##              (1) Import data
##              (2) Calc diagnostic diffs between reference groups
##              (3) Calc per locus admixture index values



# (1) Import data -------------------------------------------------------------

# We used PGDspider to convert 1_Bioinformatics/iPyrad_output_files/n93_no_missing_data/lampro_n93_noMD_usnps.vcf to a fasta file; 
# this was then converted to a data frame

sampling <- read_tsv("../5_Fixed_difference_analysis/FDA_input_files/FDA_data.txt", col_names = TRUE)

seq <- read_tsv("admix_index_input_files/seq_data.txt", col_names = TRUE, col_types = cols(.default = "c"))



# (2) Calculate diagnostic diffs between reference groups  -------------------------

# Set a threshold value and determine diagnostics differences between reference groups
threshold <- 0.8 # change to 1 to see fixed diffs (100% inds different); should be only 10 loci

seq %>%
  filter(population=="pure_gent" |
           population=="pure_sys") %>%                          # 28 rows, 3,247 vars.
  gather(locus, value, -population, -sample_ID, -read_data) %>% # 91,084 rows
  filter(value != 'N') %>%                                      # 84,486 rows
  dplyr::count(population, locus, value) %>%                    # 8,287 rows
  group_by(population, locus) %>%
  mutate(frac=n/sum(n)) %>%
  ungroup() %>%
  filter(frac>=threshold) %>%
  dplyr::select(population, locus, value) %>%
  spread(population, value) %>%
  filter(pure_gent!=pure_sys) -> pure_allele_dict               # 44 rows



# (3) Per locus allele frequencies ----------------------------------------

seq %>% 
  gather(locus, value, -population, -sample_ID, -read_data) %>% 
  filter(value != 'N') %>% 
  inner_join(pure_allele_dict, by = 'locus') %>% 
  group_by(population, locus) %>% 
  summarize(frac_gent = sum(value==pure_gent)/n(),
            frac_sys = sum(value==pure_sys)/n()) %>% 
  filter(population != 'alterna', population !='elapsoides') -> freq_data # 44 loci * 7 groups = 308 rows
