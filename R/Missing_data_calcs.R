library(devtools)
library(tidyverse)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())

## The following code calculates average amounts of missing data and produces a list of samples with
## missing data above a user-defined threshold. Amounts of missing data were calculated using the PAUP*
## missdata function.

##    FILES REQUIRED:
##          ../data_files_input_into_scripts/missing_data_snps_n105.txt



# (1) Import data -------------------------------------------------------------

miss_all <- read_tsv("../data_files_input_into_scripts/missing_data_snps_n105.txt")



# (2) Find samples with missing data >50% ----------------------------

miss_all %>% 
  filter(percent_missing>=50) # 12 samples

# Calculate average percent missing data BEFORE removing samples with >50% MD
miss_all %>%
  summarize(avg_missing = mean(percent_missing)) # 15.3

# Calculate average, min, and max percent missing data AFTER removing samples with >50% MD
miss_all %>%
  filter(!percent_missing>=50) %>% 
  summarize(avg_missing = mean(percent_missing), # 6.44
            min_missing = min(percent_missing), # 1.52
            max_missing = max(percent_missing)) # 45.0



