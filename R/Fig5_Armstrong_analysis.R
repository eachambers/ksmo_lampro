library(adegenet)
library(vcfR)
library(tidyverse)
library(ggthemes)
library(viridis)
library(cowplot)

theme_set(theme_cowplot())

## The following code build histograms (Fig. 5) for reference and transect samples based on canonical axis 1 from
## Armstrong et al. (2001).

##    FILES REQUIRED:
##          ../data_files_input_into_scripts/armstrong_data.txt

##    STRUCTURE OF CODE:
##              (1) Import data
##              (2) Generate plots for (i) reference samples and (ii) transect samples
##              (3) Export plots



# (1) Import data ---------------------------------------------------------

dat <- read_tsv("../data_files_input_into_scripts/armstrong_data.txt", col_names = TRUE)



# (2) Generate plots ------------------------------------------------------

# green (elapsoides) = #6b7d74
# blue (triangulum) = #a9bdd2

# Build reference sample plot
p_ref <-
  dat %>%
  filter(population=="ref_elapsoides" |
           population=="ref_syspila") %>% 
  ggplot(aes(CAN1)) +
  geom_histogram(aes(fill=population, group=population, y=stat(density)*2), binwidth = 2) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1)) +
  scale_x_continuous(expand=c(0,0), limits=c(-11,7)) +
  scale_fill_manual(values=c("ref_elapsoides"="#6b7d74", "ref_syspila"="#a9bdd2")) +
  theme(axis.line = element_line(color="grey"),
        axis.ticks = element_line(color="grey"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20, face="bold")) +
  ylab("Frequency of individuals") +
  xlab("Canonical axis 1")


# Build transect sample plot
p_transect <-
  dat %>%
  filter(population=="ky_elapsoides" |
           population=="ky_syspila") %>%
  ggplot(aes(CAN1)) +
  geom_histogram(aes(fill=population, group=population, y=stat(density)*2), binwidth = 2) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1)) +
  scale_x_continuous(expand=c(0,0), limits=c(-11,7)) +
  scale_fill_manual(values=c("ky_elapsoides"="#6b7d74", "ky_syspila"="#a9bdd2")) +
  theme(axis.line = element_line(color="grey"),
        axis.ticks = element_line(color="grey"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20, face="bold")) +
  ylab("Frequency of individuals") +
  xlab("Canonical axis 1")



# (3) Export plots --------------------------------------------------------

p_ref
ggsave("Fig5b_Armstrong_ref.pdf", width=11, height=5)

p_transect
ggsave("Fig5c_Armstrong_transect.pdf", width=11, height=5)
