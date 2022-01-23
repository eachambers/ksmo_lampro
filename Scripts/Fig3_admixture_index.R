setwd("~/Box Sync/Lampropeltis project/Writing_gent:tri/SuppMaterials/4_Data_visualization/Scripts")

library(tidyverse)
library(cowplot)

theme_set(theme_cowplot())

## The following code takes results from Admixture_index_analysis.R script, illustrates admixture index values 
## across loci and individuals in each group (Fig. 3).

##    FILES REQUIRED:
##    freq_data from Admixture_index_analysis.R script
##    seq from Admixture_index_analysis.R script

##    STRUCTURE OF CODE:
##              (1) Plot per locus admixture index values (Fig. 3C)
##              (2) Get per individual averages across all loci
##              (3) Build histograms with per ind data (Fig. 3B)
##              (4) Plot global average line (Fig. 3D)


# (1) Plot per locus admixture index values (Fig. 3C) ---------------------------------------------

freq_data %>% 
  left_join(sampling %>% 
              group_by(SW_onedeglong) %>% 
              summarize(long = mean(long)) %>% 
              rename(population = SW_onedeglong), by = 'population') %>% 
  mutate(population = factor(population, levels = c('pure_gent', 'ks_1', 'ks_2', 'ks_3', 'ks_4', 'ks_5', 'pure_sys'))) %>% 
  ggplot(aes(long, frac_sys, group=locus)) + geom_line(alpha = .3) +
  xlab("Longitude") +
  ylab("Admixture index") +
  theme(axis.title = element_text(size=20, face="bold"),
        axis.text = element_text(size=16),
        legend.position = "none")

ggsave("Fig3C_freq_data_lines.pdf", width=11, height=6)

# write_tsv(freq_data, "freq_data_per_locus.txt")



# (2) Get per individual averages across all loci -----------------------------

seq %>% 
  gather(locus, value, -population, -sample_ID, -read_data) %>% 
  filter(value != 'N') %>% 
  inner_join(pure_allele_dict, by = 'locus') %>% 
  group_by(sample_ID, population) %>% 
  summarize(frac_gent_ind = sum(value==pure_gent)/n(),
            frac_sys_ind = sum(value==pure_sys)/n()) %>% 
  filter(population != 'alterna', population !='elapsoides') %>% 
  left_join(sampling %>% 
              group_by(SW_onedeglong) %>% 
              summarize(avg_long = mean(long)) %>% 
              rename(population = SW_onedeglong), by = 'population') %>% 
  left_join(sampling %>% 
              select(sample_ID, long), by = "sample_ID") -> freq_data_inds

# write_csv(freq_data_inds, "freq_data_inds.csv")



# (3) Build histograms with per ind data (Fig. 3B) --------------------------------------------------------

colors <- c("ks_1"="#fcae60", "ks_2"="#89d062", "ks_3"="#35b779", "ks_4"="#22908c", "ks_5"="#443a83",
              "pure_gent"="#b67431", "pure_sys"="#a8a2ca")

freq_data_inds %>%
  mutate(population = factor(population, levels = c('pure_gent', 'ks_1', 'ks_2', 'ks_3', 'ks_4', 'ks_5', 'pure_sys'))) %>% 
  ggplot(aes(frac_sys_ind)) +
  geom_histogram(aes(fill=population, group=population, y=stat(density)*0.1), binwidth = 0.1) +
  scale_y_continuous(expand = c(0,0), limits=c(0,1.05)) +
  scale_x_continuous(expand=c(0,0), limits=c(-.1,1.1), breaks=c(0,0.5,1)) +
  theme(panel.border = element_rect(fill=NA, colour="grey", linetype="solid",size=1.5),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        axis.line = element_line(color="grey"),
        axis.ticks = element_line(color="grey"),
        panel.grid.major.y = element_line(color="gray92", linetype="dashed"),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20, face="bold")) +
  ylab("Frequency") +
  xlab("Admixture index") +
  scale_fill_manual(values = c("ks_1"="#fcae60", "ks_2"="#89d062", "ks_3"="#35b779", "ks_4"="#22908c", "ks_5"="#443a83",
                               "pure_gent"="#b67431", "pure_sys"="#a8a2ca")) +
  facet_grid(~population)

ggsave("Fig3B_freq_inds.pdf", width=11, height=5)



 # (4) Plot global average line (Fig. 3D) -------------------------------------------------------

p_avg_smoothed <- 
freq_data %>% 
  left_join(sampling %>% 
              group_by(SW_onedeglong) %>% 
              summarize(long = mean(long)) %>% 
              rename(population = SW_onedeglong), by = 'population') %>% 
  mutate(population = factor(population, levels = c('pure_gent', 'ks_1', 'ks_2', 'ks_3', 'ks_4', 'ks_5', 'pure_sys'))) %>%
  group_by(population) %>% 
  mutate(avg_pop_sys = mean(frac_sys)) %>% 
  gather(average_calc, avg_freq, -population, -locus, -frac_gent, -frac_sys, -long) %>% 
  ggplot(aes(long, avg_freq)) + 
  geom_smooth(method="auto", se=TRUE, level=0.95, color="black")

p_avg_smoothed +
  geom_point(data=freq_data_inds, aes(x=long, y=frac_sys_ind, group=population, color=population), size=3, alpha=.8) +
  scale_color_manual(values=colors) +
  xlab("Longitude") +
  ylab("Admixture index") +
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=20, face="bold"),
        legend.position = "none")

ggsave("Fig3D_average_freq_data_smoothed.pdf", width=11, height=6)

