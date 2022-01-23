setwd("~/Box Sync/Lampropeltis project/Writing_gent:tri/SuppMaterials/4_Data_visualization/Scripts")

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

## The following code takes results from Fixed_diff_analysis.R script, visualizes results using a PCoA (Fig. 4), 
## and generates heat maps based on fixed differences (Figs. S9 and S10).

##    FILES REQUIRED:
##          D_all (pre-collapsed groups; results from Fixed_diff_analysis.R)
##          D_all_3 (post-collapsed groups; results from Fixed_diff_analysis.R)
##          sampling (results from Fixed_diff_analysis.R)
##          ../data_input_files_scripts/FDA_data.txt

##    STRUCTURE OF CODE:
##              (1) Run PCoA for pre- and post-collapsed groups
##              (2) Build PCoA plots
##              (3) Export PCoA plots (Fig. 4)
##              (4) Build heat maps
##              (5) Export heat maps (Figs. S9 and S10)


# Remove hybrid individual for visualizations
sampling_nohyb <-
  sampling %>% 
  filter(!sample_ID=="alt17w4f_TX")



# (1) Run PCoA ----------------------------------------------------------------

# Post-collapsed
pca <- gl.pcoa(D_all_3)

# Pre-collapsed
pca_pre <- gl.pcoa(D_all)



# (2) Build PCoA plots -------------------------------------------------------

colors_pre <- c("ks_1"="#dd8666", "ks_2"="#fcae60", "ks_3"="#e9d207", "ks_4"="#fcae60", "ks_5"="#89d062", 
                "ks_6"="#35b779", "ks_7"="#22908c", "ks_8"="#443a83",
                "pure_gent"="#b67431", "pure_sys"="#a8a2ca", "alterna"="gray74", "elapsoides"="gray30")

colors_collapsed <- c("ks_1"="#22908c", "ks_2"="#22908c", "ks_3"="#22908c", "ks_4"="#22908c", "ks_5"="#22908c", 
                      "ks_6"="#22908c", "ks_7"="#22908c", "ks_8"="#22908c",
                      "pure_gent"="#22908c", "pure_sys"="#22908c", 
                      "alterna"="gray74", "elapsoides"="gray30")

pcoa_collapsed <-
  as_tibble(cbind(pca$scores, sampling_nohyb)) %>% 
  ggplot(aes(x=as.numeric(PC1), y=as.numeric(PC2), group=SW_onedeglong, color=SW_onedeglong)) +
  geom_vline(xintercept=0, size=1, color="lightgrey") +
  geom_hline(yintercept = 0, size=1, color="lightgrey") +
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values=colors_collapsed) +
  ylab("PC2") +
  xlab("PC1")


pcoa_pre <-
  as_tibble(cbind(pca_pre$scores, sampling_nohyb)) %>% 
  ggplot(aes(x=as.numeric(PC1), y=as.numeric(PC2), group=SW_onedeglong, color=SW_onedeglong)) +
  geom_vline(xintercept=0, size=1, color="lightgrey") +
  geom_hline(yintercept = 0, size=1, color="lightgrey") +
  geom_point(alpha=0.8, size=3) +
  scale_color_manual(values=colors_pre) +
  ylab("PC2") +
  xlab("PC1")


# (3) Export PCoA plots (Fig. 4) -------------------------------------------------------

pcoa_pre
ggsave("Fig4a_PCoA_FDA_pre-collapse.pdf", width=5, height=4)

pcoa_collapsed
ggsave("Fig4b_PCoA_FDA_collapsed.pdf", width=5, height=4)



# (4) Build heat maps ---------------------------------------------------------------

order <- 
  sampling_nohyb %>% 
  dplyr::select(sample_ID, ordered_no_long) %>% 
  rename(samples = sample_ID,
         ordered_no_long_new = ordered_no_long)

# Looking at fixed diffs among each individual
dat <-
  as_tibble(cbind(D_all_inds$fd, sampling_nohyb)) %>%
  pivot_longer(-c(sample_ID, SW_onedeglong, ordered_no_long, ordered_no_vcf), 
               names_to = "samples", values_to = "fixed_diffs")

# Build heat map (Fig. S9)
hm_inds <-
  left_join(dat, order, by="samples") %>% 
  mutate(sample_ID = fct_reorder(as.factor(sample_ID), ordered_no_long),
         samples = fct_reorder(as.factor(samples), ordered_no_long_new)) %>% 
  ggplot(aes(x=samples, y=sample_ID, fill=fixed_diffs)) +
  geom_tile() +
  scale_fill_viridis(option="mako", direction = -1) +
  theme(panel.border = element_rect(fill=NA, colour="darkgrey", linetype="solid", size=2),
        axis.title=element_blank(),
        axis.text=element_text(size=8),
        axis.text.x=element_text(angle=90, hjust=1),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        axis.line = element_line(color="darkgrey"),
        axis.ticks = element_line(color="darkgrey")) +
  coord_fixed()


### Now, let's build the heat map for populations
# For ordering samples in the heat map:
pop_order <- c("alterna", "pure_gent", "ks_1", "ks_2", "ks_3", "ks_4", "ks_5", "pure_sys", "elapsoides")

dat <- as.data.frame(cbind(D_all$fd, comp_pop=rownames(D_all$fd))) %>% 
  pivot_longer(-comp_pop, names_to = "population", values_to = "fixed_diffs")

levels(dat$comp_pop) = pop_order
levels(dat$population) = pop_order

# Build heat map for fixed diffs among groups (Fig. S10)
hm_groups <-
  dat %>% 
  ggplot(aes(x=comp_pop, y=population, fill=as.numeric(fixed_diffs))) +
  geom_tile() +
  scale_fill_viridis(option="mako", direction = -1) +
  theme(panel.border = element_rect(fill=NA, colour="darkgrey", linetype="solid", size=2),
        axis.title=element_blank(),
        axis.text.y=element_text(size=12),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        axis.line = element_line(color="darkgrey"),
        axis.ticks = element_line(color="darkgrey")) +
  coord_fixed()



# (5) Export heat maps --------------------------------------------------------

hm_inds
ggsave("FigS9_heatmap_inds.pdf", width=10, height=10)

hm_groups
ggsave("FigS10_heatmap_groups.pdf", width=10, height=10)
