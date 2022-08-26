library(tidyverse)
library(dartR)
library(SNPRelate)
library(cowplot)
library(viridis)
library(harrietr)

theme_set(theme_cowplot())

## The following code takes results from Fixed_diff_analysis.R script, visualizes results using a PCoA (Fig. 4), 
## and generates heat maps based on fixed differences (Figs. S9 and S10).

##    FILES REQUIRED:
##          4_Data_visualization/data_files_input_into_scripts/D_all.rda (pre-collapsed groups; results from Fixed_diff_analysis.R)
##          4_Data_visualization/data_files_input_into_scripts/D_all_3.rda (post-collapsed groups; results from Fixed_diff_analysis.R)
##          4_Data_visualization/data_files_input_into_scripts/D_all_inds.rda (each ind assigned a pop; results from Fixed_diff_analysis.R)
##          3_Analyses/metadata_n93.txt

##    STRUCTURE OF CODE:
##              (1) Run PCoA for pre- and post-collapsed groups
##              (2) Build PCoA plots
##              (3) Export PCoA plots (Fig. 4)
##              (4) Build heat maps
##              (5) Export heat maps (Figs. S9 and S10)



# Import and load data ----------------------------------------------------

# Import metadata and remove hybrid individual for visualizations
sampling_nohyb <- read_tsv("3_Analyses/metadata_n93.txt", col_names = TRUE) %>% 
  filter(!sample_ID=="alt17w4f_TX")

# Load results from fixed diff analysis
load("4_Data_visualization/data_files_input_into_scripts/D_all.rda")
load("4_Data_visualization/data_files_input_into_scripts/D_all_3.rda")
load("4_Data_visualization/data_files_input_into_scripts/D_all_inds.rda")


# (1) Run PCoA ----------------------------------------------------------------

# Post-collapsed
pca <- dartR::gl.pcoa(D_all_3)

# Pre-collapsed
pca_pre <- dartR::gl.pcoa(D_all)


# (2) Build PCoA plots -------------------------------------------------------

colors_pre <- c("ks_1"="#fcae60", "ks_2"="#89d062", "ks_3"="#35b779", "ks_4"="#22908c", "ks_5"="#443a83",
                "pure_gent"="#b67431", "pure_sys"="#a8a2ca", "alterna"="gray74", "elapsoides"="gray30")

colors_collapsed <- c("ks_1"="gray50", "ks_2"="gray50", "ks_3"="gray50", "ks_4"="gray50", "ks_5"="gray50", 
                      "pure_gent"="gray50", "pure_sys"="gray50", 
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

# Ordering for sample_ID_comp
order <- 
  sampling_nohyb %>% 
  dplyr::select(sample_ID, ordered_no_long) %>% 
  rename(sample_ID_comp = sample_ID,
         ordered_no_long_new = ordered_no_long)

# Tidy and process data, joining FDs with sample metadata
dat <- as.data.frame(as.matrix(D_all_inds$fd)) %>% 
  rownames_to_column("sample_ID") %>% 
  left_join(sampling_nohyb %>% 
              dplyr::select(sample_ID, ordered_no_long, SW_onedeglong)) %>% 
  as_tibble() %>% 
  pivot_longer(-c(sample_ID, SW_onedeglong, ordered_no_long), 
               names_to = "sample_ID_comp", values_to = "fixed_diffs") %>% 
  left_join(order, by="sample_ID_comp") %>% 
  mutate(sample_ID = fct_reorder(as.factor(sample_ID), ordered_no_long),
         sample_ID_comp = fct_reorder(as.factor(sample_ID_comp), ordered_no_long_new))

### Build individual-level heat map (Fig. S9)
hm_inds <-
  dat %>% 
  ggplot(aes(x=sample_ID, y=sample_ID_comp, fill=fixed_diffs)) +
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

dat_pops <- as.data.frame(as.matrix(D_all$fd)) %>% 
  rownames_to_column("SW_onedeglong") %>% 
  as_tibble() %>% 
  pivot_longer(-SW_onedeglong, 
               names_to = "comp_pop", values_to = "fixed_diffs")

dat_pops$comp_pop <- factor(dat_pops$comp_pop, levels = pop_order)
dat_pops$SW_onedeglong = factor(dat_pops$SW_onedeglong, levels = pop_order)

## Build heat map for fixed diffs among groups (Fig. S10)
hm_groups <-
  dat_pops %>% 
  ggplot(aes(x=as.factor(comp_pop), y=as.factor(SW_onedeglong), fill=as.numeric(fixed_diffs))) +
  geom_tile() +
  scale_fill_viridis(option="mako", direction = -1) +
  theme(panel.border = element_rect(fill=NA, colour="darkgrey", linetype="solid", size=2),
        axis.title=element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(angle = 90),
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
