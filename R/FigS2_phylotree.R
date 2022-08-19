# install_github("YuLab-SMU/ggtree")

# library(devtools) # just to install ggtree from github
library(ggtree)
library(tidyverse)
library(cowplot)
library(ape)
library(RColorBrewer)
library(scales)

theme_set(theme_cowplot())

## The following code generates the phylogenetic tree for the Lampropeltis dataset.
## Code written by E. Anne Chambers

##    FILES REQUIRED:
##          ../3_Analyses/1_RAxML_output_files/lampro_n93_noMD.phy.raxml.support.tree
##          ../data_files_input_into_scripts/tree_coding_data.txt



# Import trees ------------------------------------------------------------

# tree <- read.nexus("../../1_Bioinformatics/iPyrad_output_files/n93_no_missing_data/all_lampro_n93_noMD.nexus")
tree <- read.tree("../../3_Analyses/1_RAxML_output_files/lampro_n93_noMD.phy.raxml.support.tree")

## Check node numbering
ggtree(tree) +
  geom_text2(aes(subset = !isTip, label=node)) +
  geom_tiplab(size=2)

## Root the tree so alterna and elapsoides are outgroups
tree_rooted <- root(tree, node = 179, edgelabel = TRUE)

## Check node numbering again for rooted tree
ggtree(tree_rooted) +
  geom_text2(aes(subset = !isTip, label=node)) +
  geom_tiplab(size=2)

## Check bootstraps
ggtree(tree_rooted) +
  geom_text2(aes(subset = !isTip, label=label)) +
  geom_tiplab(size=2)


## Build a diverging color palette for BS values, blue is high and red is low
x <- seq(0, 1, length.out = 100)
cc <- scales::seq_gradient_pal("#c7544a", "#4365d0", "Lab")(x)



# Branch (edge) colouration -----------------------------------------------

# The first nodes to be labeled are the terminal tips, which means nodes 1-93 are tips.
# You can verify this by doing the following:
tree_edges <- c(tree_rooted$edge[,1], tree_rooted$edge[,2]) # pull out nodes connecting edges
as.data.frame(table(tree_edges)) # count frequency of nodes; tips will only occur once

# Import data for color-coding
dat <- read_tsv("../data_files_input_into_scripts/tree_coding_data.txt", col_names = TRUE)

# Tip nodes are numbered based on the order samples occur in tree$tip.label.
# Get the node numbering for samples:
node_order <-
  as.data.frame(tree_rooted$tip.label) %>% 
  rename(sample_ID = "tree_rooted$tip.label") %>% 
  mutate(node=c(1:93)) %>% 
  left_join(dat, by = "sample_ID") %>% 
  mutate(color = case_when(population == "ks_1" ~ "#fcae60",
                           population == "ks_2" ~ "#89d062",
                           population == "ks_3" ~ "#35b779",
                           population == "ks_4" ~ "#22908c",
                           population == "ks_5" ~ "#443a83",
                           population == "pure_gent" ~ "#b67431",
                           population == "pure_sys" ~ "#a8a2ca",
                           population == "alterna" ~ "gray74",
                           population == "elapsoides" ~ "gray30"))

# Also add remaining (internal nodes) so they show up colored on the tree
# max(tree_rooted$edge[,1]) # max is 184
tree_int <- 
  as.data.frame(c(94:184)) %>% 
  rename(node = "c(94:184)") %>% 
  mutate(color = "black")

# Now, combine all into a single data frame
tree_colors <-
  node_order %>% 
  select(node, color) %>% 
  rbind(., tree_int)



# Build trees ----------------------------------------

# Build tree with support as colored node points (Fig. S2)
p_nodes <-
  ggtree(tree_rooted, size=.4) +
  geom_nodepoint(aes(color=as.numeric(label))) + # size=as.numeric(label)
  scale_colour_gradientn(colors=cc) +
  geom_tiplab(size=2) +
  geom_treescale()

# Color edges according to groups
p_edges_colored <-
  ggtree(tree_rooted) %<+% tree_colors + aes(color=I(color))


### Build a simpler tree with tip circles colored according to locality
# ggtree(tree_rooted, size=.3) %<+%
#   node_order +
#   geom_tippoint(aes(color=population), size=1.5) +
#   scale_color_manual(values=c("ks_1"="#fcae60", "ks_2"="#89d062", "ks_3"="#35b779", "ks_4"="#22908c", "ks_5"="#443a83",
#                               "pure_gent"="#b67431", "pure_sys"="#a8a2ca")) +
#   geom_treescale() +
#   geom_tiplab(size=2)



# Export trees ------------------------------------------------------------

p_nodes
ggsave("lampro_tree_nodes.pdf", width=7.8, height = 8.3)

p_edges_colored
ggsave("lampro_tree_edges_colored.pdf", width=7.8, height = 8.3)

# Need to manually combine these two trees elsewhere
