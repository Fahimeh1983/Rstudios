work.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/Patchseq_vs_FACs_cre_analysis/mouse_patchseq_VISp_20181220_collapsed40_cpm/"

library(dplyr)
library(ggplot2)
library(ggExtra)
library(feather)
library(purrr)
library(cowplot)
library(tibble)
library(gridExtra)
library(ggpubr)
options(stringsAsFactors = F)

source(paste0(work.dir,"color_functions.r"))

#load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/de.summary.rda")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/V1_ALM/process_24411/de.summary.rda")
#load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/region.de.summary.rda")
#load("//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/V1_ALM/process_24411/region.de.df.rda")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/V1_ALM/process_24411/region.de.summary.rda")
#Fahimeh: de.summary has these columns: de.num (number of DE genes between each pair of cl), de.lfc(mean log2 of de genes for that pair)
#cl1 and cl2(pair of cls), ???
#Fahimeh: region.de.summary has these columns: de.num(number of DE genes between pair of cls in VISP and ALM regions), 
#de.lfc(mean log2 of de genes expression for that pair)

#anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180312/anno.feather")
anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
#Anno file contains all the cells in VISP and ALM

class_colors <- read.csv(paste0(work.dir,"class_comparison_colors.csv"))

cl_cluster1 <- anno %>%
  filter(cluster_id %in% 1:125) %>%
  select(cl, cluster_id, cluster_label, cluster_color, subclass_id, subclass_label, subclass_color, class_label) %>%
  unique()

cl_cluster2 <- cl_cluster1
names(cl_cluster1) <- c("cl1","cluster_id1","cluster_label1","cluster_color1","subclass_id1","subclass_label1","subclass_color1","class_label1")
names(cl_cluster2) <- c("cl2","cluster_id2","cluster_label2","cluster_color2","subclass_id2","subclass_label2","subclass_color2","class_label2")

pairwise <- de.summary %>%
  mutate(comp_type = "pairwise") %>%
  select(comp_type, cl1, cl2, de.num, de.lfc) %>%
  mutate(cl1 = as.numeric(as.character(cl1)),
         cl2 = as.numeric(as.character(cl2))) %>%
  left_join(cl_cluster1) %>%
  left_join(cl_cluster2)  %>%
  mutate(cluster_region1 = ifelse(grepl("ALM", cluster_label1), "ALM", "none"),
         cluster_region1 = ifelse(grepl("VISp",cluster_label1), "VISp", cluster_region1)) %>%
  mutate(cluster_region2 = ifelse(grepl("ALM", cluster_label2), "ALM", "none"),
         cluster_region2 = ifelse(grepl("VISp",cluster_label2), "VISp", cluster_region2)) %>%
  left_join(class_colors)

region <- region.de.summary %>%
  rownames_to_column("cl1_cl2") %>%
  separate(cl1_cl2, c("cl1","cl2")) %>%
  mutate(cl1 = as.numeric(sub("[A-Z|p]+","",cl1)),
         cl2 = as.numeric(sub("[A-Z|p]+","",cl2))) %>%
  mutate(comp_type = "region") %>%
  select(comp_type, cl1, cl2, de.num, de.lfc) %>%
  left_join(cl_cluster1) %>%
  left_join(cl_cluster2)  %>%
  mutate(cluster_region1 = "VISp",
         cluster_region2 = "ALM") %>%
  left_join(class_colors)

both <- rbind(pairwise, region) 

p <- ggplot() +
  geom_point(data = both,
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color)) +
  scale_color_identity() + 
  theme_bw(20)

ggMarginal(p, groupColour = TRUE, groupFill = TRUE, bw = 0.1)

# Closest match from cross-region comparisons.
p2 <- ggplot() +
  geom_point(data = both %>%
               mutate(class_color = ifelse(comp_type == "pairwise",
                                           "#808080",class_color)),
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color)) +
  scale_color_identity() +
  theme_bw(20)

ggMarginal(p2, groupColour = TRUE, groupFill = TRUE, bw = 0.1)


# Within region Gluta vs Out-of-region Gluta
p3 <- ggplot() +
  geom_point(data = both %>%
               mutate(class_color = ifelse(comp_type == "region", "#808080", class_color)) %>%
               mutate(class_color = ifelse(class_type %in% c("inreg_both_gluta",
                                                             "outreg_both_gluta"),
                                           class_color,
                                           "#808080")),
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color)) +
  scale_color_identity() +
  theme_bw(20)

ggMarginal(p3, groupColour = TRUE, groupFill = TRUE, bw = 0.1)


# Within region Gluta vs Out-of-region Gluta
# Within the same subclass
# Median values
p4_meds <- both %>%
  filter(comp_type == "pairwise") %>%
  filter(class_type %in% c("inreg_both_gluta","outreg_both_gluta")) %>%
  filter(subclass_label1 == subclass_label2) %>%
  group_by(class_type, class_color) %>%
  summarise(lfc_med = median(de.lfc),
            den_med = median(de.num))


p4 <- ggplot() +
  geom_point(data = both %>%
               mutate(class_color = ifelse(comp_type == "region", "#808080", class_color)) %>%
               mutate(class_color = ifelse(class_type %in% c("inreg_both_gluta",
                                                             "outreg_both_gluta"),
                                           class_color,
                                           "#808080")) %>%
               mutate(class_color = ifelse(subclass_label1 == subclass_label2,
                                           class_color,
                                           "#808080")),
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color),
             size = 0.1) +
  geom_hline(data = p4_meds,
             aes(yintercept = lfc_med,
                 color = class_color),
             linetype = "dashed",
             size = 0.1) +
  geom_vline(data = p4_meds,
             aes(xintercept = log10(den_med + 1),
                 color = class_color),
             linetype = "dashed",
             size = 0.1) +
  scale_color_identity() +
  theme_bw(4)

p4 <- ggMarginal(p4, 
           groupColour = TRUE, 
           groupFill = TRUE, 
           xparams = list(bw = 0.1),
           yparams = list(bw = 0.2))

plot(p4)
# Within region Gluta vs Out-of-region Gluta
# Out of the same subclass
p5_meds <- both %>%
  filter(comp_type == "pairwise") %>%
  filter(class_type %in% c("inreg_both_gluta","outreg_both_gluta")) %>%
  filter(subclass_label1 != subclass_label2) %>%
  group_by(class_type, class_color) %>%
  summarise(lfc_med = median(de.lfc),
            den_med = median(de.num))

p5 <- ggplot() +
  geom_point(data = both %>%
               mutate(class_color = ifelse(comp_type == "region", "#808080", class_color)) %>%
               mutate(class_color = ifelse(class_type %in% c("inreg_both_gluta",
                                                             "outreg_both_gluta"),
                                           class_color,
                                           "#808080")) %>%
               mutate(class_color = ifelse(subclass_label1 != subclass_label2,
                                           class_color,
                                           "#808080")) %>%
               mutate(class_color = ifelse(cluster_label1 == "CR Lhx5" | cluster_label2 == "CR Lhx5",
                                           "#808080",
                                           class_color)),
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color),
             size = 0.1) +
  geom_hline(data = p5_meds,
             aes(yintercept = lfc_med,
                 color = class_color),
             linetype = "dashed",
             size = 0.1) +
  geom_vline(data = p5_meds,
             aes(xintercept = log10(den_med + 1),
                 color = class_color),
             linetype = "dashed",
             size = 0.1) +
  scale_color_identity() +
  theme_bw(4)

p5 <- ggMarginal(p5, 
           groupColour = TRUE, 
           groupFill = TRUE, 
           xparams = list(bw = 0.04),
           yparams = list(bw = 0.1))


# Within region GABA vs GABA pairwise comparisons
# Within the same subclass
p6_meds <- both %>%
  filter(class_type %in% c("inreg_both_gaba","outreg_both_gaba")) %>%
  filter(subclass_label1 == subclass_label2) %>%
  group_by(class_type, class_color) %>%
  summarise(lfc_med = median(de.lfc),
            den_med = median(de.num))

p6 <- ggplot() +
  geom_point(data = both %>%
               mutate(class_color = ifelse(class_type %in% c("inreg_both_gaba",
                                                             "outreg_both_gaba"),
                                           class_color,
                                           "#808080")) %>%
               mutate(class_color = ifelse(subclass_label1 == subclass_label2,
                                           class_color,
                                           "#808080")),
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color),
             size = 0.1) +
  geom_hline(data = p6_meds,
             aes(yintercept = lfc_med,
                 color = class_color),
             linetype = "dashed",
             size = 0.1) +
  geom_vline(data = p6_meds,
             aes(xintercept = log10(den_med + 1),
                 color = class_color),
             linetype = "dashed",
             size = 0.1) +
  scale_color_identity() +
  theme_bw(4)

p6 <- ggMarginal(p6, 
           groupColour = TRUE, 
           groupFill = TRUE, 
           xparams = list(bw = 0.04),
           yparams = list(bw = 0.1))
plot(p6)
# Within region GABA vs GABA pairwise comparisons
# Between different subclasses
p7_meds <- both %>%
  filter(class_type %in% c("inreg_both_gaba","outreg_both_gaba")) %>%
  filter(subclass_label1 != subclass_label2) %>%
  group_by(class_type, class_color) %>%
  summarise(lfc_med = median(de.lfc),
            den_med = median(de.num))

p7 <- ggplot() +
  geom_point(data = both %>%
               mutate(class_color = ifelse(class_type %in% c("inreg_both_gaba",
                                                             "outreg_both_gaba"),
                                           class_color,
                                           "#808080")) %>%
               mutate(class_color = ifelse(subclass_label1 != subclass_label2,
                                           class_color,
                                           "#808080")) %>%
               mutate(class_color = ifelse(cluster_label1 == "Meis2 Adamts19" | cluster_label2 == "Meis2 Adamts19",
                                           "#808080",
                                           class_color)),
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color),
             size = 0.1) +
  geom_hline(data = p7_meds,
             aes(yintercept = lfc_med,
                 color = class_color),
             linetype = "dashed",
             size = 0.1) +
  geom_vline(data = p7_meds,
             aes(xintercept = log10(den_med + 1),
                 color = class_color),
             linetype = "dashed",
             size = 0.1) +
  scale_color_identity() +
  theme_bw(4)

p7 <- ggMarginal(p7, 
           groupColour = TRUE, 
           groupFill = TRUE, 
           xparams = list(bw = 0.04),
           yparams = list(bw = 0.1))
plot(p7)
all_plots <- plot_grid(p4, p5, p6, p7,
                       nrow = 1,
                       rel_widths = 1,
                       rel_heights = 1)

ggsave("de_gene_panel.pdf",
       width = 10,
       height = 2.5,
       useDingbats = F)


##########################################################
### Plot cross region comparison #########################
##########################################################

data  = both %>%
  mutate(class_color = ifelse(comp_type == "pairwise",
                              "#808080",class_color))

sp <- ggscatter(data, x = "log10(both$de.num +1)", y = "de.lfc",
                color = "class_color", palette = "jco",
                size = 1, alpha = 0.6)+
  border() 

data  = both %>%
  mutate(class_color = ifelse(comp_type == "pairwise",
                              "#FFFFFF",class_color))

xplot <- ggdensity(data, "log10(both$de.num +1)", fill = "class_color")
yplot <- ggdensity(data, "de.lfc", fill = "class_color")+
  rotate()
yplot <- yplot + clean_theme() 
xplot <- xplot + clean_theme()

ggarrange(xplot, NULL, sp, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)



####################################################################
### VISp cells only, within subclass comparison, SSt, Vip, Pvalb ###
####################################################################


cl_cluster1 <- anno %>%
  filter(cluster_id %in% 1:125, region_label=="VISp", !grepl("ALM", cluster_label)) %>%
  select(cl, cluster_id, cluster_label, cluster_color, subclass_id, subclass_label, subclass_color, class_label) %>%
  unique()

cl_cluster2 <- cl_cluster1
names(cl_cluster1) <- c("cl1","cluster_id1","cluster_label1","cluster_color1","subclass_id1","subclass_label1","subclass_color1","class_label1")
names(cl_cluster2) <- c("cl2","cluster_id2","cluster_label2","cluster_color2","subclass_id2","subclass_label2","subclass_color2","class_label2")

pairwise <- de.summary %>%
  mutate(comp_type = "pairwise") %>%
  select(comp_type, cl1, cl2, de.num, de.lfc) %>%
  mutate(cl1 = as.numeric(as.character(cl1)),
         cl2 = as.numeric(as.character(cl2))) %>%
  left_join(cl_cluster1) %>%
  left_join(cl_cluster2)  %>%
  mutate(cluster_region1 = ifelse(grepl("ALM", cluster_label1), "ALM", "none"),
         cluster_region1 = ifelse(grepl("VISp",cluster_label1), "VISp", cluster_region1)) %>%
  mutate(cluster_region2 = ifelse(grepl("ALM", cluster_label2), "ALM", "none"),
         cluster_region2 = ifelse(grepl("VISp",cluster_label2), "VISp", cluster_region2)) %>%
  left_join(class_colors)


both <- pairwise


p6_meds <- both %>%
  filter(class_type %in% c("inreg_both_gaba")) %>%
  filter(subclass_label1 %in% c("Vip", "Sst", "Pvalb")) %>%
  filter(subclass_label1 == subclass_label2) %>%
  group_by(subclass_label1, class_color) %>%
  summarise(lfc_med = median(de.lfc),
            den_med = median(de.num))

p6_meds[p6_meds$subclass_label1=="Sst","class_color"] <- "Red"
p6_meds[p6_meds$subclass_label1=="Vip","class_color"] <- "Green"
p6_meds[p6_meds$subclass_label1=="Pvalb","class_color"] <- "Blue"

  
p6 <- ggplot() +
  geom_point(data = both %>%
               mutate(class_color = ifelse(class_type %in% c("inreg_both_gaba",
                                                             "outreg_both_gaba"),
                                           class_color,
                                           "#808080")) %>%
               mutate(class_color = case_when(subclass_label1 == subclass_label2 & subclass_label1 == "Sst" ~ "Red",
                                              subclass_label1 == subclass_label2 & subclass_label1 == "Vip" ~ "Green",
                                              subclass_label1 == subclass_label2 & subclass_label1 == "Pvalb" ~ "Blue",
                                            TRUE ~ "#808080")),
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color),
             size = 0.1) +
  geom_hline(data = p6_meds,
             aes(yintercept = lfc_med,
                 color = class_color),
             linetype = "dashed",
             size = 0.1) +
  geom_vline(data = p6_meds,
             aes(xintercept = log10(den_med + 1),
                 color = class_color),
             linetype = "dashed",
             size = 0.1) +
  scale_color_identity() + xlab("log10(n DE genes + 1)") + ylab("Mean (log2 (DE gene expression))")
  theme_bw(4)

p6 <- ggMarginal(p6, 
                 groupColour = TRUE, 
                 groupFill = TRUE, 
                 xparams = list(bw = 0.04),
                 yparams = list(bw = 0.1))
plot(p6)

