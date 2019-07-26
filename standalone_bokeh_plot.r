library(dplyr)
library(feather)
library(rbokeh)

.libPaths(c("/allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Olivia/3.5", .libPaths()))

all_anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/V1_ALM/process_24411/de.summary.rda")

cluster_anno <- all_anno %>%
  select(cl, cluster_label, cluster_color) %>%
  unique()

class_anno <- all_anno %>%
  select(cl, class_id, class_label, class_color) %>%
  unique() %>%
  arrange(class_id)

class_anno2 <- class_anno %>%
  select(-cl) %>%
  unique()

class_options <- class_anno2$class_label
names(class_options) <- class_anno2$class_label

de.summary <- de.summary %>%
  mutate(cl1 = as.numeric(cl1),
         cl2 = as.numeric(cl2)) %>%
  left_join(cluster_anno, by = c("cl1" = "cl")) %>%
  rename(cl1_color = cluster_color) %>%
  left_join(cluster_anno, by = c("cl2" = "cl")) %>%
  rename(cl2_color = cluster_color) %>%
  left_join(class_anno, by = c("cl1" = "cl")) %>%
  rename(cl1_class = class_label) %>%
  left_join(class_anno, by = c("cl2" = "cl")) %>%
  rename(cl2_class = class_label)

figure(width = 800, height = 800,
       xlim = c(1, 4),
       ylim = c(1.5, 10)) %>%
  ly_wedge(data = de.summary,
           x = log10(de.num),
           y = de.lfc,
           start_angle = pi/2,
           end_angle = 3*pi/2,
           color = cl1_color,
           radius = 0.01,
           hover = list("Cluster 1" = cl1_label,
                        "Cluster 2" = cl2_label,
                        "N DEGenes" = de.num,
                        "Mean DE LFC" = round(de.lfc,2),
                        "Correlation" = round(cl.cor,2))) %>%
  ly_wedge(data = de.summary,
           x = log10(de.num),
           y = de.lfc,
           start_angle = -pi/2,
           end_angle = pi/2,
           color = cl2_color,
           radius = 0.01,
           hover = list("Cluster 1" = cl1_label,
                        "Cluster 2" = cl2_label,
                        "N DEGenes" = de.num,
                        "Mean DE LFC" = round(de.lfc,2),
                        "Correlation" = round(cl.cor,2)))


###################################
work.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/Patchseq_vs_FACs_cre_analysis/mouse_patchseq_VISp_20181220_collapsed40_cpm/"
class_colors <- read.csv(paste0(work.dir,"class_comparison_colors.csv"))

cl_cluster1 <- all_anno %>%
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


data = both %>%
  filter(subclass_label1 == subclass_label2) %>%
  #mutate(class_color = ifelse(class_type %in% c("inreg_both_gaba",
  #                                              "outreg_both_gaba"),
  #                            class_color,
  #                            "#000000")) %>%
  mutate(class_color = case_when(subclass_label1 == subclass_label2 & subclass_label1 == "Sst" ~ "Red",
                                 subclass_label1 == subclass_label2 & subclass_label1 == "Vip" ~ "Green",
                                 subclass_label1 == subclass_label2 & subclass_label1 == "Pvalb" ~ "Blue"))


figure(width = 800, height = 800,
       xlim = c(1, 4),
       ylim = c(1.5, 10),
       xlab = "log10 (n DE genes + 1)",
       ylab = "Mean (log2 (DE gene expression))") %>%
  ly_wedge(data = data,
           x = log10(de.num),
           y = de.lfc,
           start_angle = pi/2,
           end_angle = 3*pi/2,
           color = class_color,
           radius = 0.01,
           hover = list("Cluster 1" = cluster_label1,
                        "Cluster 2" = cluster_label2,
                        "N DEGenes" = de.num,
                        "Mean DE LFC" = round(de.lfc,2))) %>%
  ly_wedge(data = data,
           x = log10(de.num),
           y = de.lfc,
           start_angle = -pi/2,
           end_angle = pi/2,
           color = class_color,
           radius = 0.01,
           hover = list("Cluster 1" = cluster_label1,
                        "Cluster 2" = cluster_label2,
                        "N DEGenes" = de.num,
                        "Mean DE LFC" = round(de.lfc,2)))


data = both %>%
  filter(subclass_label1 == subclass_label2)  %>%
  filter (subclass_label1 %in% c("Sst", "Vip", "Pvalb"))


figure(width = 800, height = 800,
       xlim = c(1, 4),
       ylim = c(1.5, 10),
       xlab = "log10 (n DE genes + 1)",
       ylab = "Mean (Log2( fold change (cluster median expression)))") %>%
  ly_wedge(data = data,
           x = log10(de.num),
           y = de.lfc,
           start_angle = pi/2,
           end_angle = 3*pi/2,
           color = cluster_color1,
           radius = 0.01,
           hover = list("Cluster 1" = cluster_label1,
                        "Cluster 2" = cluster_label2,
                        "N DEGenes" = de.num,
                        "Mean DE LFC" = round(de.lfc,2))) %>%
  ly_wedge(data = data,
           x = log10(de.num),
           y = de.lfc,
           start_angle = -pi/2,
           end_angle = pi/2,
           color = cluster_color2,
           radius = 0.01,
           hover = list("Cluster 1" = cluster_label1,
                        "Cluster 2" = cluster_label2,
                        "N DEGenes" = de.num,
                        "Mean DE LFC" = round(de.lfc,2)))

