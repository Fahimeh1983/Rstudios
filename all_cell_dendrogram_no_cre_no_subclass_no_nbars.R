library(dendextend)
library(ggplot2)
library(dplyr)
library(pvclust)
library(feather)
library(scrattch.hicat)
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Rcodes/layer_and_violin_functions.R")

options(stringsAsFactors = F)

######################################################################################################
### Lockdown data-set ################################################################################
######################################################################################################

work.dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/"

facs.anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
FACS_dend <- readRDS("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/dend.RData")

#Removing ALM cells
facs.anno <- facs.anno[facs.anno$region_label == "VISp",]
#Removing lowQ cells
facs.anno <- facs.anno[facs.anno$cluster_label %in% labels(FACS_dend),]
dim(facs.anno)

REF <-  as.data.frame(unique(facs.anno[,c("cluster_label", "cluster_id", "cluster_color", 
                                          "dendcluster_id", "dendcluster_label", "dendcluster_color")]))

REF$dendcluster_label == REF$cluster_label
REF$dendcluster_color == REF$cluster_color

REF <- REF %>% arrange(dendcluster_id)
sum(cbind.data.frame(REF[,c("dendcluster_id", "dendcluster_label")], labels(FACS_dend))[,2] == 
      cbind.data.frame(REF[,c("dendcluster_id", "dendcluster_label")], labels(FACS_dend))[,3])

REF_layer_label <-  as.data.frame(unique(facs.anno[,c("layer_label", "layer_id", "layer_color")]))

######################################################################################################
### Lockdown data-set ################################################################################
######################################################################################################
res.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"
patchseq_dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/"
facs_dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520"

bp.collapse.th = 40
bp.name.add = NULL
if (!is.null(bp.collapse.th)) {
  bp.name.add = paste0(".with.bp.", bp.collapse.th)
}

patchseq.anno <- Assign_clusterlabel_clusterid_for_patchseq(patchseq_dir = patchseq_dir, 
                                                            dend = FACS_dend, 
                                                            ref_cluster_label_color = REF)

patchseq.anno <- as.data.frame(patchseq.anno)
rownames(patchseq.anno) <- patchseq.anno$sample_id

temp1 <- read_feather(paste0(patchseq_dir, "/anno.feather"))
temp1 <- as.data.frame(temp1)
rownames(temp1) <- temp1$sample_id

sum(patchseq.anno$topLeaf_label == temp1$topLeaf_label) == dim(patchseq.anno)[1]
rm(temp1)
unique(patchseq.anno$dendcluster_label) %in% labels(FACS_dend)

temp2 <- cbind.data.frame(as.data.frame(unique(patchseq.anno[,c("dendcluster_id", "dendcluster_label")])) %>% arrange(dendcluster_id), 
                          labels(FACS_dend)[labels(FACS_dend) %in% unique(patchseq.anno$cluster_label)])
temp2[,2] == temp2[,3]

Locked_sample_id <- read.csv(paste0(work.dir, "/locked_sampleids.csv"))
Locked_sample_id <- Locked_sample_id$x
patchseq.anno <- as.data.frame(patchseq.anno)

patchseq.anno <- patchseq.anno[patchseq.anno$sample_id %in% Locked_sample_id,]
core_cells <- patchseq.anno[patchseq.anno$Tree_call_label == "Core","sample_id"]
I1_cells <- patchseq.anno[patchseq.anno$Tree_call_label == "I1", "sample_id"]
I2_cells <- patchseq.anno[patchseq.anno$Tree_call_label == "I2", "sample_id"]
I3_cells <- patchseq.anno[patchseq.anno$Tree_call_label == "I3", "sample_id"]

dim(patchseq.anno)

######################################################################################################
### Subsetting GABA: #################################################################################
######################################################################################################

GABAanno <- patchseq.anno[patchseq.anno$subclass_label %in% c("Vip", "Sst", "Sncg", "Lamp5", "Pvalb", "Serpinf1"),]
dim(GABAanno)
GABAanno %>% group_by(Tree_call_label) %>% summarise(n())

######################################################################################################
### correcting layer info: ###########################################################################
######################################################################################################
#Here it will only take cells that are in this subclasses and in specific layers
dim(patchseq.anno)
dim(GABAanno)
GABAanno_Mola <- Modify_layer_label_GABAcells(GABAanno)
dim(GABAanno_Mola)  
study.cells <- setdiff(rownames(GABAanno), rownames(GABAanno_Mola))
unique(patchseq.anno[study.cells, "subclass_label"])
unique(patchseq.anno[study.cells, "structure_label"])

######################################################################################################
### Choosing Core and I1 cells: ######################################################################
######################################################################################################

GABAanno_Mola <- GABAanno_Mola%>% filter(Tree_call_label %in% c("Core", "I1"))
dim(GABAanno_Mola)

######################################################################################################
### select cres for layer comparison: ################################################################
######################################################################################################

unique(GABAanno_Mola$cre_label)

common_cres <- unique(GABAanno_Mola$cre_label)[unique(GABAanno_Mola$cre_label) %in% unique(facs.anno$cre_label)]

as.data.frame(facs.anno %>% 
                filter(cre_label %in% common_cres) %>%
                filter(layer_label %in% c("L1", "L2/3", "L4", "L5", "L6", "L6b")) %>% 
                filter(subclass_label %in% c("Sst", "Vip", "Sncg", "Pvalb", "Serpinf1", "Lamp5")) %>%
                group_by(cre_label) %>% 
                summarise(freq = n()) %>% 
                arrange(freq))

as.data.frame(GABAanno_Mola %>% 
                filter(cre_label %in% common_cres) %>% 
                filter(subclass_label %in% c("Sst", "Vip", "Sncg", "Pvalb", "Serpinf1", "Lamp5")) %>%
                group_by(cre_label) %>% 
                summarise(freq = n()) %>% 
                arrange(freq))

#All the facs cre are included in patchseq, but some of the patchseq cres are in some layers
#that facs does not have a clear info for that. So we keeo the facs cres that are in common 
#with patchseq and we have clear info about their layer label

keep.cres <- facs.anno %>% 
  filter(cre_label %in% common_cres) %>%
  filter(layer_label %in% c("L1", "L2/3", "L4", "L5", "L6", "L6b")) %>% 
  filter(subclass_label %in% c("Sst", "Vip", "Sncg", "Pvalb", "Serpinf1", "Lamp5")) %>%
  select(cre_label) %>% 
  unique()

facs.layer.figure.data <- facs.anno %>% 
  filter(cre_label %in% keep.cres$cre_label) %>%
  filter(subclass_label %in% c("Sst", "Vip", "Sncg", "Pvalb", "Serpinf1", "Lamp5")) %>%
  filter(layer_label %in% c("L1", "L2/3", "L4", "L5", "L6", "L6b"))

patchseq.layer.figure.data <- GABAanno_Mola %>% 
  filter(cre_label %in% keep.cres$cre_label)

facs.fig1.cells <- facs.layer.figure.data$sample_id
patchseq.fig1.cells <- patchseq.layer.figure.data$sample_id

rownames(GABAanno_Mola) <- GABAanno_Mola$sample_id
rownames(facs.anno) <- facs.anno$sample_id
######################################################################################################
### Figure: ##########################################################################################
######################################################################################################
GABAmola.facs.anno <- facs.anno[facs.fig1.cells,]
GABAmola.patch.anno <- GABAanno_Mola


n_clusters <- length(unique(GABAmola.facs.anno$cluster_id))
rm.labels <- setdiff(labels(FACS_dend),unique(GABAanno_Mola$cluster_label))

pruned_dend <- prune_dend(rm.labels = rm.labels, dend = FACS_dend, top.level = TRUE )
plot(pruned_dend)

# convert to ggdend
dend_gg <- as.ggdend(pruned_dend)

dend_seg <- dend_gg$segments

cluster_anno <- as.data.frame(GABAmola.facs.anno) %>%
  filter(cluster_label %in% labels(pruned_dend)) %>% 
  select(cluster_id, cluster_label, cluster_color, subclass_label) %>%
  unique()

dend_leaves <- dend_gg$labels %>%
  mutate(cluster_label = label) %>%
  left_join(cluster_anno) %>%
  rowwise() %>%
  mutate(label = sub(paste0(subclass_label," "), "", cluster_label),
         label = sub("VISp |ALM ","",label))

GABAmola.facs.anno <- GABAmola.facs.anno %>%
  left_join(dend_leaves) %>%
  mutate(cluster_id = x)

GABAmola.patch.anno <- GABAmola.patch.anno %>%
  left_join(dend_leaves) %>%
  mutate(cluster_id = x)

panel_pad <- 0.05

# Layer rectangles
layer_rects <- GABAmola.facs.anno %>%
  select(cluster_id, cluster_label, cluster_color, layer_id, layer_label, layer_color) %>%
  group_by(cluster_id, layer_id, layer_label) %>%
  mutate(ly_n = n()) %>%
  ungroup() %>%
  group_by(cluster_id) %>%
  arrange(layer_id) %>%
  mutate(cluster_n = n(),
         ly_frac = ly_n/cluster_n) %>%
  unique() %>%
  arrange(layer_id) %>%
  mutate(ly_cum_frac = cumsum(ly_frac)) %>%
  ungroup() %>%
  arrange(cluster_id, layer_id) %>%
  group_by(cluster_id) %>%
  mutate(xmin = cluster_id - 0.4,
         xmax = cluster_id ,
         ymax = -1 - lag(ly_cum_frac, default = 0)*2 - panel_pad * 2,
         ymin = -1 - ly_cum_frac*2 - panel_pad * 2)

colnames(REF_layer_label) <- c("Revisited_layer_label", "Revisited_layer_id", "Revisited_layer_color")
GABAmola.patch.anno <- GABAmola.patch.anno %>% 
  mutate(Revisited_layer_label = case_when(
    Revisited_layer_label == "VIS1_ho" | Revisited_layer_label == "VISp1" ~ "L1",
    Revisited_layer_label == "VIS23_ho" | Revisited_layer_label == "VISp2/3" ~ "L2/3",
    Revisited_layer_label == "VIS4_ho" | Revisited_layer_label == "VISp4" ~ "L4",
    Revisited_layer_label == "VIS5_ho" | Revisited_layer_label == "VISp5"~ "L5",
    Revisited_layer_label == "VIS6a_ho" | Revisited_layer_label == "VISp6a" ~ "L6",
    Revisited_layer_label == "VIS6b_ho" | Revisited_layer_label == "VISp6b" ~ "L6",
    TRUE ~ Revisited_layer_label))

drop <- c("Revisited_layer_id", "Revisited_layer_color")
GABAmola.patch.anno <- GABAmola.patch.anno[,!colnames(GABAmola.patch.anno) %in% drop]
GABAmola.patch.anno <- left_join(GABAmola.patch.anno, REF_layer_label)
  
patch.layer_rects <- GABAmola.patch.anno %>%
  select(cluster_id, cluster_label, cluster_color, Revisited_layer_id, Revisited_layer_label, Revisited_layer_color) %>%
  group_by(cluster_id, Revisited_layer_id, Revisited_layer_label) %>%
  mutate(ly_n = n()) %>%
  ungroup() %>%
  group_by(cluster_id) %>%
  arrange(Revisited_layer_id) %>%
  mutate(cluster_n = n(),
         ly_frac = ly_n/cluster_n) %>%
  unique() %>%
  arrange(Revisited_layer_id) %>%
  mutate(ly_cum_frac = cumsum(ly_frac)) %>%
  ungroup() %>%
  arrange(cluster_id, Revisited_layer_id) %>%
  group_by(cluster_id) %>%
  mutate(xmin = cluster_id ,
         xmax = cluster_id + 0.4,
         ymax = -1 - lag(ly_cum_frac, default = 0)*2 - panel_pad * 2,
         ymin = -1 - ly_cum_frac*2 - panel_pad * 2)

cluster_bar_thickness <- 0.03

cluster_color_bar <- dend_leaves %>%
  mutate(ymax = min(patch.layer_rects$ymin) - panel_pad,
         ymin = min(patch.layer_rects$ymin) - cluster_bar_thickness - panel_pad,
         xmin= x - 0.4,
         xmax = x + 0.4)

cluster_namepad_thickness <- 1

cluster_namepad_bar <- dend_leaves %>%
  mutate(ymax = min(cluster_color_bar$ymin) ,
         ymin = min(cluster_color_bar$ymin) - cluster_namepad_thickness ,
         xmin= x - 0.4,
         xmax = x + 0.4,
         namepad_color ="#CCE5FF")

# Flat version

flat_plot <- ggplot() +
  geom_segment(data = dend_seg,
               aes(x = x,
                   xend = xend,
                   y = y,
                   yend = yend,
                   size = lwd,
                   color = "black"),
               lineend = "square") +
  geom_rect(data = layer_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin+1,
                ymax = ymax+1,
                fill = layer_color)) +
  geom_rect(data = patch.layer_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin+1,
                ymax = ymax+1,
                fill = Revisited_layer_color)) +
  geom_rect(data = cluster_color_bar,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin+1,
                ymax = ymax+1,
                fill = cluster_color)) +
  geom_rect(data = cluster_namepad_bar,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin+1,
                ymax = ymax+1,
                fill = namepad_color)) +
  geom_text(data = dend_leaves,
            aes(x = x,
                y = -3.15,
                label = cluster_label,
                color = "Black"),
            angle = 90,
            hjust = 0,
            vjust = 0.3,
            size = 3) +
  scale_size(range = c(0.5, 1)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_y_continuous(expand = c(0.25,0)) +
  scale_x_continuous(limits = c(-1,n_clusters + 1)) +
  theme_void()
flat_plot

#ggsave("hierarchical_cluster_flat_plot_bars_no_subclass_t.pdf",
#       flat_plot, 
#       width = 16, height = 12, 
#       useDingbats = F)


# Rotation to match page
flat_plot_r <- ggplot() +
  geom_segment(data = dend_seg,
               aes(x = -y,
                   xend = -yend,
                   y = -x,
                   yend = -xend,
               size = lwd,
               color = col),
               lineend = "square") +
  # Annotation panels
  geom_rect(data = av_rects,
            aes(xmin = -ymin,
                xmax = -ymax,
                ymin = -xmin,
                ymax = -xmax,
                fill = region_color)) +
  geom_rect(data = layer_rects,
            aes(xmin = -ymin,
                xmax = -ymax,
                ymin = -xmin,
                ymax = -xmax,
                fill = layer_color)) +
  geom_text(data = n_rects,
            aes(x = -ymax,
                y = -(xmin + xmax)/2,
                label = n),
            hjust = 0,
            vjust = 0.3,
            size = 2) +
  # Leaf Labels
  geom_text(data = dend_leaves,
            aes(x = 4.4,
                y = -x,
                label = label,
                color = col),
            angle = 0,
            hjust = 0,
            vjust = 0.3,
            size = 2) +
  #Vertical class separators
  geom_segment(data = wedge_lines,
               aes(x = -y,
                   xend = -yend,
                   y = -x,
                   yend = -xend)) +
  #Borders between panels
  geom_segment(data = panel_segments,
               aes(x = -y, xend = -yend,
                   y = -x, yend = -xend)) +
  scale_size(range = c(0.5, 1)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0.25,0)) +
  scale_y_continuous(limits = c(-n_clusters - 2,2)) +
  theme_void()
flat_plot_r

#ggsave("hierarchical_cluster_flat_plot_bars_no_subclass_no_nbars.pdf",
#       flat_plot_r, 
#       width = 7.5, height = 10.5, 
#       useDingbats = F)



dend_nodes <- dend_gg$nodes %>%
  mutate(label = get_nodes_attr(dend,"label")) %>%
  filter(y > 0)

nodes_plot <- ggplot() +
  geom_segment(data = dend_seg,
               aes(x = x,
                   xend = xend,
                   y = y,
                   yend = yend)) +
  geom_text(data = dend_leaves,
            aes(x = x,
                y = -0.01,
                label = paste(label, cluster_id),
                color = col),
            angle = 90,
            hjust = 1,
            vjust = 0.3,
            size = 2) +
  geom_text(data = dend_nodes,
            aes(x = x + 1,
                y = y + 0.02,
                label = label),
            size = 2) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_y_continuous(expand = c(0.25,0)) +
  scale_x_continuous(limits = c(-1,n_clusters + 1)) +
  theme_void()
nodes_plot

#ggsave("hierarchical_cluster_node_labels.pdf",nodes_plot,height = 8, width = 12, useDingbats = F)




cre_anno <- anno %>%
  select(cre_id, cre_label, cre_color) %>%
  unique()

cre_legend <- ggplot() +
  geom_tile(data = cre_anno,
            aes(x = cre_id, y = 1,
                fill = cre_color)) +
  scale_fill_identity() +
  scale_x_continuous(breaks = cre_anno$cre_id, labels = cre_anno$cre_label) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust = 0.3))
#ggsave("cre_legend.pdf",cre_legend,height = 6, width = 12, useDingbats = F)


library(rbokeh)

b_av_rects <- av_rects
b_layer_rects <- layer_rects
names(b_layer_rects) <- names(b_av_rects)
b_cre_rects <- cre_rects
names(b_cre_rects) <- names(b_av_rects)

all_rects <- rbind(b_av_rects, b_layer_rects, b_cre_rects)

bokeh_plot <- figure(height = 800,
                     width = 1200,
                     padding_factor = 0,
                     xaxes = FALSE,
                     yaxes = FALSE,
                     xlab = "",
                     ylab = "") %>%
  ly_segments(data = dend_seg,
              x0 = x, x1 = xend,
              y0 = y, y1 = yend) %>%
  ly_rect(data = all_rects,
          xleft = xmin, xright = xmax,
          ybottom = ymin, ytop = ymax,
          color = region_color,
          fill_alpha = 1,
          hover = list("Cluster" = cluster_label,
                       "Annotation" = region_label,
                       "N Cells" = ly_n,
                       "Fraction" = round(ly_frac,2))) %>%
  ly_rect(data = n_rects,
          xleft = xmin, xright = xmax,
          ybottom = ymin, ytop = ymax,
          color = cluster_color,
          fill_alpha = 1,
          hover = list("Cluster" = cluster_label,
                       "N Cells" = n))
  
bokeh_plot

htmlwidgets::saveWidget(bokeh_plot, file = "all_cell_dendrogram.html")


# 
# flat_plot_av_ratio <- ggplot() +
#   geom_segment(data = dend_seg,
#                aes(x = x,
#                    xend = xend,
#                    y = y,
#                    yend = yend)) +
#   geom_segment(data = sections,
#                aes(x = x,
#                    xend = xend,
#                    y = y,
#                    yend = yend)) +
#   geom_rect(data = av_rects,
#             aes(xmin = xmin,
#                 xmax = xmax,
#                 ymin = ymin - 0.05,
#                 ymax = ymax - 0.05,
#                 fill = region_color)) +
#   geom_segment(data = wedge_lines,
#                aes(x = x,
#                    xend = xend,
#                    y = y,
#                    yend = yend - 0.3)) +
#   geom_text(data = dend_leaves,
#             aes(x = x,
#                 y = -1.07,
#                 label = label,
#                 color = col),
#             angle = 90,
#             hjust = 1,
#             vjust = 0.3,
#             size = 2) +
#   geom_point(data = region_frac,
#              aes(x = x, y = y, color = "#000000", fill = color),
#              size = 2,pch=21) +
#   scale_color_identity() +
#   scale_fill_identity() +
#   scale_y_continuous(expand = c(0.25,0)) +
#   scale_x_continuous(limits = c(-1,n_clusters + 1)) +
#   theme_void()
# flat_plot_av_ratio
# ggsave("hierarchical_cluster_flat_plot_ratio.pdf",flat_plot_av_ratio,height = 8, width = 12, useDingbats = F)

# 
# layer_heat <- anno %>%
#   left_join(layer_order) %>%
#   select(cluster_id, cluster_label, cluster_color, layer_order, layer_label, layer_color) %>%
#   group_by(cluster_id, layer_order, layer_label) %>%
#   mutate(ly_n = n()) %>%
#   ungroup() %>%
#   group_by(cluster_id) %>%
#   arrange(layer_order) %>%
#   mutate(cluster_n = n(),
#          ly_frac = ly_n/cluster_n) %>%
#   unique() %>%
#   arrange(layer_order) %>%
#   ungroup() %>%
#   arrange(cluster_id, layer_order) %>%
#   group_by(cluster_id) %>%
#   mutate(xmin = cluster_id - 0.5,
#          xmax = cluster_id + 0.5,
#          ymax = -1 - (layer_order - 1)/12,
#          ymin = -1 - layer_order/12,
#          ly_frac_color = values_to_colors(ly_frac, max = 1))
# 
# ggplot() +
#   geom_rect(data = av_rects,
#             aes(xmin = xmin, xmax = xmax,
#                 ymin = ymin, ymax = ymax,
#                 fill = region_color)) +
#   geom_rect(data = layer_heat,
#             aes(xmin = xmin, xmax = xmax,
#                 ymin = ymin, ymax = ymax,
#                 fill = ly_frac_color)) +
#   geom_text(data = av_rects,
#             aes(x = cluster_id,
#                 y = -2,
#                 label = cluster_label,
#                 color = cluster_color),
#             angle = 90, hjust = 1, vjust = 0.3) +
#   scale_fill_identity() +
#   scale_color_identity() +
#   scale_y_continuous(limits = c(-3,0)) +
#   theme_void()