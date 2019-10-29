library(dplyr)
library(feather)
library(ggplot2)
library(tibble)
options(stringsAsFactors = F)


######################################################################################################
### Reading and cleaning annos #######################################################################
######################################################################################################

facs.anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
FACS_dend <- readRDS("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/dend.RData")

#Removing ALM cells
facs.anno <- facs.anno[facs.anno$region_label == "VISp",]
#Removing lowQ cells
facs.anno <- facs.anno[facs.anno$cluster_label %in% labels(FACS_dend),]
dim(facs.anno)

patchseq.anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/anno.feather")
patchseq.anno <- as.data.frame(patchseq.anno)
rownames(patchseq.anno) <- patchseq.anno$sample_id

Locked_sample_id <- read.csv("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/locked_sampleids.csv")
Locked_sample_id <- Locked_sample_id$x

patchseq.anno <- patchseq.anno[patchseq.anno$sample_id %in% Locked_sample_id,]
dim(patchseq.anno)

######################################################################################################
### Selecting GABAergic annos ########################################################################
######################################################################################################

facs.anno <- facs.anno[facs.anno$subclass_label %in% c("Pvalb", "Sst", "Vip", "Sncg", "Lamp5", "Serpinf1"),]
patchseq.anno <- patchseq.anno[patchseq.anno$subclass_label %in% c("Pvalb", "Sst", "Vip", "Sncg", "Lamp5", "Serpinf1"),]
dim(patchseq.anno)
dim(facs.anno)

######################################################################################################
### Selecting common cres ############################################################################
######################################################################################################

common_cres <- intersect(facs.anno$cre_label, patchseq.anno$cre_label)

facs.anno <- facs.anno %>% 
  filter(cre_label %in% common_cres)

patchseq.anno <- patchseq.anno %>%
  filter(cre_label %in% common_cres)

dim(patchseq.anno)
dim(facs.anno)

######################################################################################################
### Selecting core and I1: ###########################################################################
######################################################################################################

patchseq.anno <- patchseq.anno %>%
  filter(Tree_call_label %in% c("Core", "I1"))

dim(patchseq.anno)
dim(facs.anno)

######################################################################################################
### preparing plotting data ##########################################################################
######################################################################################################

cluster_anno <- facs.anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  arrange(cluster_id)%>%
  unique() %>%
  rownames_to_column("id") %>%
  mutate(cluster_id = id) %>%
  column_to_rownames("id") 

#cluster_anno <- patchseq.anno %>%
#  select(topLeaf_id, topLeaf_label, topLeaf_color) %>%
#  unique()

cre_anno <- facs.anno %>%
  select(cre_id, cre_label, cre_color) %>%
  unique() %>%
  arrange(cre_id) %>% 
  rownames_to_column("id") %>%
  mutate(cre_id = id) %>% 
  column_to_rownames("id")

# change the cre_id in facs data
facs.anno <- facs.anno %>%
  select(-cre_id) %>%
  left_join(cre_anno)

facs.anno <- facs.anno %>%
  select(-cluster_id) %>%
  left_join(cluster_anno)

ref_plot_anno <- facs.anno %>%
  filter(inj_type_label == "No Injection") %>%
  filter(facs_label != "RFP-negative")

ref_plot_data <- ref_plot_anno %>%
  group_by(cre_id, cre_label, cre_color,
           cluster_id, cluster_label, cluster_color) %>%
  summarise(n_cells = n())


patchseq.anno$cluster_label <- patchseq.anno$topLeaf_label
patchseq.anno <- patchseq.anno %>% 
  select(-cluster_color) %>%
  select(-cluster_id) %>% 
  left_join(cluster_anno)

patchseq.anno <- patchseq.anno %>%
  select(-cre_id) %>%
  select(-cre_color) %>%
  left_join(cre_anno)

patch_plot_data <- patchseq.anno %>%
  group_by(cre_id, cre_label, cre_color,
           cluster_id, cluster_label, cluster_color) %>%
  summarise(n_cells = n())

max_cre <-  length(unique(facs.anno$cre_id))
max_cluster <- max(as.integer(facs.anno$cluster_id))

bg_rects <- data.frame(xmin = 1:(max_cre/2 +1)*2 - 1.5,
                       xmax = (1:(max_cre/2 + 1) + 1)*2 - 2.5,
                       ymin = 0.5,
                       ymax = max_cluster + 0.5)



class_breaks <- facs.anno %>%
  group_by(subclass_id) %>%
  summarise(yintercept = max(as.integer(cluster_id)) + 0.5)

ggplot() +
  geom_rect(data = bg_rects,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            fill = "#CAD7D7") +
  geom_hline(data = class_breaks,
             aes(yintercept = yintercept),
             size = 0.2, linetype="dashed") +
  geom_point(data = ref_plot_data,
             aes(x = as.integer(cre_id)-0.25, 
                 y = as.integer(cluster_id),
                 size = n_cells,
                 #fill = cluster_color,
                 fill = "#0000CC",
                 color = "#FFFFFF"),
             pch = 21) + 
  geom_point(data = patch_plot_data,
             aes(x = as.integer(cre_id) + 0.25, 
                 y = as.integer(cluster_id),
                 size = n_cells,
                 #fill = cluster_color,
                 fill = "#FF0000",
                 color = "#FFFFFF"),
             pch = 21) +
  geom_text(data = cluster_anno,
            aes(x = -0.1,
                y = as.integer(cluster_id),
                label = cluster_label,
                #color = cluster_color,
                color = "#0000CC"
                ),
            color = "Black",
            hjust = 1,
            vjust = 0.3,
            size = 2.5) +
  scale_color_manual(values = c("#FFFFFF","#FFFFFF"))+
  scale_fill_manual(values = c("#0000CC","#FF0000"))+
  scale_y_reverse("",
                  limits = c(62.5, -0.5),
                  expand = c(0,0)) +
  scale_x_continuous("",
                     expand = c(0,0),
                     limits = c(-5, 28),
                     breaks = as.integer(cre_anno$cre_id),
                     labels = cre_anno$cre_label,
                     position = "top") +
  scale_size_area(max_size = 5,
                  breaks = c(1,10,50,100,200,500)) +
  theme_classic(base_size = 7) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 0, size=8))

ggsave("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Figure/cre_dotplot.pdf",width = 8.25, height = 10.5, useDingbats = F)


# Pan-broad_specific ordering
pan_cre <- c("Gad2-IRES-Cre","Snap25-IRES2-Cre","Slc32a1-IRES-Cre","Slc17a7-IRES2-Cre")
broad_cre <- c("Rbp4-Cre_KL100","Pvalb-IRES-Cre","Vip-IRES-Cre","Sst-IRES-Cre","Ctgf-T2A-dgCre","Htr3a-Cre_NO152","Trib2-F2A-CreERT2","Ntsr1-Cre_GN220","Tlx3-Cre_PL56","Ndnf-IRES2-dgCre","Rorb-IRES2-Cre")

plot_anno2 <- plot_anno %>%
  mutate(cre_class_id = ifelse(cre_label %in% pan_cre, 1, 
                               ifelse(cre_label %in% broad_cre,2,3)))


cre_anno2 <- plot_anno2 %>%
  select(cre_id, cre_label, cre_color, cre_class_id) %>%
  unique() %>%
  arrange(cre_class_id, cre_id) %>%
  mutate(new_cre_id = 2:(n()+1))

plot_anno2 <- plot_anno2 %>%
  left_join(cre_anno2)

plot_data2 <- plot_anno2 %>%
  group_by(new_cre_id, cre_label, cre_color,
           cluster_id, cluster_label, cluster_color) %>%
  summarise(n_cells = n())

ggplot() +
  geom_rect(data = bg_rects,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            fill = "#CAD7D7") +
  geom_point(data = plot_data2,
             aes(x = new_cre_id, 
                 y = cluster_id,
                 size = n_cells,
                 color = cluster_color)) +
  geom_text(data = cluster_anno,
            aes(x = 1,
                y = cluster_id,
                label = cluster_label,
                color = cluster_color),
            hjust = 1,
            vjust = 0.3,
            size = 2*5/6) +
  geom_hline(data = class_breaks,
             aes(yintercept = yintercept),
             size = 0.2) +
  scale_color_identity() +
  scale_y_reverse("",
                  limits = c(116.5, -0.5),
                  expand = c(0,0)) +
  scale_x_continuous("",
                     expand = c(0,0),
                     limits = c(-5, 37),
                     breaks = cre_anno2$new_cre_id,
                     labels = cre_anno2$cre_label,
                     position = "top") +
  scale_size_area(max_size = 3,
                  breaks = c(10,50,100,200,500)) +
  theme_classic(base_size = 7) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 0))

ggsave("cre_grouped_dotplot.pdf",width = 8.25, height = 10.5, useDingbats = F)
