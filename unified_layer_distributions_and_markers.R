.libPaths("/home/fahimehb/R/x86_64-redhat-linux-gnu-library/3.5")
#.libPaths(c("/allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Olivia/3.5", .libPaths()))
library(dplyr)
library(ggplot2)
library(cowplot)
library(feather)
library(dendextend)
library(scrattch.vis)
library(scrattch.io)
library(tibble)
library("cowplot")

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Rcodes/color_functions.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Rcodes/prune_leaf_custom.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Rcodes/layer_and_violin_functions.R")
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Plot_utils.R")

options(stringsAsFactors = F)

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
REF <- REF %>% filter(!is.na(dendcluster_id)) %>% arrange(dendcluster_id)
sum(cbind.data.frame(REF[,c("dendcluster_id", "dendcluster_label")], labels(FACS_dend))[,2] == 
  cbind.data.frame(REF[,c("dendcluster_id", "dendcluster_label")], labels(FACS_dend))[,3])


######################################################################################################
### Locked down data-set #############################################################################
######################################################################################################
res.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"
patchseq_dir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/"

bp.collapse.th = 40
bp.name.add = NULL
if (!is.null(bp.collapse.th)) {
  bp.name.add = paste0(".with.bp.", bp.collapse.th)
}

patchseq_anno <- Assign_clusterlabel_clusterid_for_patchseq(patchseq_dir = patchseq_dir, 
                                                            dend = FACS_dend, 
                                                            ref_cluster_label_color = REF)

patchseq_anno <- as.data.frame(patchseq_anno)
rownames(patchseq_anno) <- patchseq_anno$sample_id

temp1 <- read_feather(paste0(patchseq_dir, "/anno.feather"))
temp1 <- as.data.frame(temp1)
rownames(temp1) <- temp1$sample_id

sum(patchseq_anno$topLeaf_label == temp1$topLeaf_label) == dim(patchseq_anno)[1]
unique(patchseq_anno$dendcluster_label) %in% labels(FACS_dend)

temp2 <- cbind.data.frame(as.data.frame(unique(patchseq_anno[,c("dendcluster_id", "dendcluster_label")])) %>% arrange(dendcluster_id), 
                 labels(FACS_dend)[labels(FACS_dend) %in% unique(patchseq_anno$cluster_label)])
temp2[,2] == temp2[,3]

Locked_sample_id <- read.csv(paste0(work.dir, "/locked_sampleids.csv"))
Locked_sample_id <- Locked_sample_id$x
patchseq_anno <- as.data.frame(patchseq_anno)

patchseq_anno <- patchseq_anno[patchseq_anno$sample_id %in% Locked_sample_id,]
core_cells <- patchseq_anno[patchseq_anno$Tree_call_label == "Core","sample_id"]
I1_cells <- patchseq_anno[patchseq_anno$Tree_call_label == "I1", "sample_id"]
I2_cells <- patchseq_anno[patchseq_anno$Tree_call_label == "I2", "sample_id"]
I3_cells <- patchseq_anno[patchseq_anno$Tree_call_label == "I3", "sample_id"]

dim(patchseq_anno)

######################################################################################################
### Subsetting GABA: #################################################################################
######################################################################################################

GABAanno <- patchseq_anno[patchseq_anno$subclass_label %in% c("Vip", "Sst", "Sncg", "Lamp5", "Pvalb"),]
dim(GABAanno)
GABAanno %>% group_by(Tree_call_label) %>% summarise(n())

######################################################################################################
### correcting layer info: ###########################################################################
######################################################################################################
#Here it will only take cells that are in this subclasses and in specific layers
dim(patchseq_anno)
dim(GABAanno)
GABAanno_Mola <- Modify_layer_label_GABAcells(GABAanno)
dim(GABAanno_Mola)  
study.cells <- setdiff(rownames(GABAanno), rownames(GABAanno_Mola))
unique(patchseq_anno[study.cells, "subclass_label"])
unique(patchseq_anno[study.cells, "structure_label"])

######################################################################################################
### Initialization: ##################################################################################
######################################################################################################

#cluster_anno <- anno %>%
#  select(dendcluster_id, cluster_id, cluster_label, cluster_color) %>%
#  unique()


#patchseq_cluster_anno <- patchseq_anno %>%
#  select(dendcluster_id, cluster_id, cluster_label, cluster_color) %>%
#  unique()

fdir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"
pdir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current"
#work.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/"


data_file <- paste0(fdir, "/data.feather")
REF.data.feather <- feather::feather(data_file)

anno_file <- paste0(fdir, "/anno.feather")
REF.anno.feather <- feather::read_feather(anno_file) 

data_file <- paste0(pdir, "/data.feather")
map.data.feather <- feather::feather(data_file)


#sst_pvalb_layers <- build_layer_plot(anno,
#                              FACS_dend,
#                              dendcluster_ids = 85:115)
#p <- plot(sst_pvalb_layers)
#save_plot(filename = paste0(work.dir, "/sst_pvalb_layers.pdf"), p)

dendcluster_ids <- sort(unique(GABAanno_Mola[,c("cluster_label", "dendcluster_id")]) [,2])
sst_pvalb_dendcluster_ids <- sort(unique(GABAanno_Mola[GABAanno_Mola$subclass_label %in% c("Sst", "Pvalb"), "dendcluster_id"]))
Lamp5_vip_sncg_dendcluster_ids <- sort(unique(GABAanno_Mola[GABAanno_Mola$subclass_label %in% c("Lamp5", "Sncg", "Vip"), "dendcluster_id"]))

patchseq_sst_pvalb_layers <- build_layer_plot(GABAanno_Mola,
                                     FACS_dend,
                                     dendcluster_ids = dendcluster_ids, 
                                     patch = TRUE,
                                     modify_layer_label = TRUE)
p <- plot(patchseq_sst_pvalb_layers)
p

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Rcodes/layer_and_violin_functions.R")

build_layer_comparison_FACS_patch_plot(facs.anno,
                                      GABAanno_Mola %>% 
                                        mutate(Revisited_layer_label = case_when(Revisited_layer_label == "VIS1_ho" ~ "VISp1",
                                                                                                 Revisited_layer_label == "VIS23_ho" ~ "VISp2/3",
                                                                                                 Revisited_layer_label == "VIS4_ho" ~ "VISp4",
                                                                                                 Revisited_layer_label == "VIS5_ho" ~ "VISp5",
                                                                                                 Revisited_layer_label == "VIS6a_ho" ~ "VISp6a",
                                                                                                 Revisited_layer_label == "VIS6b_ho" ~ "VISp6b",
                                                                                                 TRUE ~ Revisited_layer_label)),
                                        dend,
                                      Lamp5_vip_sncg_dendcluster_ids) 
save_plot(filename = paste0(work.dir, "/patchseq_sst_pvalb_layers.pdf"), p)

lamp5_sncg_vip_layers <- build_layer_plot(anno,
                                     FACS_dend,
                                     dendcluster_ids = 56:84)
p <- plot(lamp5_sncg_vip_layers)
save_plot(filename = paste0(work.dir, "/lamp5_sncg_vip_layers.pdf"), p)
p
patchseq_lamp5_sncg_vip_layers <- build_layer_plot(patchseq_anno,
                                          patchseq_dend,
                                          #cocl,
                                          dendcluster_ids = 33:61,
                                          patch = TRUE,
                                          modify_layer_label = TRUE)
p <- plot(patchseq_lamp5_sncg_vip_layers)
save_plot(filename = paste0(work.dir, "/patchseq_lamp5_sncg_vip_layers.pdf"), p)



sst_pvalb_genes <- c("Sst","Chodl","Nos1","Mme","Tac1","Tacr3","Calb2","Nr2f2",
                     "Myh8","Tac2","Hpse","Crhr2","Crh","Esm1","Rxfp1","Nts",
                     "Pvalb","Gabrg1","Th","Calb1","Akr1c18","Sema3e","Gpr149",
                     "Reln","Tpbg","Cpne5","Vipr2","Nkx2-1")

lamp5_sncg_vip_genes <- c("Lamp5","Ndnf","Krt73","Fam19a1","Pax6","Ntn1","Plch2",
                          "Lsp1","Lhx6","Nkx2-1","Vip","Sncg","Slc17a8","Nptx2",
                          "Gpr50","Itih5","Serpinf1","Igfbp6","Gpc3","Lmo1",
                          "Ptprt","Rspo4","Chat","Crispld2","Col15a1","Pde1a")

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Rcodes/layer_and_violin_functions.R")

 group_violin_FACS_patch_plot(data_source = fdir, 
                                       group_by = "dendcluster", 
                                       clusters = sst_pvalb_dendcluster_ids,
                                       genes = sst_pvalb_genes,
                                       logscale = TRUE,
                                       labelheight = 2,
                                       max_width = 10,
                                       fontsize = 10,
                                       showcounts = FALSE, 
                                       REF.anno.feather = REF.anno.feather,
                                       REF.data.feather = REF.data.feather,
                                       map.anno.feather = GABAanno_Mola,
                                       map.data.feather = map.data.feather,
                                       dend = FACS_dend)
 
 group_violin_FACS_patch_plot(data_source = fdir, 
                              group_by = "dendcluster", 
                              clusters = Lamp5_vip_sncg_dendcluster_ids,
                              genes = lamp5_sncg_vip_genes,
                              logscale = TRUE,
                              labelheight = 2,
                              max_width = 10,
                              fontsize = 10,
                              showcounts = FALSE, 
                              REF.anno.feather = REF.anno.feather,
                              REF.data.feather = REF.data.feather,
                              map.anno.feather = GABAanno_Mola,
                              map.data.feather = map.data.feather,
                              dend = FACS_dend)

 
save_plot(filename = paste0(work.dir, "/sst_pvalb_markers.pdf"), p)

source("layer_and_violin_functions.R")
patchseq_sst_pvalb_markers <- group_violin_plot2(data_source = pdir, 
                                        group_by = "dendcluster", 
                                        clusters = unique(GABAanno[,c("cluster_label", "dendcluster_id")]) [,2],
                                        genes = c(sst_pvalb_genes, lamp5_sncg_vip_genes),
                                        logscale = TRUE,
                                        labelheight = 2,
                                        max_width = 10,
                                        fontsize = 5,
                                        showcounts = FALSE, 
                                        anno.feather = GABAanno[GABAanno$Revisited_layer_label %in% c("VISp1", "VISp2/3", "VISp4", "VISp5", "VISp6a", "VISp6b"),], 
                                        data.feather = data)

p <- plot(patchseq_sst_pvalb_markers)
save_plot(filename = paste0(work.dir, "/patchseq_sst_pvalb_markers.pdf"), p)## Violin plots
lamp5_sncg_vip_markers <- group_violin_plot2(data_source = fdir, 
                                        group_by = "dendcluster", 
                                        clusters = 56:84,
                                        genes = lamp5_sncg_vip_genes,
                                        logscale = TRUE,
                                        labelheight = 2,
                                        max_width = 10,
                                        fontsize = 5,
                                        showcounts = FALSE)
p <- plot(lamp5_sncg_vip_markers)
save_plot(filename = paste0(work.dir, "/lamp5_sncg_vip_markers.pdf"), p)
patchseq_lamp5_sncg_vip_markers <- group_violin_plot2(data_source = pdir, 
                                        group_by = "dendcluster", 
                                        clusters = unique(patchseq_anno[patchseq_anno$subclass_label %in% c("Vip", "Sncg", "Lamp5"), "dendcluster_id"]),
                                        genes = lamp5_sncg_vip_genes,
                                        logscale = TRUE,
                                        labelheight = 2,
                                        max_width = 10,
                                        fontsize = 5,
                                        showcounts = FALSE, 
                                        anno.feather = patchseq_anno, 
                                        data.feather = patchseq_data_feather)
p <- plot(patchseq_lamp5_sncg_vip_markers)
save_plot(filename = paste0(work.dir, "/patchseq_lamp5_sncg_vip_markers.pdf"), p)
all_plots <- plot_grid(sst_pvalb_layers,
                       patchseq_sst_pvalb_layers,
                       sst_pvalb_markers,
                       patchseq_sst_pvalb_markers,
                       align = "v",
                       nrow = 2,
                       rel_widths = 1,
                       rel_heights = 1,
                       labels = c("c","d","e","f"))

plot(all_plots)
save_plot(paste0(work.dir, "Gird_FACS_Patchseq_Sst_Pvalb.pdf"),
          all_plots,
          ncol = 2,
          nrow = 2,
          base_width = 7.5/2,
          base_height = 6/2)


all_plots <- plot_grid(lamp5_sncg_vip_layers,
                       patchseq_lamp5_sncg_vip_layers,
                       lamp5_sncg_vip_markers,
                       patchseq_lamp5_sncg_vip_markers,
                       align = "v",
                       nrow = 2,
                       rel_widths = 1,
                       rel_heights = 1,
                       labels = c("c","d","e","f"))

plot(all_plots)
save_plot(paste0(work.dir, "Grid_FACS_Patchseq_lamp5_sncg_vip.pdf"),
          all_plots,
          ncol = 2,
          nrow = 2,
          base_width = 7.5/2,
          base_height = 6/2)
