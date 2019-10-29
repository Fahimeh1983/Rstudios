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

#source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Rcodes/color_functions.R")
#source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Rcodes/prune_leaf_custom.R")
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

REF$dendcluster_label == REF$cluster_label
REF$dendcluster_color == REF$cluster_color

REF <- REF %>% arrange(dendcluster_id)
sum(cbind.data.frame(REF[,c("dendcluster_id", "dendcluster_label")], labels(FACS_dend))[,2] == 
  cbind.data.frame(REF[,c("dendcluster_id", "dendcluster_label")], labels(FACS_dend))[,3])


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
unique(GABAanno_Mola[,c("structure_label", "Revisited_layer_label")])
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

######################################################################################################
### Initialization: ##################################################################################
######################################################################################################
anno_file <- paste0(facs_dir, "/anno.feather")
REF.anno.feather <- feather::read_feather(anno_file) 
REF.anno.feather <- REF.anno.feather[REF.anno.feather$sample_id %in% facs.fig1.cells,]


data_file <- paste0(facs_dir, "/data.feather")
REF.data.feather <- feather::feather(data_file)
REF.data.feather <- REF.data.feather[REF.data.feather$sample_id %in% facs.fig1.cells,]
#REF.data.feather <- REF.data.feather[REF.data.feather$sample_id %in% REF.anno.feather$sample_id,]

data_file <- paste0(patchseq_dir, "/data.feather")
map.data.feather <- feather::feather(data_file)
map.data.feather <- map.data.feather[map.data.feather$sample_id %in% patchseq.fig1.cells,]

dendcluster_ids <- 
  sort(unique(patchseq.layer.figure.data[,c("cluster_label", "dendcluster_id")]) [,2])

sst_pvalb_dendcluster_ids <- 
  sort(unique(patchseq.layer.figure.data[patchseq.layer.figure.data$subclass_label %in% c("Sst", "Pvalb"), "dendcluster_id"]))

Lamp5_vip_sncg_dendcluster_ids <- 
  sort(unique(patchseq.layer.figure.data[patchseq.layer.figure.data$subclass_label %in% c("Lamp5", "Sncg", "Vip", "Serpinf1"), "dendcluster_id"]))

fig1a <- build_layer_comparison_FACS_patch_plot(facs.layer.figure.data %>%
                                                  mutate(layer_label = ifelse(layer_label == "L6b", "L6", layer_label)),
                                                patchseq.layer.figure.data %>% 
                                                  mutate(Revisited_layer_label = case_when(
                                                    Revisited_layer_label == "VIS1_ho" ~ "VISp1",
                                                    Revisited_layer_label == "VIS23_ho" ~ "VISp2/3",
                                                    Revisited_layer_label == "VIS4_ho" ~ "VISp4",
                                                    Revisited_layer_label == "VIS5_ho" ~ "VISp5",
                                                    Revisited_layer_label == "VIS6a_ho" ~ "VISp6",
                                                    Revisited_layer_label == "VIS6b_ho" ~ "VISp6",
                                                    TRUE ~ Revisited_layer_label)) %>%
                                                  mutate(Revisited_layer_label = 
                                                           ifelse(Revisited_layer_label %in% c("VISp6b", "VISp6a"), 
                                                                  "VISp6", Revisited_layer_label)),
                                                FACS_dend,
                                                sst_pvalb_dendcluster_ids) 

fig1a <- build_dotplot_comparison_FACS_patch_plot(facs.layer.figure.data,
                                         patchseq.layer.figure.data %>% 
                                           mutate(Revisited_layer_label = case_when(
                                             Revisited_layer_label == "VIS1_ho" ~ "VISp1",
                                             Revisited_layer_label == "VIS23_ho" ~ "VISp2/3",
                                             Revisited_layer_label == "VIS4_ho" ~ "VISp4",
                                             Revisited_layer_label == "VIS5_ho" ~ "VISp5",
                                             Revisited_layer_label == "VIS6a_ho" ~ "VISp6",
                                             Revisited_layer_label == "VIS6b_ho" ~ "VISp6",
                                             TRUE ~ Revisited_layer_label)) %>%
                                           mutate(Revisited_layer_label = 
                                                    ifelse(Revisited_layer_label %in% c("VISp6b", "VISp6a"), 
                                                           "VISp6", Revisited_layer_label)),
                                         FACS_dend,
                                         sst_pvalb_dendcluster_ids, right_pad=10) 

fig1b <- build_layer_comparison_FACS_patch_plot(facs.layer.figure.data,
                                       patchseq.layer.figure.data %>% 
                                        mutate(Revisited_layer_label = case_when(
                                          Revisited_layer_label == "VIS1_ho" ~ "VISp1",
                                          Revisited_layer_label == "VIS23_ho" ~ "VISp2/3",
                                          Revisited_layer_label == "VIS4_ho" ~ "VISp4",
                                          Revisited_layer_label == "VIS5_ho" ~ "VISp5",
                                          Revisited_layer_label == "VIS6a_ho" ~ "VISp6",
                                          Revisited_layer_label == "VIS6b_ho" ~ "VISp6",
                                          TRUE ~ Revisited_layer_label)) %>%
                                         mutate(Revisited_layer_label = 
                                                  ifelse(Revisited_layer_label %in% c("VISp6b", "VISp6a"), 
                                                         "VISp6", Revisited_layer_label)),
                                          FACS_dend,
                                          Lamp5_vip_sncg_dendcluster_ids) 
fig1b <- build_dotplot_comparison_FACS_patch_plot(facs.layer.figure.data,
                                                patchseq.layer.figure.data %>% 
                                                  mutate(Revisited_layer_label = case_when(
                                                    Revisited_layer_label == "VIS1_ho" ~ "VISp1",
                                                    Revisited_layer_label == "VIS23_ho" ~ "VISp2/3",
                                                    Revisited_layer_label == "VIS4_ho" ~ "VISp4",
                                                    Revisited_layer_label == "VIS5_ho" ~ "VISp5",
                                                    Revisited_layer_label == "VIS6a_ho" ~ "VISp6",
                                                    Revisited_layer_label == "VIS6b_ho" ~ "VISp6",
                                                    TRUE ~ Revisited_layer_label)) %>%
                                                  mutate(Revisited_layer_label = 
                                                           ifelse(Revisited_layer_label %in% c("VISp6b", "VISp6a"), 
                                                                  "VISp6", Revisited_layer_label)),
                                                FACS_dend,
                                                Lamp5_vip_sncg_dendcluster_ids) 


sst_pvalb_genes <- c("Sst", "Chodl", "Nos1", "Mme", "Tac1", "Tacr3", "Calb2", "Nr2f2",
                     "Myh8", "Tac2", "Hpse", "Crhr2", "Crh", "Esm1", "Rxfp1", "Nts",
                     "Pvalb", "Gabrg1", "Th", "Calb1", "Akr1c18", "Sema3e", "Gpr149",
                     "Reln", "Tpbg", "Cpne5", "Vipr2", "Nkx2-1")

lamp5_sncg_vip_genes <- c("Lamp5", "Ndnf", "Krt73", "Fam19a1", "Pax6", "Ntn1", "Plch2",
                          "Lsp1", "Lhx6", "Nkx2-1", "Vip", "Sncg", "Slc17a8", "Nptx2",
                          "Gpr50", "Itih5", "Serpinf1", "Igfbp6", "Gpc3", "Lmo1",
                          "Ptprt", "Rspo4", "Chat", "Crispld2", "Col15a1", "Pde1a")

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/Rcodes/layer_and_violin_functions.R")

fig1c<- group_violin_FACS_patch_plot(data_source = facs_dir, 
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
                                       map.anno.feather = GABAanno_Mola[patchseq.fig1.cells,],
                                       map.data.feather = map.data.feather,
                                       dend = FACS_dend)
 
 fig1d <- group_violin_FACS_patch_plot(data_source = fdir, 
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
                              map.anno.feather = GABAanno_Mola[patchseq.fig1.cells,],
                              map.data.feather = map.data.feather,
                              dend = FACS_dend)

 
all_plots <- plot_grid(fig1a,
                       fig1b,
                       fig1c,
                       fig1d,
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


all_plots <- plot_grid(fig1a,
                       fig1b,
                       fig1c,
                       fig1d,
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

#patchseq_sst_pvalb_layers <- build_layer_plot(GABAanno_Mola,
#                                     FACS_dend,
#                                     dendcluster_ids = dendcluster_ids, 
#                                     patch = TRUE,
#                                     modify_layer_label = TRUE)
#p <- plot(patchseq_sst_pvalb_layers)

# patchseq_sst_pvalb_markers <- group_violin_plot2(data_source = pdir, 
#                                                  group_by = "dendcluster", 
#                                                  clusters = unique(GABAanno[,c("cluster_label", "dendcluster_id")]) [,2],
#                                                  genes = c(sst_pvalb_genes, lamp5_sncg_vip_genes),
#                                                  logscale = TRUE,
#                                                  labelheight = 2,
#                                                  max_width = 10,
#                                                  fontsize = 5,
#                                                  showcounts = FALSE, 
#                                                  anno.feather = GABAanno[GABAanno$Revisited_layer_label %in% c("VISp1", "VISp2/3", "VISp4", "VISp5", "VISp6a", "VISp6b"),], 
#                                                  data.feather = data)
