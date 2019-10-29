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

tmp.load2 = load(file.path(file=res.dir, file=paste0("V1.dend", bp.name.add,".rda"))) 
tmp.load3 = load(file.path(res.dir, file=paste0("V1.dend.list", bp.name.add,".rda"))) 
patchseq_dend <- dend
rm(dend)

patchseq_anno <- Assign_clusterlabel_clusterid_for_patchseq(patchseq_dir = patchseq_dir, 
                                                            dend = patchseq_dend, 
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
### Layer info: ######################################################################################
######################################################################################################

GABAanno <- patchseq_anno[patchseq_anno$subclass_label %in% c("Vip", "Sst", "Sncg", "Lamp5", "Pvalb"),]
dim(GABAanno)
dim(patchseq_anno)
GABAanno_Mola <- Modify_layer_label_GABAcells(GABAanno)
dim(GABAanno_Mola)  

tmp <- GABAanno %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id ,layer_id, layer_label, layer_color, roi_id, roi_label, roi_color) %>%
  `colnames<-` (c("sample_id", "map_cluster_id", "map_cluster_label", 
                  "map_cluster_color", "cluster_id", "cluster_label", "cluster_color"))

p <- river_plot(tmp, min.cells=0, min.frac=0)
p
save_plot(filename = paste0(work.dir, "/river_roi_layer_label.pdf"), p)

tmp <- GABAanno %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id ,layer_id, layer_label, layer_color, structure_id, structure_label, structure_color) %>%
  `colnames<-` (c("sample_id", "map_cluster_id", "map_cluster_label", 
                  "map_cluster_color", "cluster_id", "cluster_label", "cluster_color"))

p <- river_plot(tmp, min.cells=0, min.frac=0)

tmp<- GABAanno_Mola%>% 
  rownames_to_column("id") %>% 
  dplyr::select(id ,layer_id, layer_label, layer_color, Revisited_layer_id, Revisited_layer_label, Revisited_layer_color) %>%
  `colnames<-` (c("sample_id", "map_cluster_id", "map_cluster_label", 
                  "map_cluster_color", "cluster_id", "cluster_label", "cluster_color"))

p <- river_plot(tmp, min.cells=0, min.frac=0)
p

######################################################################################################
### VISp vs VIS_ho barplots: #########################################################################
######################################################################################################

plot_data <- GABAanno_Mola %>% 
  group_by(Revisited_layer_label, Tree_call_label) %>% 
  summarise(size = n()) %>%
  mutate(Layer = case_when(Revisited_layer_label %in% c("VISp1", "VIS1_ho") ~ "L1",
                           Revisited_layer_label %in% c("VISp2/3", "VIS23_ho") ~ "L2/3",
                           Revisited_layer_label %in% c("VISp4", "VIS4_ho") ~ "L4",
                           Revisited_layer_label %in% c("VISp5", "VIS5_ho") ~ "L5",
                           Revisited_layer_label %in% c("VISp6a", "VIS6a_ho") ~ "L6a",
                           Revisited_layer_label %in% c("VISp6b", "VIS6b_ho") ~ "L6b",
                           TRUE ~ "Missing")) %>%
  mutate(Area = ifelse(Revisited_layer_label %in% c("VISp1", "VISp2/3", "VISp4", "VISp5", "VISp6a", "VISp6b"), "VISp", "VIS_ho"))


p1 <- ggplot(plot_data[plot_data$Tree_call_label=="Core",] ,aes(x=Layer,y=size,fill=Area))+
  geom_bar(stat="identity",position="dodge") +
  xlab("Core")+ylab("Size") + theme(axis.text.x = element_text(angle = 60, hjust = 1))

p2 <-ggplot(plot_data[plot_data$Tree_call_label=="I1",] ,aes(x=Layer,y=size,fill=Area))+
  geom_bar(stat="identity",position="dodge") +
  xlab("I1")+ylab("Size") + theme(axis.text.x = element_text(angle = 60, hjust = 1))

p3 <-ggplot(plot_data[plot_data$Tree_call_label=="I2",] ,aes(x=Layer,y=size,fill=Area))+
  geom_bar(stat="identity",position="dodge") +
  xlab("I2")+ylab("Size") + theme(axis.text.x = element_text(angle = 60, hjust = 1))

p4 <-ggplot(plot_data[plot_data$Tree_call_label=="I3",] ,aes(x=Layer,y=size,fill=Area))+
  geom_bar(stat="identity",position="dodge") +
  xlab("I3")+ylab("Size") + theme(axis.text.x = element_text(angle = 60, hjust = 1))

p5 <-ggplot(plot_data[plot_data$Tree_call_label=="PoorQ",] ,aes(x=Layer,y=size,fill=Area))+
  geom_bar(stat="identity",position="dodge") +
  xlab("PoorQ")+ylab("Size") + theme(axis.text.x = element_text(angle = 60, hjust = 1))


plot_grid(p1,p2,p3,p4,p5)

######################################################################################################
### VISp vs VIS_ho barplots: #########################################################################
######################################################################################################
lamp5_sncg_vip_genes <- c("Lamp5", "Ndnf", "Krt73", "Fam19a1", "Pax6", "Ntn1", "Plch2", "Lsp1", "Lhx6", 
                          "Nkx2-1", "Vip", "Sncg", "Slc17a8", "Nptx2", "Gpr50", "Itih5", "Serpinf1",
                          "Igfbp6", "Gpc3", "Lmo1", "Ptprt", "Rspo4", "Chat", "Crispld2", "Col15a1", "Pde1a")

sst_pvalb_genes <- c("Sst", "Chodl", "Nos1", "Mme", "Tac1", "Tacr3", "Calb2", "Nr2f2", "Myh8", "Tac2",
                     "Hpse", "Crhr2", "Crh", "Esm1", "Rxfp1", "Nts", "Pvalb", "Gabrg1", "Th", "Calb1",
                     "Akr1c18", "Sema3e", "Gpr149", "Reln", "Tpbg", "Cpne5", "Vipr2", "Nkx2-1")

query.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
batch_date="20190722_BT014-RSC-215"

tmp<-load(paste0(query.dir,batch_date,"_mouse_patchseq_star2.0_cpm.Rdata"))
query.dat = cpmR

# loading samp.dat object
tmp<-load(paste0(query.dir,batch_date,"_mouse_patchseq_star2.0_samp.dat.Rdata"))

keepcells = which(samp.dat$Region=="VISp" & samp.dat$Type=="patch_seq")
samp.dat = samp.dat[c(keepcells, which(samp.dat$Region=="TCx"),which(samp.dat$Region=="FCx"),
                      which(samp.dat$Region=="MOp"),which(samp.dat$Region=="TEa")   ),]  #FCx is for Brian.  Rat samples mapped in mouse

query.dat = query.dat[,as.character(samp.dat$exp_component_name)]
colnames(query.dat)=as.character(samp.dat$patched_cell_container)

query.dat.norm = log2(as.matrix(query.dat+1))
query.dat.norm <- query.dat.norm[,Locked_sample_id]

VISp.cells <- rownames(GABAanno_Mola[GABAanno_Mola$Revisited_layer_label %in% c("VISp1", "VISp2/3", "VISp4", "VISp5", "VISp6a", "VISp6b"),])
VIS.ho.cells <- setdiff(rownames(GABAanno_Mola), VISp.cells)

tmp1 <- do.call("cbind", tapply(GABAanno_Mola[VISp.cells, "sample_id"], 
                                GABAanno_Mola[VISp.cells, "cluster_label"], 
                                function(x)rowMeans(query.dat.norm[c(lamp5_sncg_vip_genes, sst_pvalb_genes), x, drop = F])))

tmp2 <- do.call("cbind",tapply(patchseq_anno[VIS.ho.cells, "sample_id"], 
                               patchseq_anno[VIS.ho.cells, "topLeaf_label"],
                               function(x)rowMeans(query.dat.norm[c(lamp5_sncg_vip_genes, sst_pvalb_genes), x, drop = F])))
library(reshape2)
select.cl <- intersect(colnames(tmp1), colnames(tmp2))
tmp1<- tmp1[,select.cl]
tmp2<- tmp2[,select.cl]
tmp1 <- melt(tmp1)
colnames(tmp1) <- c("Gene", "cl", "VISp.med")

tmp2 <- melt(tmp2)
colnames(tmp2) <- c("Gene", "cl", "VISp.ho.med")

tmp3 <- left_join(tmp1, tmp2)

ref.cl.color <- unique(patchseq_anno %>% select(cluster_label, cluster_color))
tmp3$cl_color <- sapply(tmp3$cl, function(x)ref.cl.color[ref.cl.color$cluster_label==x, "cluster_color"])

library(purrr)
library(rbokeh)
.libPaths(c("/allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Olivia/3.5", .libPaths()))

figure(width =800, height = 800,
       xlim = c(0, 20),
       ylim = c(0, 20),
       xlab = "VISp median gene expression",
       ylab = "ho mean gene expression") %>%
  ly_wedge(data = tmp3,
           x = VISp.med,
           y = VISp.ho.med,
           #color = Gene_color,
           start_angle = 0,
           end_angle = 3.99*pi/2,
           radius = 0.1,
           hover = list("Cluster label" = cl,
                        "Gene label" = Gene)) 
