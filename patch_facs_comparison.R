##########################################################################################
### Setting up the libraries: ############################################################
##########################################################################################
.libPaths("/home/fahimehb/R/x86_64-redhat-linux-gnu-library/3.5")
.libPaths(c("/allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Olivia/3.5", .libPaths()))
library(rbokeh)
library(ggplot2)
library(dendextend)
script_rep="//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/"
library(reshape2)
library(matrixStats)
library(feather)
library(tibble)
library(dplyr)
library(purrr)
library(cowplot)
library(scrattch.hicat)
source(file.path(script_rep,"patchseq/heatmap.R"))
source(file.path(script_rep,"patchseq/de.genes.R"))
source(file.path(script_rep,"patchseq/dendro.R"))
source(file.path(script_rep,"patchseq/patchseq.R"))
source(file.path(script_rep,"patchseq/Mapping_helper_functions.R"))
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")
options(stringsAsFactors = F)


### TODO: The batch date should be updated everytime
#just update batch_date and source it
batch_date="20190722_BT014-RSC-215"

######################################################################################################
### Setting up some paths ############################################################################
######################################################################################################
### TODO: Change these if you changed them in the mapping.R
ref.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"
res.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"
query.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
work.dir = "/allen//programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/"
patchseq.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/"
latest.mapping.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20190729_collapsed40_cpm/"

######################################################################################################
### Reading ref data #################################################################################
######################################################################################################
### TODO: the following two lines should be modified if the number of FACS cells has changed  
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/REF_mapping_probability.rda")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/select.markers.rda")

bp.collapse.th = 40
bp.name.add = NULL
if (!is.null(bp.collapse.th)) {
  bp.name.add = paste0(".with.bp.", bp.collapse.th)
}

###load reference data and tree
tmp.load1 = load(file=file.path(res.dir, "ref.data.rda")) # should include cl, cl.df, norm.dat. # The loaded cl is not used because it only includes cluster ids but not cluster labels 
tmp.load2 = load(file.path(file=res.dir, file=paste0("V1.dend", bp.name.add,".rda"))) # should include the pruned V1 tree
tmp.load3 = load(file.path(res.dir, file=paste0("V1.dend.list", bp.name.add,".rda"))) # should include dend.list

plot(dend)

rownames(cl.df)=cl.df$cluster_id
cltmp=cl.df[as.character(cl),"cluster_label"]
names(cltmp)=names(cl)
cl=factor(cltmp)

FACS.cells <- colnames(norm.dat)
FACS.anno <- read_feather(paste0(ref.dir, "/anno.feather"))
FACS.anno <- as.data.frame(FACS.anno)
rownames(FACS.anno) <- FACS.anno$sample_id
FACS.anno <- FACS.anno[FACS.cells,]


select.genotypes <- unique(FACS.anno$genotype_label)[grepl("Gad2", unique(FACS.anno$genotype_label)) |
                                                       grepl("32", unique(FACS.anno$genotype_label))]

select.subclass <- c("Sst", "Pvalb", "Vip", "Sncg", "Lamp5")

FACS.anno <- FACS.anno %>% 
  filter(subclass_label %in% select.subclass) #%>%
  #filter(genotype_label %in% select.genotypes)

dim(FACS.anno)

######################################################################################################
### Loading the query data ###########################################################################
######################################################################################################

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
idx=match(rownames(norm.dat), rownames(query.dat.norm))
query.dat.norm=query.dat.norm[idx,]

patchseq_anno <- read_feather(paste0(patchseq.dir, "/anno.feather"))
patchseq_anno <- as.data.frame(patchseq_anno)
rownames(patchseq_anno) <- patchseq_anno$sample_id
#Patchseq Cells of interests
locked_cells_spec_id = rownames(read.csv(paste0(work.dir, "mouse_met_Jun_14.csv"), check.names=FALSE, row.names = 1 ))
locked_cells_sample_id = patchseq_anno[patchseq_anno$spec_id_label %in% locked_cells_spec_id, "sample_id"]
length(locked_cells_spec_id) == length(locked_cells_sample_id)
patchseq_anno <- patchseq_anno[locked_cells_sample_id,]
dim(query.dat.norm[select.markers, locked_cells_sample_id])

Core.cells <- patchseq_anno[patchseq_anno$Tree_call_label == "Core", "sample_id"]
I1.cells <- patchseq_anno[patchseq_anno$Tree_call_label == "I1", "sample_id"]
I2.cells <- patchseq_anno[patchseq_anno$Tree_call_label == "I2", "sample_id"]
I3.cells <- patchseq_anno[patchseq_anno$Tree_call_label == "I3", "sample_id"]
patchseq_anno <- patchseq_anno[patchseq_anno$subclass_label %in% select.subclass,]

######################################################################################################
### Comparing gene expression for FACS and patchseq: #################################################
######################################################################################################
lamp5_sncg_vip_genes <- c("Lamp5", "Ndnf", "Krt73", "Fam19a1", "Pax6", "Ntn1", "Plch2", "Lsp1", "Lhx6", 
                          "Nkx2-1", "Vip", "Sncg", "Slc17a8", "Nptx2", "Gpr50", "Itih5", "Serpinf1",
                          "Igfbp6", "Gpc3", "Lmo1", "Ptprt", "Rspo4", "Chat", "Crispld2", "Col15a1", "Pde1a")

sst_pvalb_genes <- c("Sst", "Chodl", "Nos1", "Mme", "Tac1", "Tacr3", "Calb2", "Nr2f2", "Myh8", "Tac2",
                     "Hpse", "Crhr2", "Crh", "Esm1", "Rxfp1", "Nts", "Pvalb", "Gabrg1", "Th", "Calb1",
                     "Akr1c18", "Sema3e", "Gpr149", "Reln", "Tpbg", "Cpne5", "Vipr2", "Nkx2-1")


tmp1 <- do.call("cbind",tapply(FACS.anno[FACS.anno$subclass_label %in% c("Sst", "Pvalb"), "sample_id"], 
                               FACS.anno[FACS.anno$subclass_label %in% c("Sst", "Pvalb"), "cluster_label"], 
                               function(x)rowMedians(norm.dat[sst_pvalb_genes, x, drop = F])))

tmp2 <- do.call("cbind",tapply(patchseq_anno[patchseq_anno$subclass_label %in% c("Sst", "Pvalb"), "sample_id"], 
                               patchseq_anno[patchseq_anno$subclass_label %in% c("Sst", "Pvalb"), "topLeaf_label"],
                               function(x)rowMedians(query.dat.norm[sst_pvalb_genes, x, drop = F])))

cormat <- cor(tmp1[,colnames(tmp2)], tmp2)
tmp3 <- unique((FACS.anno[FACS.anno$subclass_label %in% c("Sst", "Pvalb"),c("cluster_id", "cluster_label", "cluster_color")]))
fake.cl <- setNames(tmp3$cluster_id, tmp3$cluster_label)
col <- t(as.data.frame(setNames(tmp3$cluster_color, tmp3$cluster_label)))
plot_co_matrix(co.ratio =cormat , cl = fake.cl, max.cl.size = 1, col = col)


tmp1 <- do.call("cbind",tapply(FACS.anno[FACS.anno$subclass_label %in% c("Vip", "Lamp5"), "sample_id"], 
                             FACS.anno[FACS.anno$subclass_label %in% c("Vip",  "Lamp5"), "cluster_label"], 
                               function(x)rowMedians(norm.dat[lamp5_sncg_vip_genes, x, drop = F])))

tmp2 <- do.call("cbind",tapply(patchseq_anno[patchseq_anno$subclass_label %in% c("Vip", "Lamp5"), "sample_id"], 
                               patchseq_anno[patchseq_anno$subclass_label %in% c("Vip",  "Lamp5"), "topLeaf_label"],
                               function(x)rowMedians(query.dat.norm[lamp5_sncg_vip_genes, x, drop = F])))

tmp1 <- tmp1[,colnames(tmp2)]
dim(tmp1)
dim(tmp2)
cormat <- cor(tmp1[,colnames(tmp2)], tmp2)
tmp3 <- unique((FACS.anno[FACS.anno$subclass_label %in% c("Vip",  "Lamp5"),c("cluster_id", "cluster_label", "cluster_color")]))
fake.cl <- setNames(tmp3$cluster_id, tmp3$cluster_label)
fake.cl <- fake.cl[rownames(cormat)]
col <- t(as.data.frame(setNames(tmp3$cluster_color, tmp3$cluster_label)))
plot_co_matrix(co.ratio =cormat , cl = fake.cl, max.cl.size = 1)

dim(cormat)
ggplot(data = melt(cor(tmp1[,colnames(tmp2)], tmp2)), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Clustering lable") + ylab("NN Mapping lables") +  
  scale_fill_gradient(low = "white", high = "red")
plot(norm.dat["Lamp5",])


######################################################################################################
### Comparing gene expression for FACS and patchseq: #################################################
######################################################################################################
tmp1 <- do.call("cbind", tapply(FACS.anno[, "sample_id"], 
                               FACS.anno[, "cluster_label"], 
                               function(x)rowMeans(norm.dat[c(lamp5_sncg_vip_genes, sst_pvalb_genes), x, drop = F])))

tmp2 <- do.call("cbind",tapply(patchseq_anno[c(Core.cells, I1.cells), "sample_id"], 
                               patchseq_anno[c(Core.cells, I1.cells), "topLeaf_label"],
                               function(x)rowMeans(query.dat.norm[c(lamp5_sncg_vip_genes, sst_pvalb_genes), x, drop = F])))

select.cl <- intersect(colnames(tmp1), colnames(tmp2))
tmp1<- tmp1[,select.cl]
tmp2<- tmp2[,select.cl]
rownames(tmp1) <- c(lamp5_sncg_vip_genes, sst_pvalb_genes)
rownames(tmp2) <- c(lamp5_sncg_vip_genes, sst_pvalb_genes)
tmp1 <- melt(tmp1)
colnames(tmp1) <- c("Gene", "cl", "FACS.med")
#tmp1$color <- c("Red")

tmp2 <- melt(tmp2)
colnames(tmp2) <- c("Gene", "cl", "Patchseq.med")
#tmp2$color <- c("Green")

tmp3 <- left_join(tmp1, tmp2)

ref.cl.color <- unique(FACS.anno %>% select(cluster_label, cluster_color))
tmp3$cl_color <- sapply(tmp3$cl, function(x)ref.cl.color[ref.cl.color$cluster_label==x, "cluster_color"])
tmp3 <- tmp3 %>% mutate(Gene_color = case_when(Gene == "Sst" ~ "#F54607",
                                               Gene == "Pvalb" ~"#0707F5",
                                               Gene == "Lamp5" ~ "#48F214",
                                               TRUE ~ "#CBCAC6"))

figure(width =800, height = 800,
       xlim = c(0, 20),
       ylim = c(0, 20),
       xlab = "Patchseq mean gene expression",
       ylab = "FACS mean gene expression") %>%
  ly_wedge(data = tmp3,
           x = Patchseq.med,
           y = FACS.med,
           color = Gene_color,
           start_angle = 0,
           end_angle = 3.99*pi/2,
           radius = 0.1,
           hover = list("Cluster label" = cl,
                        "Gene label" = Gene)) 
  

p <- ggplot()+
     geom_point(data = tmp3,
           aes(y = FACS.med,
               x = Patchseq.med,
               colour = Gene_color),
           size = 2,
           show.legend = FALSE) +
     scale_color_identity() + 
     geom_abline(slope = 1) + 
     xlim(0,15) + 
     ylim(0,15) +
     xlab("Patchseq mean gene expression") +
     ylab("FACS mean gene expression")

















