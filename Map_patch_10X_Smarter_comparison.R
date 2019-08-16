##########################################################################################
### Setting up the libraries: ############################################################
##########################################################################################

library(dendextend)
library(matrixStats)
library(feather)
library(dplyr)
library(scrattch.io)
library(scrattch.vis)
library(feather)
library(dendextend)
library(tibble)
library(Matrix)
mydir = "/allen//programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/"
source("/allen//programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/heatmap.R")
source("/allen//programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/de.genes.R")
source("/allen//programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/dendro.R")
source("/allen//programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/patchseq.R")
source("/allen//programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Plot_tree_functions.R")
source("/allen//programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")
setwd(mydir)
#devtools::install_github("AllenInstitute/scrattch.io")
options(stringsAsFactors = F)

##########################################################################################
### Input paths: #########################################################################
##########################################################################################

#cl and cl.df for GABA cells only
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/V1/inh.cl.df.rda")
load("/allen//programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/select.markers.rda")

#Tree has 93 types from V1
dend <- readRDS("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/dend.RData")

#norm.dat
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/V1/dat.list.rda")

TenX_cells_norm <- dat.list$`10X_cells`
Smarter_cells_norm <- dat.list$Smartseq_cells

TenX.cells <- colnames(TenX_cells_norm)
Smarter.cells <- unique(colnames(Smarter_cells_norm))

smaller.select.markers <- select.markers[select.markers %in% rownames(TenX_cells_norm)]

rownames(cl.df)=cl.df$cl
cltmp=cl.df[as.character(cl),"cluster_label"]
names(cltmp)=names(cl)
cl=factor(cltmp)
cl <- cl[cl!="Meis2 Adamts19"]

GABA.TenX.cl <- cl[names(cl) %in% TenX.cells]
GABA.Smarter.cl <- cl[unique(names(cl)[names(cl) %in% Smarter.cells])]

GABA.TenX.cells <- names(GABA.TenX.cl)
GABA.Smarter.cells <- names(GABA.Smarter.cl)

rm(TenX.cells)
rm(Smarter.cells)
rm(cl)

GABA_TenX_dat_norm <- TenX_cells_norm[smaller.select.markers,GABA.TenX.cells]
GABA_Smarter_dat_norm <- Smarter_cells_norm[smaller.select.markers, GABA.Smarter.cells]


##########################################################################################
### Patchseq: ############################################################################
##########################################################################################

robject.dir = "/allen//programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
batch_date="20190722_BT014-RSC-215"
patchseq.dir =  "/allen//programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/"

tmp<-load(paste0(robject.dir,batch_date,"_mouse_patchseq_star2.0_cpm.Rdata"))
query.dat = cpmR

# loading samp.dat object
tmp<-load(paste0(robject.dir,batch_date,"_mouse_patchseq_star2.0_samp.dat.Rdata"))

keepcells = which(samp.dat$Region=="VISp" & samp.dat$Type=="patch_seq")
samp.dat = samp.dat[c(keepcells, which(samp.dat$Region=="TCx"),which(samp.dat$Region=="FCx"),which(samp.dat$Region=="MOp"),which(samp.dat$Region=="TEa")   ),]   #FCx is for Brian.  Rat samples mapped in mouse

query.dat = query.dat[,as.character(samp.dat$exp_component_name)]
colnames(query.dat)=as.character(samp.dat$patched_cell_container)

query.dat.norm = log2(as.matrix(query.dat+1))


patchseq_anno <- Read_patchseq_anno(paste0(patchseq.dir,"anno.feather"))
sum(select.markers %in% rownames(query.dat.norm))

dend.list = dend_list(dend)
inh.dend = dend.list[["n59"]]
patch.GABA.cells <- patchseq_anno[patchseq_anno$topLeaf_label %in% labels(inh.dend), "sample_id"]
patch.GABA.cells <- patch.GABA.cells[patch.GABA.cells %in% colnames(query.dat.norm)]
length(patch.GABA.cells)
dim(query.dat.norm)
GABA.query.dat.norm <- query.dat.norm[smaller.select.markers, patch.GABA.cells]

##########################################################################################
### Mapping patchseq on 10X: #############################################################
##########################################################################################

labels(inh.dend)

GABApatch_10X_memb <- map_dend_membership(inh.dend, TenX.cls, GABA_TenX_cells_norm[, GABA.TenX.cells], 
                    GABA.query.dat.norm, patch.GABA.cells, bs.num=100, p=0.7, low.th=0.15, Method = "Mean")

save(GABApatch_10X_memb, file = paste0(mydir, "/GABApatch_10X_memb.rda"))
load(paste0(mydir, "/GABApatch_10X_memb.rda"))

mapped.cl <- colnames(GABApatch_10X_memb)[colnames(GABApatch_10X_memb) %in% labels(inh.dend)]
sorted_GABApatch_10X_memb_cl = t(apply(GABApatch_10X_memb[ , mapped.cl], 1, function(x) names(x)[order(x, decreasing = TRUE)]))


##########################################################################################
### Mapping patchseq on Smarter : ########################################################
##########################################################################################

GABApatch_Smarter_memb <- map_dend_membership(inh.dend, GABA.Smarter.cl, GABA_Smarter_dat_norm[smaller.select.markers,GABA.Smarter.cells], 
                    GABA.query.dat.norm[smaller.select.markers, patch.GABA.cells], 
                    patch.GABA.cells, bs.num=100, p=0.7, low.th=0.15)

save(GABApatch_Smarter_memb, file = paste0(mydir, "/GABApatch_Smarter_memb.rda"))
#load(paste0(mydir, "/GABApatch_Smarter_memb.rda" ))

mapped.cl <- colnames(GABApatch_Smarter_memb)[colnames(GABApatch_Smarter_memb) %in% labels(inh.dend)]
sorted_GABApatch_Smarter_memb_cl = t(apply(GABApatch_Smarter_memb[ , mapped.cl], 1, function(x) names(x)[order(x, decreasing = TRUE)]))

##########################################################################################
### Analysis: ############################################################################
##########################################################################################

sum(sorted_GABApatch_Smarter_memb_cl[patch.GABA.cells,1] == sorted_GABApatch_10X_memb_cl[patch.GABA.cells,1])/length(patch.GABA.cells)


library(reshape2)
patchseq_10X_cl <- as.data.frame(table(sorted_GABApatch_10X_memb_cl[patch.GABA.cells, 1]))
patchseq_Smarter_cl <- as.data.frame(table(sorted_GABApatch_Smarter_memb_cl[patch.GABA.cells, 1]))

tmp <- merge(patchseq_10X_cl, patchseq_Smarter_cl,  by="Var1", all.x = TRUE)
colnames(tmp) <- c("Var1", "patchseq_10X_cl", "patchseq_Smarter_cl")
tmp <- melt(tmp)
colnames(tmp) <- c("Variable", "Method", "Freq")

ggplot(tmp,aes(x=Variable,y=Freq,fill=Method))+
  geom_bar(stat="identity",position="dodge") +
  xlab("Type")+ylab("#number of cells") + theme(axis.text.x = element_text(angle = 90, hjust = 1))



TenX <- as.data.frame(table(GABA.TenX.cl))
colnames(TenX) <- c("cluster_label", "size.Tenx")
TenX[, "size.Tenx"] <- TenX[, "size.Tenx"] /sum(TenX$size.Tenx)
Smarter <- as.data.frame(table(GABA.Smarter.cl))
colnames(Smarter) <- c("cluster_label", "size.Smarter")
Smarter[, "size.Smarter"] <- Smarter[, "size.Smarter"] /sum(Smarter$size.Smarter)
tmp <- merge(TenX, Smarter,  by="cluster_label", all.x = TRUE)
colnames(tmp) <- c("cluster_label", "Tenx", "Smarter")

tmp <- melt(tmp)
colnames(tmp) <- c("Variable", "Method", "Freq")

ggplot(tmp,aes(x=Variable,y=Freq,fill=Method))+
  geom_bar(stat="identity",position="dodge") +
  xlab("Type")+ylab("Porportion of the cells") + theme(axis.text.x = element_text(angle = 90, hjust = 1))



results <- cbind.data.frame(sorted_GABApatch_10X_memb_cl[patch.GABA.cells,c(1,2,3,4,5)], 
                            sorted_GABApatch_Smarter_memb_cl[patch.GABA.cells, c(1,2,3,4,5)])

colnames(results) <- c("10X_first_cl", "10X_second_cl", "10X_third_cl", "10X_forth_cl", "10X_fifth_cl", 
                       "Smarter_first_cl", "Smarter_second_cl", "Smarter_third_cl", "Smarter_forth_cl", "Smarter_fifth_cl")

save(results, file = paste0(mydir, "/confusion.rda"))

  

