library(dplyr)
library(feather)
library(ggplot2)
library("philentropy")
library(reshape2)

options(stringsAsFactors = F)
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")
work.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/"
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
### Selecting Core and I1 cells for patchseq annos ###################################################
######################################################################################################

patchseq.anno <- patchseq.anno[patchseq.anno$Tree_call_label %in% c("Core", "I1"),]
dim(patchseq.anno)
dim(facs.anno)
facs.cells <- facs.anno$sample_id
patch.cells <- patchseq.anno$sample_id

######################################################################################################
### Gene lists #######################################################################################
######################################################################################################

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/select.markers.rda")
gene_list_4k <- select.markers
rm(select.markers)
gene_list <- read.csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/good_genes_beta_score.csv")
gene_list_12k <- gene_list$Gene
gene_list_1k <- gene_list[gene_list$BetaScore >= 0.425,"Gene"]
gene_list_45k <- rownames(facs.prob.mat)
######################################################################################################
### reading and Converting norm.dat to prob ##########################################################
######################################################################################################

ref.data.rda.path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"
load(paste0(ref.data.rda.path, "ref.data.rda"))
dim(norm.dat)
facs.prob.mat <- (2 ^ norm.dat) - 1
dim(facs.prob.mat)
facs.prob.mat <- facs.prob.mat[,facs.cells]
facs.prob.mat[facs.prob.mat==0] <- 0.00000001
facs.prob.mat <- t(t(facs.prob.mat)/colSums(facs.prob.mat))

query.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
batch_date="20190722_BT014-RSC-215"
tmp<-load(paste0(query.dir,batch_date,"_mouse_patchseq_star2.0_cpm.Rdata"))
query.dat = cpmR
# loading samp.dat object
tmp<-load(paste0(query.dir,batch_date,"_mouse_patchseq_star2.0_samp.dat.Rdata"))
keepcells = which(samp.dat$Region=="VISp" & samp.dat$Type=="patch_seq")
samp.dat = samp.dat[c(keepcells, which(samp.dat$Region=="TCx"),which(samp.dat$Region=="FCx"),which(samp.dat$Region=="MOp"),which(samp.dat$Region=="TEa")   ),]   #FCx is for Brian.  Rat samples mapped in mouse
query.dat = query.dat[,as.character(samp.dat$exp_component_name)]
colnames(query.dat)=as.character(samp.dat$patched_cell_container)
colnames(query.dat)
patch.prob.mat <- query.dat[,patch.cells]
dim(patch.prob.mat)
rm(samp.dat,  tmp)
patch.prob.mat[patch.prob.mat==0] <- 0.00000001
patch.prob.mat <- t(t(patch.prob.mat)/colSums(patch.prob.mat))

######################################################################################################
### Computing JSD using all the genes ################################################################
######################################################################################################
select.markers <- gene_list_1k
facs.prob.cl.med <- do.call("cbind",tapply(names(cl[facs.cells]), 
                                     cl[facs.cells],
                                     function(x){rowMeans(as.matrix(facs.prob.mat[select.markers,x,drop=F]))}))

dim(facs.prob.cl.med)
facs.prob.cl.med <- t(t(facs.prob.cl.med)/colSums(facs.prob.cl.med))
sum(colSums(facs.prob.cl.med) ==1)
sum(colSums(facs.prob.cl.med) >0.9999999 & colSums(facs.prob.cl.med) <1.0000001)
dim(facs.prob.cl.med)

facs.prob.cl.med <- t(facs.prob.cl.med)
dim(facs.prob.cl.med)
facs.jsd <- JSD(facs.prob.cl.med, test.na = TRUE, unit = "log2", est.prob = NULL)

dim(facs.jsd)

row.and.colnames <- Renew_list(ls = as.numeric(rownames(facs.prob.cl.med)), 
                               ref.df = cl.df, 
                               label = "cluster_id", 
                               new.label = "cluster_label")

colnames(facs.jsd) <- row.and.colnames
rownames(facs.jsd) <- row.and.colnames


ggplot(data = melt(facs.jsd), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()+ theme(axis.text = element_text(size=7)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("") +
  scale_fill_gradient(low = "blue", high = "red" )

#Lets sort each row:
sorted_facs.jsd <- t(apply(facs.jsd, 1, function(x) names(x)[order(x, decreasing = FALSE)]))
sorted_facs.jsd.value <- t(apply(facs.jsd, 1, function(x) x[order(x, decreasing = FALSE)]))

i = 49; head(sorted_facs.jsd[i,]) ; head(sorted_facs.jsd.value[i,])

facs.jsd.1kgenes <- facs.jsd 
write.csv(facs.jsd.1kgenes, file = paste0(work.dir, "facs_jsd_1kgenes.csv"))

######################################################################################################
### Computing cor using all the genes ################################################################
######################################################################################################
norm.dat <- norm.dat[,facs.GABA.cells]

facs.cl.med <- do.call("cbind",tapply(names(cl[facs.GABA.cells]), 
                                           cl[facs.GABA.cells],
                                           function(x){rowMeans(as.matrix(norm.dat[select.markers,x,drop=F]))}))
dim(facs.cl.med)
facs.cor.cl <- cor(facs.cl.med, facs.cl.med)

row.and.colnames <- Renew_list(ls = as.numeric(colnames(facs.cl.med)), 
                               ref.df = cl.df, 
                               label = "cluster_id", 
                               new.label = "cluster_label")
colnames(facs.cor.cl) <- row.and.colnames
rownames(facs.cor.cl) <- row.and.colnames


ggplot(data = melt(facs.cor.cl), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()+ theme(axis.text = element_text(size=7)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("") +
  scale_fill_gradient(low = "white", high = "red")

#Lets sort each row:
sorted_facs.cor.cl <- t(apply(facs.cor.cl, 1, function(x) names(x)[order(x, decreasing = TRUE)]))
sorted_facs.cor.cl.value <- t(apply(facs.cor.cl, 1, function(x) x[order(x, decreasing = TRUE)]))

i = 49; head(sorted_facs.cor.cl[i,]) ; head(sorted_facs.cor.cl.value[i,])

######################################################################################################
### Computing JSD patchseq data ######################################################################
######################################################################################################
query.dat <- query.dat[,patch.cells]
select.markers <- gene_list_45k

patch.cl <- setNames(patchseq.anno$topLeaf_label, patchseq.anno$sample_id)

patch.prob.cl.med <- do.call("cbind",tapply(names(patch.cl[patch.cells]), 
                                       patch.cl[patch.cells],
                                      function(x){rowMeans(as.matrix(patch.prob.mat[select.markers,x,drop=F]))}))

dim(patch.prob.cl.med)
patch.prob.cl.med <- t(t(patch.prob.cl.med)/colSums(patch.prob.cl.med))
sum(colSums(patch.prob.cl.med) ==1)
sum(colSums(patch.prob.cl.med) >0.9999999 & colSums(patch.prob.cl.med) <1.0000001)
dim(patch.prob.cl.med)

patch.prob.cl.med <- t(patch.prob.cl.med)
dim(patch.prob.cl.med)
patch.jsd <- JSD(patch.prob.cl.med, test.na = TRUE, unit = "log2", est.prob = NULL)

dim(patch.jsd)
rownames(patch.jsd) <- rownames(patch.prob.cl.med)
colnames(patch.jsd) <- rownames((patch.prob.cl.med))

reorder.cl <- sort(match(rownames(patch.prob.cl.med), labels(dend)))
row.and.colnames <- labels(dend)[reorder.cl]
patch.jsd <- patch.jsd[row.and.colnames, row.and.colnames]


ggplot(data = melt(patch.jsd), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()+ theme(axis.text = element_text(size=7)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("") +
  scale_fill_gradient(low = "blue", high = "red")

#Lets sort each row:
sorted_patch.jsd.cl <- t(apply(patch.jsd, 1, function(x) names(x)[order(x, decreasing = FALSE)]))
sorted_patch.jsd.value <- t(apply(patch.jsd, 1, function(x) x[order(x, decreasing = FALSE)]))

i = 48 ; head(sorted_patch.jsd.cl[i,]) ; head(sorted_patch.jsd.value[i,])

######################################################################################################
### Does JSD correspond with how mapping confusion looks like? #######################################
######################################################################################################
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/mapping.memb.rda")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/REF_mapping_probability.rda")

ggplot(data = melt(Tree_mapping_probability[33:93,33:93]), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Real label") + ylab("Mapped labels") +  
  scale_fill_gradient(low = "dark blue", high = "yellow")

tail(sort(Tree_mapping_probability[81, ]))



