library(tibble)
library(feather)
library(dplyr)
library(matrixStats)
library(scrattch.hicat)
library(cowplot)
library(feather)
library(dendextend)
library(scrattch.vis)
library(scrattch.io)
library("cowplot")
library(ggbeeswarm)
.libPaths("/home/fahimehb/R/x86_64-redhat-linux-gnu-library/3.5")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")


##########################################################################################
### Functions: ###########################################################################
##########################################################################################

flat.mapping <- function(test.dat, train.cl.med, train.markers, sample.markers.prop=0.8, iter=100)
{
  train.markers=intersect(train.markers, row.names(test.dat))
  
  map.freq = sapply(1:iter, function(i){
    tmp.markers=sample(train.markers, round(length(train.markers)*sample.markers.prop))
    test.cl.cor = cor(as.matrix(test.dat[tmp.markers,]), train.cl.med[tmp.markers,])
    test.cl.cor[is.na(test.cl.cor)]=0
    best.cl= setNames(colnames(test.cl.cor)[apply(test.cl.cor, 1, which.max)], row.names(test.cl.cor))
  })
  tmp=  as.data.frame(as.table(as.matrix(map.freq)))
  map.freq <- as.matrix(table(tmp$Var1, tmp$Freq))/iter
  
  return(map.freq)
}

##########################################################################################
### Reading inputs: ######################################################################
##########################################################################################

path.cortex.hp <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/Cortex_HIP/Cortex_HIP/"
path.anno <- "//allen/programs/celltypes/workgroups/mct-t200/Manuscript/2019_Yao/common/" 
work.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/smartseq_10X_combined"

#markers
load(paste0(path.cortex.hp, "select.markers.rda"))
length(select.markers)


#cl file
load(paste0(path.anno, "cl.final.rda"))
length(cl)
cl.consensus <- cl
cl.df.consensus <- cl.df
rm(cl)
rm(cl.df)

#cluster centroid
load(paste0(path.anno,"/cl.means.list.rda"))
cl.means_10X <- cl.means.list$`10X_cells`
cl.means_smartseq <- cl.means.list$Smartseq_cells
dim(cl.means_10X)
dim(cl.means_smartseq)

consensus_types <- union(colnames(cl.means_smartseq), colnames(cl.means_10X))
length(consensus_types)

##########################################################################################
### Mapping Cembrwoski dataset on smartseq ref: ###############################################
##########################################################################################
Cembrowski.folder <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/smartseq_10X_combined/public_data/Cembrowski/"

select.cl <- consensus_types
anno <- read.table(paste0(Cembrowski.folder, "/idents.txt"), sep = "")
data <- read.table(paste0(Cembrowski.folder, "/data.txt"), sep = "")

dim(data)
colnames(data)

Genes <- rownames(data)[rownames(data) %in% select.markers]
sum(select.markers %in% rownames(data))
Cemb.cells <- colnames(data)


#library(scrattch.io)
data <- as.matrix(data)
cpm.data <- logCPM(data)
cpm.data <- as.data.frame(cpm.data)
dim(cl.means_10X)
dim(cpm.data)


select.cl <- colnames(cl.means_smartseq)

memb <- flat.mapping(test.dat = cpm.data, 
                     train.cl.med = cl.means_smartseq[, select.cl], 
                     train.markers = select.markers, 
                     sample.markers.prop = 0.8, 
                     iter = 100)

load(paste0(Cembrowski.folder, "/memb.Cembrowski.rda"))

sort.cl <- t(apply(memb, 1, function(x) names(x)[order(x, decreasing = TRUE)]))
best.cl <- sort.cl[,1]
best.label <- Renew_list(ls = best.cl, ref.df = cl.df.consensus, label = "cl", new.label = "cluster_label")

# Remove those with less than 10 cells
small.cls <- names(table(best.cl))[(table(best.cl) < 10)]
best.cl <- best.cl[!best.cl %in% small.cls]
best.cl.cells <- names(best.cl)
best.label <- best.label[best.cl.cells]

Cemb.label <- setNames(anno$x, rownames(anno))
Cemb.label <- Cemb.label[best.cl.cells]

save(Cemb.label, file=paste0(Cembrowski.folder, "/Cemb.label.rda"))
save(best.cl, file=paste0(Cembrowski.folder, "/best.cl.Cembrowski.rda"))

levels <- cl.df.consensus$cl
rownames(cl.df.consensus) <-cl.df.consensus$cl
best.cl = droplevels(factor(best.cl, levels = levels))

Cembrowski.results <- compare_annotate(cl= Cemb.label, 
                                       ref.cl = best.cl, 
                                       ref.cl.df = cl.df.consensus, 
                                       reorder = TRUE, 
                                       rename = FALSE)

#save(Cembrowski.results, file=paste0(Cembrowski.folder, "/Cembrowski.results.rda"))
load(paste0(Cembrowski.folder, "/Cembrowski.results.rda"))