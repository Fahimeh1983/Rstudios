library(feather)
library(dplyr)
library(matrixStats)
library(scrattch.hicat)
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/patchseq_10X.R")
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/dend.markers_10X.R")
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/")
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/patchseq/de.genes.R")
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/scrattch.hicat/R/markers.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")
.libPaths(c("/allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Olivia/3.5", .libPaths()))
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/patchseq/Mapping_helper_functions.R")

##########################################################################################
### Reading inputs: ######################################################################
##########################################################################################
path.cortex.hp <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/Cortex_HIP/Cortex_HIP/"
path.anno <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/Cortex_HIP/Cortex_HIP/correct/" 
work.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/smartseq_10X_combined"

#Downsampled 10X and smartseq norm.data 
load(paste0(path.cortex.hp, "ref.dat.list.rda"))
dim(ref.dat.list$`10X_cells`)
dim(ref.dat.list$Smartseq_cells)
norm.dat <- cbind(ref.dat.list$Smartseq_cells, ref.dat.list$`10X_cells`)

Smartseq_cells <- colnames(norm.dat)[grepl("Smartseq_cells", colnames(norm.dat))]
TenX_cells <- colnames(norm.dat)[grepl("10X_cells", colnames(norm.dat))]

#markers
load(paste0(path.cortex.hp, "select.markers.rda"))
length(select.markers)

#anno file 
load(paste0(path.anno, "anno.df.rda"))
dim(anno.df)
anno.df <- as.data.frame(anno.df)
rownames(anno.df) <- anno.df$sample_name

#cl file
load(paste0(path.anno, "cl.final.rda"))
length(cl)

#DE genes for both 10X and smartseq
load(paste0(path.cortex.hp, "comb.de.genes.rda"))
length(comb.de.genes)

#dend file
load(paste0(path.anno, "dend.rda"))
#plot(dend)

#cluster centroid
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/Cortex_HIP/Cortex_HIP/cl.means.list.rda")
cl.means_10X <- cl.means.list$`10X_cells`
cl.means_smartseq <- cl.means.list$Smartseq_cells

#load(paste0(work.dir, "/dend_pruned_markers_attached.rda"))

##########################################################################################
### Adding markers to the dend which is based on 10X types:: #############################
##########################################################################################
#Rm low quality cells
unique(anno.df[TenX_cells, "class_label"])
unique(anno.df[Smartseq_cells, "class_label"])
LQ_TenX_cells <- TenX_cells[anno.df[TenX_cells, "class_label"] =="Low Quality"]
LQ_Smartseq_cells <- Smartseq_cells[anno.df[Smartseq_cells, "class_label"] == "Low Quality"]

#Number of types in each platform with more than 10 cells
sum(sort(table(cl[setdiff(Smartseq_cells, LQ_Smartseq_cells)]), decreasing = TRUE) > 10)
sum(sort(table(cl[setdiff(TenX_cells, LQ_TenX_cells)]), decreasing = TRUE) > 10)

#preparing 10X training data
train_cl <- names(table(cl))[table(cl[setdiff(TenX_cells, LQ_TenX_cells)]) >= 10]
missing_cl <- setdiff(names(table(cl)), train_cl)
dend_pruned <- prune_dend(dend = dend, rm.labels = missing_cl)

#We are going to add the markers based on both 10X and smartseq
cl <- cl[c(setdiff(Smartseq_cells, LQ_Smartseq_cells), setdiff(TenX_cells, LQ_TenX_cells))]
cl <- cl[cl %in% labels(dend_pruned)]
cl <- droplevels(cl)

tmp <- get_gene_score(comb.de.genes, top.n = 100, max.num = 2000, bin.th = 4)
up.gene.score <- tmp$up.gene.score
down.gene.score <- tmp$down.gene.score
select.genes <- row.names(up.gene.score)
n.markers = 30

dend_pruned <- select_dend_markers(dend_pruned, norm.dat = norm.dat, cl = cl,
                                   de.genes = comb.de.genes, up.gene.score = up.gene.score[select.genes,],
                                   down.gene.score = down.gene.score[select.genes,], n.markers = n.markers)
temp <- dend_pruned
dend_pruned = select_pos_dend_markers(dend= dend_pruned, norm.dat = norm.dat, cl = cl)


#save(dend, file =paste0(work.dir, "dend_markers_attached.rda"))
save(dend_pruned, file =paste0(work.dir, "/dend_pruned_markers_attached_373types.rda"))
#load(paste0(work.dir, "/dend_pruned_markers_attached.rda"))
#save(dend, file="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2/Fahimeh_dend_markers_attached.rda")


##########################################################################################
### Adding markers to the dend which is based on smartseq types: #########################
##########################################################################################
#Rm low quality cells
unique(anno.df[TenX_cells, "class_label"])
unique(anno.df[Smartseq_cells, "class_label"])
LQ_TenX_cells <- TenX_cells[anno.df[TenX_cells, "class_label"] =="Low Quality"]
LQ_Smartseq_cells <- Smartseq_cells[anno.df[Smartseq_cells, "class_label"] == "Low Quality"]

#Number of types in each platform with more than 10 cells
sum(sort(table(cl[setdiff(Smartseq_cells, LQ_Smartseq_cells)]), decreasing = TRUE) >= 10)
sum(sort(table(cl[setdiff(TenX_cells, LQ_TenX_cells)]), decreasing = TRUE) > 10)

#preparing 10X training data
train_cl <- names(table(cl))[table(cl[setdiff(Smartseq_cells, LQ_Smartseq_cells)]) >= 10]
missing_cl <- setdiff(names(table(cl)), train_cl)
dend_pruned <- prune_dend(dend = dend, rm.labels = missing_cl)

#We are going to add the markers based on both 10X and smartseq
cl <- cl[c(setdiff(Smartseq_cells, LQ_Smartseq_cells), setdiff(TenX_cells, LQ_TenX_cells))]
cl <- cl[cl %in% labels(dend_pruned)]
cl <- droplevels(cl)

tmp <- get_gene_score(comb.de.genes, top.n = 100, max.num = 2000, bin.th = 4)
up.gene.score <- tmp$up.gene.score
down.gene.score <- tmp$down.gene.score
select.genes <- row.names(up.gene.score)
n.markers = 30

dend_pruned <- select_dend_markers(dend_pruned, norm.dat = norm.dat, cl = cl,
                                   de.genes = comb.de.genes, up.gene.score = up.gene.score[select.genes,],
                                   down.gene.score = down.gene.score[select.genes,], n.markers = n.markers)
temp <- dend_pruned
dend_pruned = select_pos_dend_markers(dend= dend_pruned, norm.dat = norm.dat, cl = cl)


#save(dend, file =paste0(work.dir, "dend_markers_attached.rda"))
save(dend_pruned, file =paste0(work.dir, "/dend_pruned_markers_attached_314types.rda"))
#load(paste0(work.dir, "/dend_pruned_markers_attached.rda"))
#save(dend, file="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2/Fahimeh_dend_markers_attached.rda")

