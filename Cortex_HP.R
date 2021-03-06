#Mapping whole cortex:
#source("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20190125_collapsed40_cpm/patchseq/patchseq.R")
#source("//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/patchseq/patchseq.R")
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/patchseq_10X.R")
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/dend.markers_10X.R")
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/")
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/patchseq/de.genes.R")
source("//allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/scrattch.hicat/R/markers.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")
.libPaths(c("/allen/programs/celltypes/workgroups/rnaseqanalysis/Script_Repository/Olivia/3.5", .libPaths()))
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/patchseq/Mapping_helper_functions.R")

library(scrattch.hicat)
library(scrattch.io)
library(scrattch.vis)
library(dendextend)
library(matrixStats)
library(feather)
library(Matrix)

path.cortex.hp <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/Cortex_HIP/Cortex_HIP/"
path.anno <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/Cortex_HIP/Cortex_HIP/correct/" 
work.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/smartseq_10X_combined"

##########################################################################################
### Reading FACS calls: ##################################################################
##########################################################################################
# res.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"
# tmp.load1 = load(file=file.path(res.dir, "ref.data.rda")) 
# 
# rownames(cl.df)=cl.df$cluster_id
# cltmp=as.data.frame(cl.df[as.character(cl),c("cluster_label", "cluster_color", "cluster_id")])
# rownames(cltmp)=paste0("Smartseq_cells.", names(cl))
# FACS <- cltmp
# #cl=factor(cltmp)
# #FACS.cl <- cl
# #FACS.cl.df <- cl.df
# #FACS.dend <- dend
# rm(cl)
# rm(dend)
# rm(cl.df)


##########################################################################################
### Comparing FACS with concensus clustering: ############################################
##########################################################################################
#cl file
load(paste0(path.anno, "cl.final.rda"))
length(cl)

#anno file 
load(paste0(path.anno, "anno.df.rda"))
dim(anno.df)
anno.df <- as.data.frame(anno.df)
rownames(anno.df) <- anno.df$sample_name

# FACS.cells <- names(FACS.cl)
# FACS.cells <- paste0("Smartseq_cells.", FACS.cells)
# All.smartseq.cells <- names(cl)[grepl("Smartseq", names(cl))]
# 
# FACS.cells <- FACS.cells[FACS.cells %in% All.smartseq.cells]
# rownames(cl.df) = cl.df$cl
# cltmp = as.data.frame(cl.df[as.character(cl),c("cluster_label", "cluster_id", "cluster_color")])
# rownames(cltmp) = names(cl)
# Smartseq.FACS <- cltmp[FACS.cells,]
# FACS <- FACS[FACS.cells,]
# colnames(FACS) <- c("old_cluster_label", "old_cluster_color", "old_cluster_id")
# colnames(Smartseq.FACS) <- c("new_cluster_label","new_cluster_id", "new_cluster_color")
# results <- cbind.data.frame(FACS, Smartseq.FACS)
# anno.df <- anno.df[FACS.cells, ]
# anno.df <- anno.df[anno.df$class_label == "GABAergic",]
# GABA.cells <- rownames(anno.df)
# unique(anno.df[GABA.cells, "subclass_label"])
# Pvalb_Sst_cells <- GABA.cells[anno.df$subclass_label %in% c("Pvalb", "Sst", "Sst Chodl")]
# Vip_Lamp5_cells <- GABA.cells[!anno.df$subclass_label %in% c("Pvalb", "Sst", "Sst Chodl")]
# 
# dim(results)
# source("/allen//programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Plot_utils.R")
# library(dplyr)
# 
# 
# tmp <- results[Vip_Lamp5_cells,] %>% 
#   rownames_to_column("id") %>% 
#   dplyr::select(id, old_cluster_label, old_cluster_color, old_cluster_id, 
#                 new_cluster_label, new_cluster_id, new_cluster_color) %>%
#   `colnames<-` (c("sample_id", "map_cluster_label", "map_cluster_color", 
#                   "map_cluster_id", "cluster_label", "cluster_id", "cluster_color"))
# 
# river_plot(tmp, min.cells=0, min.frac=0)

##########################################################################################
### Reading inputs: ######################################################################
##########################################################################################
#Downsampled 10X and smartseq norm.data 
load(paste0(path.cortex.hp, "ref.dat.list.rda"))
dim(ref.dat.list$`10X_cells`)
dim(ref.dat.list$Smartseq_cells)
norm.dat <- cbind(ref.dat.list$Smartseq_cells, ref.dat.list$`10X_cells`)

Smartseq_cells <- colnames(norm.dat)[grepl("Smartseq_cells", colnames(norm.dat))]
TenX_cells <- colnames(norm.dat)[grepl("10X_cells", colnames(norm.dat))]

#DE genes for both 10X and smartseq
load(paste0(path.cortex.hp, "comb.de.genes.rda"))
length(comb.de.genes)

#markers
load(paste0(path.cortex.hp, "select.markers.rda"))
length(select.markers)

#dend file
#load(paste0(path.anno, "dend.rda"))
#plot(dend)

#cluster centroid
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/joint_analysis/Cortex_HIP/Cortex_HIP/cl.means.list.rda")
cl.means_10X <- cl.means.list$`10X_cells`
cl.means_smartseq <- cl.means.list$Smartseq_cells

##########################################################################################
### Some initialization: #################################################################
##########################################################################################

TenX.norm.dat <- ref.dat.list$`10X_cells`
Smartseq.norm.dat <- ref.dat.list$Smartseq_cells

cl_10X_only <- cl[TenX_cells]
cl_10X_only <- droplevels(cl_10X_only)

cl_smartseq_only <- cl[Smartseq_cells]
cl_smartseq_only <- droplevels(cl_smartseq_only)

##########################################################################################
### Adding markers to the dend: ##########################################################
##########################################################################################
#We prune the dend based on 10X data
missing_cl <- names(table(cl))[!names(table(cl)) %in% names(table(cl_10X_only))] 

#Also all the types have less than 10cells should be removed???
dend_pruned <- prune_dend(dend = dend, rm.labels = missing_cl)

#We are going to add the markers based on both 10X and smartseq
cl <- cl[colnames(norm.dat)]
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
save(dend_pruned, file =paste0(work.dir, "/dend_pruned_markers_attached.rda"))
#load(paste0(work.dir, "/dend_pruned_markers_attached.rda"))
#save(dend, file="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2/Fahimeh_dend_markers_attached.rda")

##########################################################################################
### Loading the prunded dend: ############################################################
##########################################################################################

#load(paste0(work.dir, "/dend_pruned_markers_attached.rda"))

##########################################################################################
### Initialization: ######################################################################
##########################################################################################

select.cl <- labels(dend_pruned)
TenX_cells <- names(cl_10X_only[cl_10X_only %in% labels(dend_pruned)]) 

cl_10X_only <- cl_10X_only[TenX_cells]
cl_10X_only <- droplevels(cl_10X_only)

Smartseq_cells <- names(cl_smartseq_only[cl_smartseq_only %in% labels(dend_pruned)])
cl_smartseq_only <- cl_smartseq_only[Smartseq_cells]
cl_smartseq_only <- droplevels(cl_smartseq_only)

#write.csv(as.matrix(TenX.norm.dat[select.markers,]), file = paste0(work.dir, "/TenX_norm_dat.csv"))
#source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/patchseq_10X.R")


##########################################################################################
### Mapping 10X on 10X: ##################################################################
##########################################################################################
TenX_test_cells <- sample(TenX_cells, 35000)
TenX_train_cells <- setdiff(TenX_cells, TenX_test_cells)

train.cl.med= do.call("cbind",tapply(names(cl_10X_only[TenX_test_cells]), cl_10X_only[TenX_test_cells], function(x){
  rowMeans(as.matrix(norm.dat[,x,drop=F]))
}))

memb_tenx2tenx = map_dend_membership(dend_pruned, cl=cl_10X_only[TenX_train_cells], 
                                     TenX.norm.dat[select.markers, TenX_train_cells], 
                                     TenX.norm.dat[select.markers, TenX_test_cells], 
                                     TenX_test_cells, bs.num=100, 
                                     p=0.7, low.th=0.15, 
                                     cl.mean = cl.means_10X[select.markers, select.cl])

save(memb_tenx2tenx, file = paste0(work.dir, "/memb_tenx2tenx.rda"))


##########################################################################################
### 10X on 10X mapping performance: ######################################################
##########################################################################################
rownames(anno.df) <- anno.df$sample_name
anno.df <- as.data.frame(anno.df)

labels(dend_pruned)[!labels(dend_pruned) %in% colnames(memb_tenx2tenx)]
cl.df[cl.df$cl == labels(dend_pruned)[!labels(dend_pruned) %in% colnames(memb_tenx2tenx)], "cluster_label"]
mapped.cl <- colnames(memb_tenx2tenx)[colnames(memb_tenx2tenx) %in% labels(dend_pruned)]

sorted_cl <- t(apply(memb_tenx2tenx[, mapped.cl], 1, function(x) names(x)[order(x, decreasing = TRUE)]))
sorted_bt <- t(apply(memb_tenx2tenx[, mapped.cl], 1, function(x) x[order(x, decreasing = TRUE)]))

sum(sorted_cl[TenX_test_cells,1] == cl_10X_only[TenX_test_cells]) / length(TenX_test_cells)

Confusion_mat <- table(sorted_cl[TenX_test_cells,1], cl_10X_only[TenX_test_cells])
Confusion_mat <- Confusion_mat/rowMaxs(Confusion_mat)
  
library(reshape2)
ggplot(data = melt(Confusion_mat), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Clustering lable") + ylab("Mapping lables") +  
  scale_fill_gradient(low = "white", high = "red")

compare_annotate(cl=sorted_cl[TenX_test_cells,1], ref.cl = cl_10X_only[TenX_test_cells],ref.cl.df = cl.df , reorder = TRUE)

##########################################################################################
### Mapping Smartseq using 10X as ref: ###################################################
##########################################################################################

memb_smart2tenx = map_dend_membership(dend_pruned, cl=cl_10X_only, 
                                      TenX.norm.dat[select.markers,], 
                                      Smartseq.norm.dat[select.markers,], 
                                      Smartseq_cells, bs.num=100, 
                                      p=0.7, low.th=0.15, 
                                      cl.mean = cl.means_10X[select.markers, select.cl])

save(memb_smart2tenx, file = paste0(work.dir, "/memb_smart2tenx.rda"))
#No cell was mapped to this cluster:
labels(dend_pruned)[!labels(dend_pruned) %in% colnames(memb_smart2tenx)]
cl.df[cl.df$cl == labels(dend_pruned)[!labels(dend_pruned) %in% colnames(memb_smart2tenx)], "cluster_label"]
mapped.cl <- colnames(memb_smart2tenx)[colnames(memb_smart2tenx) %in% labels(dend_pruned)]

#rownames(anno.df) <- anno.df$sample_name
#anno.df <- as.data.frame(anno.df)
sorted_cl <- t(apply(memb_smart2tenx[, mapped.cl], 1, function(x) names(x)[order(x, decreasing = TRUE)]))
sorted_bt <- t(apply(memb_smart2tenx[, mapped.cl], 1, function(x) x[order(x, decreasing = TRUE)]))

sum(as.character(anno.df[Smartseq_cells, "cl"]) == as.character(sorted_cl[Smartseq_cells,1])) / length(Smartseq_cells)
sum(as.character(anno.df[Smartseq_cells, "cl"]) == as.character(sorted_cl[Smartseq_cells,2])) / length(Smartseq_cells)

compare_annotate(cl=sorted_cl[Smartseq_cells,1], ref.cl = cl_smartseq_only[Smartseq_cells],ref.cl.df = cl.df , reorder = TRUE)

##########################################################################################
### Region comparison: ###################################################################
##########################################################################################

TenX_regions_freq <- as.data.frame(table(anno.df[anno.df$sample_name %in% TenX_cells,"joint_region_label"]))
TenX_regions_freq$platform <- "10X"

Smartseq_regions_freq <- as.data.frame(table(anno.df[anno.df$sample_name %in% Smartseq_cells,"joint_region_label"]))
Smartseq_regions_freq$platform <- "Smartseq"

left_join(TenX_regions_freq, Smartseq_regions_freq, by= "Var1")

tmp <- rbind.data.frame(TenX_regions_freq, Smartseq_regions_freq)

ggplot(tmp,aes(x=Var1,y=Freq,fill=platform))+
  geom_bar(stat="identity",position="dodge") +
  xlab("Region")+ylab("Number of cells") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



##########################################################################################
### cor: #################################################################################
##########################################################################################
anno.df <- as.data.frame(anno.df)
rownames(anno.df) <- anno.df$sample_name
anno.df_TenX <- anno.df[TenX_cells,]
anno.df_TenX <- anno.df_TenX[anno.df_TenX$cl %in% labels(dend_pruned),]

cl.mean <- do.call("cbind", tapply(anno.df_TenX$sample_name, 
                                   anno.df_TenX$cl, function(x)rowMeans(TenX.norm.dat[select.markers, x, drop=F])))

save(cl.mean, file=paste0(work.dir, "/cl.mean.rda"))
load(paste0(work.dir, "/cl.mean.rda"))

matrix_cor <- cor(cl.mean[select.markers,], cl.mean[select.markers,])
sorted_cl <- as.character(sort(as.numeric(rownames(matrix_cor))))
matrix_cor <- matrix_cor[sorted_cl, sorted_cl]
save(matrix_cor, file=paste0(work.dir, "/matrix_cor.rda"))

library(reshape2)
ggplot(data = melt(matrix_cor), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Clustering lable") + ylab("NN Mapping lables") +  
  scale_fill_gradient(low = "white", high = "red")

##########################################################################################
### JSD: #################################################################################
##########################################################################################
library("philentropy")
t_prob <- TenX.norm.dat
#norm.dat is log2(cpm+1) and we want the cpm data
t_prob@x <- (2 ^ TenX.norm.dat@x) - 1
sf = Matrix::colSums(t_prob)
sep = t_prob@p
sep = sep[-1] - sep[-length(sep)]
j = S4Vectors::Rle(1:length(sep), sep)
t_prob@x = t_prob@x/sf[as.integer(j)]

cl.t.prob.mean <- do.call("cbind", tapply(anno.df_TenX$sample_name, 
                                          anno.df_TenX$cl, function(x)rowMeans(t_prob[, x, drop=F])))

save(cl.t.prob.mean, file=paste0(work.dir, "/cl.t.prob.mean.rda"))
load(paste0(work.dir, "/cl.t.prob.mean.rda"))
dim(cl.t.prob.mean)

matrix_JSD <- JSD(cl.t.prob.mean[,select.markers],  test.na = TRUE, unit = "log2", est.prob = NULL)
rownames(matrix_JSD) <- rownames(cl.t.prob.mean)
colnames(matrix_JSD) <- rownames(cl.t.prob.mean)
matrix_JSD <- matrix_JSD[sorted_cl, sorted_cl]
save(matrix_JSD, file=paste0(work.dir, "/matrix_JSD_markergenes.rda"))

library(reshape2)
ggplot(data = melt(matrix_JSD), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Clustering lable") + ylab("NN Mapping lables") +  
  scale_fill_gradient(low = "white", high = "red")

load(paste0(work.dir, "/matrix_JSD_markergenes.rda"))
load(paste0(work.dir, "/matrix_cor.rda"))


temp <-  matrix_cor[upper.tri(matrix_cor, diag = FALSE)]
tt1 <- as.vector(temp)
sorted_tt1 <- sort(tt1, decreasing = TRUE)
tt2 <- sapply(tt1, function(x)match(x, sorted_tt1))
ranked_cor_mat <- matrix(0, nrow = 383, ncol = 383)
ranked_cor_mat[upper.tri(ranked_cor_mat, diag = FALSE)] <- tt2
tt3 <- t(ranked_cor_mat)
ranked_cor_mat[lower.tri(ranked_cor_mat, diag = FALSE)] <- tt3[lower.tri(tt3, diag = FALSE)]


temp <- matrix_JSD[upper.tri(matrix_JSD, diag = FALSE)]
tt1 <- as.vector(temp)
sorted_tt1 <- sort(tt1, decreasing = FALSE)
tt2 <- sapply(tt1, function(x)match(x, sorted_tt1))
ranked_JSD_mat <- matrix(0, nrow = 383, ncol = 383)
ranked_JSD_mat[upper.tri(ranked_JSD_mat, diag = FALSE)] <- tt2
tt3 <- t(ranked_JSD_mat)
ranked_JSD_mat[lower.tri(ranked_JSD_mat, diag = FALSE)] <- tt3[lower.tri(tt3, diag = FALSE)]
rownames(ranked_JSD_mat) <- seq(1:383)
colnames(ranked_JSD_mat) <- seq(1:383)

dist <- abs(ranked_JSD_mat-ranked_cor_mat)
norm.dist <- dist/max(dist)

ggplot(data = melt(norm.dist), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("") + ylab("") +  
  scale_fill_gradient(low = "white", high = "red")


##########################################################################################
### Mapping smartseq on smartseq: ########################################################
##########################################################################################

Smartseq_test_cells <- sample(Smartseq.cells, 10000)
Smartseq_train_cells <- setdiff(Smartseq_test_cells, TenX_test_cells)

memb_smart2smart = map_dend_membership(dend_pruned, cl=cl_smartseq_only[Smartseq_train_cells], 
                                     Smartseq.norm.dat[select.markers, Smartseq_train_cells], 
                                     Smartseq.norm.dat[select.markers, Smartseq_test_cells], 
                                     Smartseq_test_cells, bs.num=100, 
                                     p=0.7, low.th=0.15, 
                                     cl.mean = cl.means_smartseq[select.markers, ])

save(memb_smart2smart, file = paste0(work.dir, "/memb_smart2smart.rda"))

##########################################################################################
### Mapping Patchseq on 10X: #############################################################
##########################################################################################

batch_date="20190227_BT014-RSC-195"
patchseq.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/"
robject.dir = "//allen//programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
query.dat.norm = Get_patchseq_norm(cpm_rdata_path = paste0(robject.dir,batch_date,"_mouse_patchseq_star2.0_cpm.Rdata"),
                                   samp_rdata_path = paste0(robject.dir,batch_date,"_mouse_patchseq_star2.0_samp.dat.Rdata"))
patchseq_anno <- Read_patchseq_anno(paste0(patchseq.dir,"anno.feather"))
patchseq.cells <- colnames(query.dat.norm)
select.cl <- labels(dend_pruned)

memb_patchseq2tenx= map_dend_membership(dend_pruned, cl=cl_10X_only[TenX.cells], 
                                       TenX.norm.dat[select.markers, TenX.cells], 
                                       query.dat.norm[select.markers, patchseq.cells], 
                                       patchseq.cells, bs.num=100, 
                                       p=0.7, low.th=0.15, 
                                       cl.mean = cl.means_10X[select.markers, select.cl])

#save(memb_patchseq2tenx, file = paste0(work.dir, "/memb_patchseq2tenx.rda"))

mapped.cl <- colnames(memb_patchseq2tenx)[colnames(memb_patchseq2tenx) %in% labels(dend_pruned)]

sorted_cl <- t(apply(memb_patchseq2tenx[, mapped.cl], 1, function(x) names(x)[order(x, decreasing = TRUE)]))
sorted_bt <- t(apply(memb_patchseq2tenx[, mapped.cl], 1, function(x) x[order(x, decreasing = TRUE)]))

first_cl <- sorted_cl[,1]
second_cl <- sorted_cl[,2]
first_cluster_label <- Renew_list(first_cl, ref.df = cl.df, label = "cl", new.label = "cluster_label")
second_cluster_label <- Renew_list(second_cl, ref.df = cl.df, label = "cl", new.label = "cluster_label")

sum(as.character(patchseq_anno[patchseq.cells, "cluster_label"]) == as.character(first_cluster_label[patchseq.cells])) / length(patchseq.cells)
sum(as.character(anno.df[, "cluster_label"]) == as.character(second_cluster_label[Smartseq.cells])) / length(Smartseq.cells)
sort(table(first_cluster_label))


##########################################################################################
### River Plots: #########################################################################
##########################################################################################

source("/allen//programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Plot_utils.R")
library(dplyr)

select.cells <- Smartseq.cells[Smartseq.cells %in% rownames(anno.df)[anno.df$class_label == "GABAergic"]]
select.cells <- Smartseq.cells[Smartseq.cells %in% rownames(anno.df)[anno.df$class_label == "Glutamatergic"]]

tmp <- cbind.data.frame(select.cells , first_cluster_label[select.cells], as.character(anno.df[select.cells, "cluster_label"]))
colnames(tmp) <- c("sample_id","cluster_label", "cluster_label_anno")
tmp <- left_join(tmp, cl.df) %>% select(sample_id, cluster_label, cluster_id, cluster_color, cluster_label_anno)
colnames(tmp) <- c("sample_id", "map_cluster_label", "map_cluster_id", "map_cluster_color", "cluster_label")

tmp <- left_join(tmp, cl.df) %>% select(sample_id, map_cluster_label, map_cluster_id, map_cluster_color, 
                                        cluster_label, cluster_id, cluster_color)
colnames(tmp) <- c( "sample_id", "map_cluster_label", "map_cluster_id", "map_cluster_color",
                    "cluster_label", "cluster_id", "cluster_color")

river_plot(tmp, min.cells=0, min.frac=0)

##########################################################################################
### Bar Plots: ###########################################################################
##########################################################################################

library(reshape2)
mapping <- as.data.frame(table(first_cl))
colnames(mapping) <- c("cluster_label", "mapping")
clustering <- as.data.frame(table(anno.df[Smartseq.cells, "cl"])) 
colnames(clustering) <- c("cluster_label", "clustering")
tmp <- left_join(mapping, clustering, by= "cluster_label")
tmp <- melt(tmp)
colnames(tmp) <- c("cluster_label", "Method", "value" )

ggplot(tmp,aes(x=cluster_label,y=value,fill=Method))+
  geom_bar(stat="identity",position="dodge") +
  xlab("Cl")+ylab("Size") 
