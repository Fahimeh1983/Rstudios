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

#DE genes for both 10X and smartseq
load(paste0(path.cortex.hp, "comb.de.genes.rda"))
length(comb.de.genes)

#markers
load(paste0(path.cortex.hp, "select.markers.rda"))
length(select.markers)

#anno file 
load(paste0(path.anno, "anno.df.rda"))
dim(anno.df)

#cl file
load(paste0(path.anno, "cl.final.rda"))
length(cl)


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

TenX.cells <- colnames(TenX.norm.dat)
Smartseq.cells <- colnames(Smartseq.norm.dat)

cl_10X_only <- cl[TenX.cells]
cl_10X_only <- droplevels(cl_10X_only)

cl_smartseq_only <- cl[Smartseq.cells]
cl_smartseq_only <- droplevels(cl_smartseq_only)

##########################################################################################
### cor: #################################################################################
##########################################################################################
anno.df <- as.data.frame(anno.df)
rownames(anno.df) <- anno.df$sample_name
anno.df_TenX <- anno.df[TenX.cells,]
tmp <- unique(anno.df_TenX$cl)
labels <- as.character(Renew_list(tmp, ref.df = cl.df, label = "cl", new.label = "cluster_label"))
tmp <- setdiff(tmp, tmp[which(is.na(labels))])
labels <- setdiff(labels, labels[which(is.na(labels))])
ref.labels <- setNames(labels, tmp)

cl.mean <- do.call("cbind", tapply(anno.df_TenX$sample_name, 
                                   anno.df_TenX$cl, function(x)rowMeans(TenX.norm.dat[select.markers, x, drop=F])))
matrix_cor <- cor(cl.mean[select.markers,], cl.mean[select.markers,])
ref.labels[rownames(matrix_cor)]
rownames(matrix_cor) <- ref.labels[rownames(matrix_cor)]
colnames(matrix_cor) <- ref.labels[colnames(matrix_cor)]
organized.cl <- cl.df$cluster_label
organized.labels <- as.character(organized.cl[organized.cl %in% labels])
matrix_cor <- matrix_cor[organized.labels, organized.labels]
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


#t_prob <- as.matrix(t_prob)
anno.df <- as.data.frame(anno.df)
rownames(anno.df) <- anno.df$sample_name
anno.df_TenX <- anno.df[TenX.cells,]


cl.t.prob.mean <- do.call("cbind", tapply(anno.df_TenX$sample_name, 
                                          anno.df_TenX$cl, function(x)rowMeans(t_prob[, x, drop=F])))
cl.t.prob.mean <- t(cl.t.prob.mean)

matrix_JSD <- JSD(cl.t.prob.mean,  test.na = TRUE, unit = "log2", est.prob = NULL)
rownames(matrix_JSD) <- ref.labels[rownames(cl.t.prob.mean)] 
colnames(matrix_JSD) <- ref.labels[rownames(cl.t.prob.mean)]
matrix_JSD <- matrix_JSD[organized.labels, organized.labels]
save(matrix_JSD, file=paste0(work.dir, "/matrix_JSD_allgenes.rda"))

library(reshape2)
ggplot(data = melt(matrix_JSD), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Clustering lable") + ylab("NN Mapping lables") +  
  scale_fill_gradient(low = "white", high = "red")

load(paste0(work.dir, "/matrix_JSD_allgenes.rda"))
matrix_JSD_all <- matrix_JSD
rm(matrix_JSD)
load(paste0(work.dir, "/matrix_JSD_markergenes.rda"))
matrix_JSD_markers <- matrix_JSD
rm(matrix_JSD)
load(paste0(work.dir, "/matrix_cor.rda"))

sorted_cor <- t(apply(matrix_cor, 1, function(x)x[order(x, decreasing = TRUE)]))
sorted_cor_cl <- t(apply(matrix_cor, 1, function(x)names(x)[order(x, decreasing = TRUE)]))

sorted_JSD_all <- t(apply(matrix_JSD_all, 1, function(x)x[order(x, decreasing = FALSE)]))
sorted_JSD_all_cl <- t(apply(matrix_JSD_all, 1, function(x)names(x)[order(x, decreasing = FALSE)]))

sorted_JSD_markers <- t(apply(matrix_JSD_markers, 1, function(x)x[order(x, decreasing = FALSE)]))
sorted_JSD_markers_cl <- t(apply(matrix_JSD_markers, 1, function(x)names(x)[order(x, decreasing = FALSE)]))

results <- cbind.data.frame(sorted_cor[,c(1,2,3)], sorted_cor_cl[,c(1,2,3)], 
                            sorted_JSD_markers[,c(1,2,3)], sorted_JSD_markers_cl[,c(1,2,3)],
                            sorted_JSD_all[,c(1,2,3)], sorted_JSD_all_cl[,c(1,2,3)])

colnames(results) <- c("cor_score1", "cor_score2", "cor_score3",
                       "cor_cl1", "cor_cl2", "cor_cl3",
                       "JSD_markers_score1", "JSD_markers_score2",  "JSD_markers_score3",
                       "JSD_markers_cl1", "JSD_markers_cl2", "JSD_markers_cl3",
                       "JSD_all_score1", "JSD_all_score2",  "JSD_all_score3",
                       "JSD_all_cl1", "JSD_all_cl2", "JSD_all_cl3")

sum(results$cor_cl1 == results$JSD_markers_cl1) 
sum(results$cor_cl2 == results$JSD_markers_cl2) 
sum(results$cor_cl3 == results$JSD_markers_cl3) 

sum(results$cor_cl1 == results$JSD_all_cl1) 
sum(results$cor_cl2 == results$JSD_all_cl2) 
sum(results$cor_cl3 == results$JSD_all_cl3) 

temp <- matrix(matrix_cor, nrow = 404, ncol = 404) 
tt1 <- as.vector(temp)
tt2 <- sapply(tt1, function(x)match(x, sort(tt1, decreasing = TRUE)))
ranked_cor_mat <- matrix(tt2, nrow = 404, ncol = 404)

temp <- matrix(matrix_JSD, nrow = 404, ncol = 404) 
tt1 <- as.vector(temp)
tt2 <- sapply(tt1, function(x)match(x, sort(tt1)))
ranked_JSD_mat <- matrix(tt2, nrow = 404, ncol = 404)




##########################################################################################
### Adding markers to the dend: ##########################################################
##########################################################################################
#We prune the dend based on 10X data
missing_cl <- names(table(cl))[!names(table(cl)) %in% names(table(cl_10X_only))] 
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
#save(dend_pruned, file =paste0(work.dir, "dend_pruned_markers_attached.rda"))
load(paste0(work.dir, "/dend_pruned_markers_attached.rda"))
#save(dend, file="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2/Fahimeh_dend_markers_attached.rda")

##########################################################################################
### Initialization: ######################################################################
##########################################################################################

select.cl <- labels(dend_pruned)
TenX.cells <- names(cl_10X_only[cl_10X_only %in% labels(dend_pruned)]) 

cl_10X_only <- cl_10X_only[TenX.cells]
cl_10X_only <- droplevels(cl_10X_only)

Smartseq.cells <- names(cl_smartseq_only[cl_smartseq_only %in% labels(dend_pruned)])
cl_smartseq_only <- cl_smartseq_only[Smartseq.cells]
cl_smartseq_only <- droplevels(cl_smartseq_only)

#write.csv(as.matrix(TenX.norm.dat[select.markers,]), file = paste0(work.dir, "/TenX_norm_dat.csv"))
#source("//allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/patchseq_10X.R")

##########################################################################################
### Mapping Smartseq using 10X as ref: ###################################################
##########################################################################################

memb_smart2tenx = map_dend_membership(dend_pruned, cl=cl_10X_only, 
                                      TenX.norm.dat[select.markers,], 
                                      Smartseq.norm.dat[select.markers,], 
                                      Smartseq.cells, bs.num=100, 
                                      p=0.7, low.th=0.15, 
                                      cl.mean = cl.means_10X[select.markers, select.cl])

#save(memb_smart2tenx, file = paste0(work.dir, "memb_smart2tenx.rda"))
#No cell was mapped to this cluster:
labels(dend_pruned)[!labels(dend_pruned) %in% colnames(memb_smart2tenx)]
cl.df[cl.df$cl == labels(dend_pruned)[!labels(dend_pruned) %in% colnames(memb_smart2tenx)], "cluster_label"]
mapped.cl <- colnames(memb_smart2tenx)[colnames(memb_smart2tenx) %in% labels(dend_pruned)]

rownames(anno.df) <- anno.df$sample_name
anno.df <- as.data.frame(anno.df)
sorted_cl <- t(apply(memb_smart2tenx[, mapped.cl], 1, function(x) names(x)[order(x, decreasing = TRUE)]))
sorted_bt <- t(apply(memb_smart2tenx[, mapped.cl], 1, function(x) x[order(x, decreasing = TRUE)]))

first_cl <- sorted_cl[,1]
second_cl <- sorted_cl[,2]
first_cluster_label <- Renew_list(first_cl, ref.df = cl.df, label = "cl", new.label = "cluster_label")
second_cluster_label <- Renew_list(second_cl, ref.df = cl.df, label = "cl", new.label = "cluster_label")

sum(as.character(anno.df[Smartseq.cells, "cluster_label"]) == as.character(first_cluster_label[Smartseq.cells])) / length(Smartseq.cells)
sum(as.character(anno.df[Smartseq.cells, "cl"]) == as.character(first_cl[Smartseq.cells])) / length(Smartseq.cells)
sum(as.character(anno.df[Smartseq.cells, "cluster_label"]) == as.character(second_cluster_label[Smartseq.cells])) / length(Smartseq.cells)
sort(table(first_cluster_label))

##########################################################################################
### Mapping 10X on 10X: ##################################################################
##########################################################################################
TenX_test_cells <- sample(TenX.cells, 35000)
TenX_train_cells <- setdiff(TenX.cells, TenX_test_cells)

memb_tenx2tenx = map_dend_membership(dend_pruned, cl=cl_10X_only[TenX_train_cells], 
                                      TenX.norm.dat[select.markers, TenX_train_cells], 
                                      TenX.norm.dat[select.markers, TenX_test_cells], 
                                      TenX_test_cells, bs.num=100, 
                                      p=0.7, low.th=0.15, 
                                      cl.mean = cl.means_10X[select.markers, select.cl])

save(memb_tenx2tenx, file = paste0(work.dir, "/memb_tenx2tenx.rda"))

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
