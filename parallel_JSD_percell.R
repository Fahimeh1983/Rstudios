.libPaths("/home/fahimehb/R/x86_64-redhat-linux-gnu-library/3.5")
library(matrixStats)
library(feather)
library(Matrix)
library("philentropy")
library(emdist)

##########################################################################################
### Reading inputs from command line: ####################################################
##########################################################################################

args = commandArgs(trailingOnly=TRUE)
n_iter = as.numeric(args[1])
n_loop = as.numeric(args[2])

##########################################################################################
### Reading inputs: ######################################################################
##########################################################################################
work.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/"
ref.data.rda.path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"
facs.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"

load(paste0(ref.data.rda.path, "ref.data.rda"))
FACS.anno <- read_feather(paste0(facs.dir, "/anno.feather"))
FACS.anno <- as.data.frame(FACS.anno)
rownames(FACS.anno) <- FACS.anno$sample_id
FACS.anno <- FACS.anno[names(cl),]
Sst.cells <- rownames(FACS.anno)[(FACS.anno$subclass_label == "Sst" & FACS.anno$region_label == "VISp")]
Sst.cl <- setNames(FACS.anno[Sst.cells, "cluster_id"], Sst.cells)


#Sst.norm.dat <- norm.dat[, Sst.cells]
#rm(norm.dat)
#t_prob <- (2 ^ Sst.norm.dat) - 1
#t_prob <- t(t_prob)
#pairs <- t(combn(Sst.cells,2))

#save(t_prob, file = paste0(work.dir, "/t_prob.rda"))
#save(pairs, file = paste0(work.dir, "/pairs.rda"))
#load("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/select.markers.rda")
#Sst.norm.dat <- Sst.norm.dat[select.markers,]
#save(Sst.norm.dat, file=paste0(work.dir, "Sst.norm.dat.rda"))

##########################################################################################
### Compute JSD per cell on hpc: #########################################################
##########################################################################################
load (paste0(work.dir, "/t_prob.rda"))
print("Loaded t_prob")
load(paste0(work.dir, "/pairs.rda"))
print("Loaded pairs")
t_prob <- t(t_prob)
t_prob[t_prob==0] <- 0.000001
t_prob <- t_prob/colSums(t_prob)
t_prob <- t(t_prob)

file_name <- paste0("JSD_results_", n_iter)
begin = (n_iter-1) * n_loop + 1
end = begin + n_loop -1
if (end > dim(pairs)[1]) {end <- dim(pairs)[1]}
results_JSD <- matrix(nrow = n_loop, ncol = 3)

for (i in begin:end){
  select.pair <- pairs[i,]
  t_prob_pair <- t_prob[select.pair, ]
  pair_JSD <- JSD(t_prob_pair,  test.na = TRUE, unit = "log2", est.prob = NULL)
  results_JSD[i-begin+1, 1] <- select.pair[1]
  results_JSD[i-begin+1, 2] <- select.pair[2]
  results_JSD[i-begin+1, 3] <- pair_JSD
}

save(results_JSD, file= paste0(work.dir, file_name, ".rda"))
#load(paste0(work.dir, "JSD_results_1.rda"))

##########################################################################################
### Compute EMD per cell on hpc: #########################################################
##########################################################################################
# load (paste0(work.dir, "/t_prob.rda"))
# print("Loaded t_prob")
# load(paste0(work.dir, "/pairs.rda"))
# print("Loaded pairs")
# 
# file_name <- paste0("EMD_results_", n_iter)
# begin = (n_iter-1) * n_loop + 1
# end = begin + n_loop -1
# if (end > dim(pairs)[1]) {end <- dim(pairs)[1]}
# results_EMD <- matrix(nrow = n_loop, ncol = 3)
# 
# 
# for (i in begin:end){
#     select.pair <- pairs[i,]
#     t_prob_pair1 <- t_prob[select.pair[1], ]
#     t_prob_pair2 <- t_prob[select.pair[2], ]
#     pair_EMD <- emd(t_prob_pair1 , t_prob_pair2)
# 
#     results_EMD[i-begin+1, 1] <- select.pair[1]
#     results_EMD[i-begin+1, 2] <- select.pair[2]
#     results_EMD[i-begin+1, 3] <- pair_EMD
# }
# 
# save(results_EMD, file= paste0(work.dir, file_name, ".rda"))
#load(paste0(work.dir, "JSD_results_1.rda"))

#########################################################################################
### Compute cor per cell on hpc: #########################################################
##########################################################################################
# load (paste0(work.dir, "/Sst.norm.dat.rda"))
# print("Loaded Sst.norm.dat.rda")
# load(paste0(work.dir, "/pairs.rda"))
# print("Loaded pairs")
# Sst.norm.dat <- t(Sst.norm.dat)
# 
# file_name <- paste0("cor_results_", n_iter)
# begin = (n_iter-1) * n_loop + 1
# end = begin + n_loop -1
# if (end > dim(pairs)[1]) {end <- dim(pairs)[1]}
# results_cor <- matrix(nrow = n_loop, ncol = 3)
# 
# for (i in begin:end){
#   select.pair <- pairs[i,]
#   Sst.norm.dat.pair1 <- Sst.norm.dat[select.pair[1], ]
#   Sst.norm.dat.pair2 <- Sst.norm.dat[select.pair[2], ]
#   pair_cor <- cor(Sst.norm.dat.pair1, Sst.norm.dat.pair2)
#   results_cor[i-begin+1, 1] <- select.pair[1]
#   results_cor[i-begin+1, 2] <- select.pair[2]
#   results_cor[i-begin+1, 3] <- pair_cor
# }
# 
# save(results_cor, file= paste0(work.dir, file_name, ".rda"))
#load(paste0(work.dir, "JSD_results_1.rda"))

##########################################################################################
### Reading JSD output: ##################################################################
##########################################################################################
get_pair_matrix_coor <- function(m, rows, cols)
{
  if(!is.numeric(rows)){
    rows <- match(rows, row.names(m))
  }

  if(!is.numeric(cols)){
    cols <- match(cols, colnames(m))
  }

  coor <- (cols - 1) * nrow(m) + rows

  return(coor)
}

set_pair_matrix <- function(m, rows, cols, vals)
{
  coor <- get_pair_matrix_coor(m, rows, cols)
  m[coor] <- vals
  return(m)
}


load(paste0(work.dir, "/pairs.rda"))
print("Loaded pairs")
JSD_all <- matrix(nrow = 16*100000, ncol = 3)

for (i in 1:16){
  file_name = paste0("/JSD_results_", i, ".rda")
  load(paste0(work.dir, file_name))
  begin = (i-1) * dim(results_JSD)[1] + 1
  end = begin + dim(results_JSD)[1] -1
  JSD_all[begin:end, ] <- results_JSD
}

JSD_all <- JSD_all[1:dim(pairs)[1],]
JSD_all <- as.data.frame(JSD_all)
colnames(JSD_all) <- c("p1", "p2", "dist")
JSD_all$p1 <- as.character(JSD_all$p1)
JSD_all$p2 <- as.character(JSD_all$p2)
JSD_all$dist <- as.numeric(JSD_all$dist)

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")

JSD_mat <- matrix(0, ncol=length(Sst.cells), nrow=length(Sst.cells))
colnames(JSD_mat) <- Sst.cells
rownames(JSD_mat) <- Sst.cells

JSD_mat = set_pair_matrix(JSD_mat, JSD_all$p1, JSD_all$p2, JSD_all$dist)
tmp <- t(JSD_mat)
JSD_mat[lower.tri(JSD_mat, diag = FALSE)] <- tmp[lower.tri(tmp, diag = FALSE)]

JSD_mat <- JSD_mat/max(JSD_mat)

col <- t(as.data.frame(setNames(FACS.anno[Sst.cells, "cluster_color"], Sst.cells)))
plot_co_matrix(co.ratio = (1-JSD_mat), cl = Sst.cl, max.cl.size = 50)
save(JSD_mat, file = paste0(work.dir, "/jsd_markers_genes/JSD_makers_mat.rda"))
#dev.off()



##########################################################################################
### Reading cor output: ##################################################################
##########################################################################################
# load(paste0(work.dir, "/pairs.rda"))
# print("Loaded pairs")
# cor_all <- matrix(nrow = 16*100000, ncol = 3)
# 
# for (i in 1:16){
#   file_name = paste0("cor_markers_genes/cor_results_", i, ".rda")
#   load(paste0(work.dir, file_name))
#   begin = (i-1) * dim(results_cor)[1] + 1
#   end = begin + dim(results_cor)[1] -1
#   cor_all[begin:end, ] <- results_cor
# }
# 
# cor_all <- cor_all[1:dim(pairs)[1],]
# cor_all <- as.data.frame(cor_all)
# colnames(cor_all) <- c("p1", "p2", "dist")
# cor_all$p1 <- as.character(cor_all$p1)
# cor_all$p2 <- as.character(cor_all$p2)
# cor_all$dist <- as.numeric(as.character(cor_all$dist))
# 
# cor_mat <- matrix(1, ncol=length(Sst.cells), nrow=length(Sst.cells))
# colnames(cor_mat) <- Sst.cells
# rownames(cor_mat) <- Sst.cells
# 
# cor_mat = set_pair_matrix(cor_mat, cor_all$p1, cor_all$p2, cor_all$dist)
# tmp <- t(cor_mat)
# cor_mat[lower.tri(cor_mat, diag = FALSE)] <- tmp[lower.tri(tmp, diag = FALSE)]
# 
# col <- t(as.data.frame(setNames(FACS.anno[Sst.cells, "cluster_color"], Sst.cells)))
# load("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/V1_ALM/process_24411/all.col.rda")
# plot_co_matrix(co.ratio = cor_mat, cl = Sst.cl, max.cl.size = 25, col= col2)
# plot_co_matrix(co.ratio = cor_mat, cl = Sst.cl, max.cl.size = 25, col= all.col[9,names(Sst.cl),drop=F])
# 
# save(cor_mat, file = paste0(work.dir, "/cor_markers_genes/cor_makers_mat.rda"))
# 
##########################################################################################
### Reading cor output: ##################################################################
##########################################################################################
pairs <- facs.pairs
subsamp_FACS_cl <- cl[subsamp_FACS_cells]
EMD_all <- matrix(nrow = 44*100000, ncol = 3)

for (i in 1:44){
  print(i)
  file_name = paste0("/JSD/EMD_results_", i, ".csv")
  results_EMD <- read.csv(paste0(work.dir, file_name))
  begin = (i-1) * dim(results_EMD)[1] + 1
  end = begin + dim(results_EMD)[1] -1
  EMD_all[begin:end, 1] <- results_EMD$p1
  EMD_all[begin:end, 2] <- results_EMD$p2
  EMD_all[begin:end, 3] <- results_EMD$dist
}

EMD_all <- EMD_all[1:4400000,]
EMD_all <- as.data.frame(EMD_all)
colnames(EMD_all) <- c("p1", "p2", "dist")
EMD_all$p1 <- as.character(EMD_all$p1)
EMD_all$p2 <- as.character(EMD_all$p2)
EMD_all$dist <- as.numeric(as.character(EMD_all$dist))

EMD_mat <- matrix(1, ncol=length(subsamp_FACS_cells), nrow=length(subsamp_FACS_cells))
colnames(EMD_mat) <- subsamp_FACS_cells
rownames(EMD_mat) <- subsamp_FACS_cells

EMD_mat = set_pair_matrix(EMD_mat, EMD_all$p1, EMD_all$p2, EMD_all$dist)
tmp <- t(EMD_mat)
EMD_mat[lower.tri(EMD_mat, diag = FALSE)] <- tmp[lower.tri(tmp, diag = FALSE)]

library(reshape2)
ggplot(data = melt(EMD_mat), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()+ theme(axis.text = element_text(size=7)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Clustering lable") + ylab("NN Mapping lables") +
  scale_fill_gradient(low = "white", high = "red")

#freqs <- scale(EMD_mat, center = FALSE, scale = colSums(EMD_mat))
#class(EMD_mat)
#EMD_mat <- EMD_mat/max(EMD_mat)
col <- t(as.data.frame(setNames(FACS.anno[subsamp_FACS_cells, "cluster_color"], subsamp_FACS_cells)))
#load("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/SmartSeq_cells/V1_ALM/process_24411/all.col.rda")
plot_co_matrix(co.ratio = 1-EMD_mat, cl = subsamp_FACS_cl, max.cl.size = 5, col=col)
plot_co_matrix(co.ratio = EMD_mat, cl = subsamp_FACS_cl , max.cl.size = 25, col= all.col[9,names(Sst.cl),drop=F])

# save(EMD_mat, file = paste0(work.dir, "/emd/cor_makers_mat.rda"))


##########################################################################################
### JSD and cor for cl.med: ##############################################################
##########################################################################################

# Sst.cl.norm.med <- do.call("cbind", tapply(FACS.anno[Sst.cells, "sample_id"], FACS.anno[Sst.cells, "cluster_label"], function(x)rowMedians(norm.dat[select.markers, x, drop=F])))
# t_prob <- t(t_prob)
# Sst.cl.t_prob.med <- do.call("cbind", tapply(FACS.anno[Sst.cells, "sample_id"], FACS.anno[Sst.cells, "cluster_label"], function(x)rowMedians(t_prob[select.markers, x, drop=F])))
# 
# Sst.cl.t_prob.med <- t(Sst.cl.t_prob.med)
# 
# JSD.cl.med <- JSD(Sst.cl.t_prob.med,  test.na = TRUE, unit = "log2", est.prob = NULL)
# cor.cl.med <- cor(Sst.cl.norm.med, Sst.cl.norm.med)
# JSD.cl.med <- JSD.cl.med/max(JSD.cl.med)
# rownames(JSD.cl.med) <- rownames(cor.cl.med)
# colnames(JSD.cl.med) <- colnames(cor.cl.med)
# ord_lab <- FACS.anno %>% select(cluster_label, dendcluster_id) %>% unique() %>% arrange(dendcluster_id) %>% select(cluster_label)
# ord_lab <- ord_lab$cluster_label
# ord_lab <- ord_lab[ord_lab %in% rownames(JSD.cl.med)]
# cor.cl.med <- cor.cl.med[ord_lab, ord_lab]
# JSD.cl.med <- JSD.cl.med[ord_lab, ord_lab]
# rotate <- function(x) t(apply(x, 2, rev))
# rotated.cor.cl.med <- rotate(cor.cl.med)
# rotated.JSD.cl.med <- rotate(JSD.cl.med)
# 
# library(reshape2)
# ggplot(data = melt(rotated.cor.cl.med), aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile()+ theme(axis.text = element_text(size=7)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#   xlab("Clustering lable") + ylab("NN Mapping lables") +  
#   scale_fill_gradient(low = "white", high = "red")
# 
# library(reshape2)
# ggplot(data = melt(1-rotated.JSD.cl.med), aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile()+ theme(axis.text = element_text(size=7)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#   xlab("Clustering lable") + ylab("NN Mapping lables") +  
#   scale_fill_gradient(low = "white", high = "red")

load(paste0(work.dir, "/facs_memb.rda"))
write.csv(FACS.memb, file = paste0(work.dir, "/facs_memb.csv"))
FACS.cells <- rownames(FACS.memb)
facs.pairs <- t(combn(FACS.cells,2))
write.csv(facs.pairs, file = paste0(work.dir, "/facs_pairs.csv"))
subsamp_FACS_cells <- sample(FACS.cells, 3000)
facs.pairs <- t(combn(subsamp_FACS_cells,2))
write.csv(facs.pairs, file = paste0(work.dir, "/facs_pairs.csv"))
length(table(cl[subsamp_FACS_cells]))

