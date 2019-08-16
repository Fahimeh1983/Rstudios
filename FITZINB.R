##########################################################################################
### Reading the counts: ##################################################################
##########################################################################################
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")
#devtools::install_github("AllenInstitute/scrattch.io", ref = "dev")
require(ggplot2)
require(pscl)
require(MASS)
require(boot)
library(dplyr)
library(tidyr)
library(feather)
library(scrattch.io)

work.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/Patchseq_vs_FACs_cre_analysis/mouse_patchseq_VISp_20181220_collapsed40_cpm/"
tome <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/tomes/facs/mouse_V1_ALM_20180520/transcrip.tome"
#markers.path <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/select.markers.rda"
markers.path <- paste0(work.dir, "/1000_marker_genes_include47.rda")

# Sample annotations
# Will throw Warnings, but these are OK - just because of how NAs are stored in HDF5 files.
#FACS.anno <- read_tome_anno(tome)
# Read all counts as sparse matrix
# These are stored in samples (rows) x genes (columns) format
#exon_counts <- read_tome_dgCMatrix(tome = tome, target = "/data/exon")
#intron_counts <- read_tome_dgCMatrix(tome = tome, target = "/data/intron")
#FACS.counts <- exon_counts + intron_counts
#colnames(FACS.counts) <- read_tome_gene_names(tome)
#rownames(FACS.counts) <- read_tome_sample_names(tome)
# See everything stored in the tome
#h5ls(tome)
#print("Done!")
#save(FACS.anno,file = paste0(work.dir, "/Counts/FACS_anno.rda"))
#save(FACS.counts,file = paste0(work.dir, "/Counts/FACS_counts.rda"))
load(paste0(work.dir, "/Counts/FACS_anno.rda"))
load(paste0(work.dir, "/Counts/FACS_counts.rda"))

load(markers.path)
#Long_markers_list <- select.markers
#Long_markers_list <- sample(select.markers, 1000)
#short_markers_list <- c("Npy", "Npy1r", "Npy2r", "Npy5r",
#                        "Sst" ,"Sstr1", "Sstr2", "Sstr3", "Sstr4", "Cort",
#                        "Vip", "Vipr1", "Vipr2",
#                        "Tac2", "Tacr3",
#                        "Cck", "Cckbr",
#                        "Penk", "Oprd1", "Oprm1",
#                        "Crh", "Crhr1", "Crhr2",
#                        "Tac1", "Tacr1",
#                        "Pdyn", "Oprk1",
#                        "Pthlh","Pth1r",
#                        "Pnoc", "Oprl1",
#                        "Trh", "Trhr", "Trhr2",
#                        "Grp", "Grpr",
#                        "Rln1", "Rxfp1", "Rxfp2", "Rxfp3",
#                        "Adcyap1", "Adcyap1r1",
#                        "Nts", "Ntsr1", "Ntsr2",
#                        "Nmb", "Nmbr")
#Long_markers_list <- sample(setdiff(Long_markers_list, short_markers_list), 953)
#Long_markers_list <- union(Long_markers_list, short_markers_list)
#Long_markers_list[grepl( "Rik" , Long_markers_list) ]  <- paste0("rename",Long_markers_list[grepl( "Rik" , Long_markers_list) ])
#Long_markers_list <- gsub("-", "_", Long_markers_list)
#rm(select.markers)

#Removing all low quality cells from FACS data
LowQ_types <- c("Low Quality VISp L5 PT Ctxn3 2", "Batch Grouping VISp L5 PT Chrna6",
                "Batch Grouping VISp L5 PT Ctxn3", "Low Quality VISp L6 CT Ptprt_2",
                "Low Quality VISp L5 PT Ctxn3 1", "Doublet SMC and Glutamatergic",
                "Doublet Astro Aqp4 Ex", "Low Quality ALM L6 CT Cpa6",
                "Low Quality Meis2 Adamts19 ", "Doublet Endo and Peri_1",
                "Doublet VISp L5 NP and L6 CT", "Low Quality Sst Chodl",
                "Low Quality Astro Aqp4" , "Doublet Endo Peri SMC",
                "Low Quality L4 Rspo1", "High Intronic VISp L5 Endou",
                "Low Quality VISp L6 CT Ptprt_1")

highQ_type <- setdiff(unique(FACS.anno$cluster_label), LowQ_types)
FACS.anno <- as.data.frame(FACS.anno) 
rownames(FACS.anno) <- FACS.anno$sample_name
#test.BAD.cells <- FACS.anno[FACS.anno$cluster_label %in% LowQ_types, "sample_name"]
#all.GOOD.cells <- FACS.anno[!FACS.anno$cluster_label %in% LowQ_types, "sample_name"]
#test.GOOD.cells <- sample(all.GOOD.cells, 2000)
#validation.cells <- sample(setdiff(all.GOOD.cells, test.GOOD.cells), 2500)
#tmp<- table(FACS.anno[validation.cells, "cluster_label"])
#sort(tmp, decreasing = TRUE)
#length(tmp)
#train.cells <- setdiff(rownames(FACS.anno), c(test.BAD.cells, test.GOOD.cells, validation.cells))
#save(test.GOOD.cells, file = paste0(work.dir, "/Test_GOOD_cells.rda"))
#save(test.BAD.cells, file = paste0(work.dir, "/Test_BAD_cells.rda"))
#save(validation.cells, file = paste0(work.dir, "/validation_cells.rda"))
#save(train.cells, file = paste0(work.dir, "/train_cells.rda"))
load(paste0(work.dir, "/Test_GOOD_cells.rda"))
load(paste0(work.dir, "/Test_BAD_cells.rda"))
load(paste0(work.dir, "/validation_cells.rda"))
load(paste0(work.dir, "/train_cells.rda"))
sum(validation.cells %in% c(test.GOOD.cells, test.BAD.cells, train.cells)) == 0
sum(test.BAD.cells %in% c(test.GOOD.cells, validation.cells, train.cells)) == 0
sum(test.GOOD.cells %in% c(validation.cells, test.BAD.cells, train.cells)) == 0
sum(train.cells %in% c(test.GOOD.cells, test.BAD.cells, validation.cells)) == 0
dim(FACS.anno)
#cluster_lable_id <- as.data.frame(unique(FACS.anno[,c("cluster_label", "cluster_id")]))
FACS.anno <- FACS.anno[train.cells,] #We removed the test cells from the anno file
FACS.counts <- FACS.counts[train.cells, ]#We removed the test cells from the anno file
dim(FACS.anno)

temp <- table(FACS.anno[FACS.anno$cluster_label %in% highQ_type,"cluster_id"]) >= 10 
Good_types <- as.numeric(names(temp[temp]))
Good_pairs <- t(combn(Good_types, 2))
#Adding pure types as 1_1, 2_2 and ...
for (t in Good_types) {
  Good_pairs <- rbind(Good_pairs, c(t, t))
}
save(Good_types, file = paste0(work.dir, "/Good_trained_types.rda"))
save(Good_pairs, file = paste0(work.dir, "/Good_trained_pairs.rda"))
colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ]  <- paste0("rename",colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ])
colnames(FACS.counts) <- gsub("-", "_", colnames(FACS.counts))
FACS.counts <- FACS.counts[train.cells, Long_markers_list]
sum(rownames(FACS.anno) == rownames(FACS.counts)) == length(train.cells)
df <- cbind(FACS.anno$cluster_id, as.data.frame.matrix(FACS.counts))
colnames(df) <- c("Type", colnames(FACS.counts))
dim(df)

##########################################################################################
### Temporary: ###########################################################################
##########################################################################################
tmp <- list()
load(paste0(work.dir,"/1_10.rda"))
tmp <- All_fit_values
load(paste0(work.dir, "/11_20.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/21_30.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/31_40.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/41_50.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/51_60.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/61_70.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/71_80.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/81_90.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/91_100.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/101_110.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/111_120.rda"))
tmp <- c(tmp, All_fit_values)
load(paste0(work.dir, "/121_130.rda"))
tmp <- c(tmp, All_fit_values)
All_fit_values <- tmp
save(All_fit_values, file=paste0(work.dir, "/All_fit_values_4020.rda"))
##########################################################################################
### Some initialization: #################################################################
##########################################################################################

#genes <- short_markers_list
genes <- Long_markers_list
#save(Long_markers_list, file = paste0(work.dir,"/1000_markers.rda"))

##########################################################################################
### fit ZINB or NB or logit per gene per Type: ###########################################
##########################################################################################

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/ZINB_helper_functions.R")
if (FIT_ZINB){
  All_fit_values <- list()
  for (t in Good_types[1]) {
    start_time = Sys.time()
    print(c("Type:", t))
    new_df <- df[df$Type == t ,]
    Fit_values <- tapply(1:nrow(new_df), as.character(new_df$Type), function(x) sapply(colnames(new_df)[-1], function(g)Fit_model(new_df[x,], g, t)))
    All_fit_values <- c(All_fit_values, Fit_values)
    end_time = Sys.time()
    print(end_time - start_time)
  }
  save(All_fit_values, file = paste0(work.dir, "/All_fit_values_4020.rda"))
} else {
  load(paste0(work.dir, "/All_fit_values_1000_include47genes.rda"))
}

##########################################################################################
### Compute median gene expression in each type: #########################################
##########################################################################################
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

median_expression <- data_frame()
row_names <- c()
t_type <- c()
for (t in Good_types){
  print(t)
  new_df <- df[df$Type == t ,]
  row_names <- c(row_names, paste0("median_of_type_", t))
  t_type <- c(t_type, t)
  median_expression <- rbind(median_expression, colMedians(as.matrix(new_df[-1])))
}

colnames(median_expression) <- colnames(df[-1])
rownames(median_expression) <- row_names
median_expression <- cbind(median_expression, t_type)

save(median_expression, file = paste0(work.dir, "/median_expression.rda"))
df <- t(median_expression)


