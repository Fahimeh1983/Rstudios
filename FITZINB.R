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
#GOOD.cells <- sample(rownames(FACS.anno)[FACS.anno$cluster_label %in% highQ_type], 1000)
#BAD.cells <- rownames(FACS.anno)[FACS.anno$cluster_label %in% LowQ_types]
#save(GOOD.cells, file=paste0(work.dir, "Test_GOOD_cells.rda"))
#save(BAD.cells, file=paste0(work.dir, "Test_BAD_cells.rda"))
#save(Long_markers_list, file = paste0(work.dir, "1000_marker_genes_include47.rda"))

load(paste0(work.dir, "Test_GOOD_cells.rda"))
load(paste0(work.dir, "Test_BAD_cells.rda"))
rownames(FACS.anno) <- FACS.anno$sample_name
train.cells <- setdiff(rownames(FACS.anno), c(GOOD.cells, BAD.cells))
cluster_lable_id <- as.data.frame(unique(FACS.anno[,c("cluster_label", "cluster_id")]))
FACS.anno <- FACS.anno[train.cells,] #We removed the test cells from the anno file
FACS.counts <- FACS.counts[train.cells, ]#We removed the test cells from the anno file

#if (ANALYSIS_FOR_GOOD_CELLS){
#  print("Analysis is being done for highQ cells!")
#  select.cl <- setdiff(unique(FACS.anno$cluster_label), LowQ_types)
#} else {
#  print("Analysis is being done for lowQ cells!")
#  select.cl <- LowQ_types
#  FIT_ZINB <- FALSE
#}

temp <- table(FACS.anno[FACS.anno$cluster_label %in% highQ_type,"cluster_id"]) >= 10 
Good_types <- as.numeric(names(temp[temp]))
Good_pairs <- t(combn(Good_types, 2))
#Adding pure types as 1_1, 2_2 and ...
for (t in Good_types) {
  Good_pairs <- rbind(Good_pairs, c(t, t))
}

rownames(FACS.anno) <- FACS.anno$sample_name
colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ]  <- paste0("rename",colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ])
colnames(FACS.counts) <- gsub("-", "_", colnames(FACS.counts))
FACS.counts <- FACS.counts[train.cells, Long_markers_list]
sum(rownames(FACS.anno) == rownames(FACS.counts)) == length(train.cells)
df <- cbind(FACS.anno$cluster_id, as.data.frame.matrix(FACS.counts))
colnames(df) <- c("Type", colnames(FACS.counts))
dim(df)


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
  for (t in Good_types) {
    start_time = Sys.time()
    print(c("Type:", t))
    new_df <- df[df$Type == t ,]
    Fit_values <- tapply(1:nrow(new_df), as.character(new_df$Type), function(x) sapply(colnames(new_df)[-1], function(g)Fit_model(new_df[x,], g, t)))
    All_fit_values <- c(All_fit_values, Fit_values)
    end_time = Sys.time()
    print(end_time - start_time)
  }
  save(All_fit_values, file = paste0(work.dir, "/All_fit_values_1000_include47genes.rda"))
} else {
  load(paste0(work.dir, "/All_fit_values_1000_include47genes.rda"))
}


