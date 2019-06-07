.libPaths("/home/fahimehb/R/x86_64-redhat-linux-gnu-library/3.5")
FIT_ZINB =  FALSE 
COMPUTE_BETA = TRUE 
COMPUTE_LIKELIHOOD = TRUE 
AGGREGATE_RESULTS = TRUE
args = commandArgs(trailingOnly=TRUE)
run_iter = as.numeric(args[1]) 

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
markers.path <- paste0(work.dir, "/1000_markers.rda")

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

load(markers.path)
print("Done loading markers list!!")
#Long_markers_list <- select.markers
#Long_markers_list <- sample(select.markers, 500)
#Long_markers_list[grepl( "Rik" , Long_markers_list) ]  <- paste0("rename",Long_markers_list[grepl( "Rik" , Long_markers_list) ])
#Long_markers_list <- gsub("-", "_", Long_markers_list)
#rm(select.markers)
short_markers_list <- c("Npy", "Npy1r", "Npy2r", "Npy5r", 
                        "Sst" ,"Sstr1", "Sstr2", "Sstr3", "Sstr4", "Cort",
                        "Vip", "Vipr1", "Vipr2",
                        "Tac2", "Tacr3", 
                        "Cck", "Cckbr", 
                        "Penk", "Oprd1", "Oprm1", 
                        "Crh", "Crhr1", "Crhr2",
                        "Tac1", "Tacr1", 
                        "Pdyn", "Oprk1", 
                        "Pthlh","Pth1r",
                        "Pnoc", "Oprl1",
                        "Trh", "Trhr", "Trhr2", 
                        "Grp", "Grpr",
                        "Rln1", "Rxfp1", "Rxfp2", "Rxfp3", 
                        "Adcyap1", "Adcyap1r1", 
                        "Nts", "Ntsr1", "Ntsr2", 
                        "Nmb", "Nmbr")


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

#save(FACS.anno,file = paste0(work.dir, "/Counts/FACS_anno.rda"))
#save(FACS.counts,file = paste0(work.dir, "/Counts/FACS_counts.rda"))
#print("Done saving FACS.anno!!")
load(paste0(work.dir, "/Counts/FACS_anno.rda"))
load(paste0(work.dir, "/Counts/FACS_counts.rda"))
select.cl <- setdiff(unique(FACS.anno$cluster_label), LowQ_types)
FACS.anno <- FACS.anno[FACS.anno$cluster_label %in% select.cl,]
rownames(FACS.anno) <- FACS.anno$sample_name
FACS.cells <- rownames(FACS.anno)
#FACS.counts <- FACS.counts[FACS.cells,short_markers_list]
colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ]  <- paste0("rename",colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ])
colnames(FACS.counts) <- gsub("-", "_", colnames(FACS.counts))
FACS.counts <- FACS.counts[FACS.cells,Long_markers_list]
df <- cbind(FACS.anno$cluster_id, as.data.frame.matrix(FACS.counts))
colnames(df) <- c("Type", colnames(FACS.counts))
print(dim(df))
##########################################################################################
### Some initialization: #################################################################
##########################################################################################

temp <- table(FACS.anno$cluster_id) >= 10 
Good_types <- as.numeric(names(temp[temp]))
Good_pairs <- t(combn(Good_types, 2))
#Adding pure types as 1_1, 2_2 and ...
for (t in Good_types) {
  Good_pairs <- rbind(Good_pairs, c(t, t))
}
#genes <- short_markers_list
genes <- Long_markers_list
print("Done initialization!!!")
##########################################################################################
### fit ZINB or NB or logit per gene per Type: ###########################################
##########################################################################################

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/ZINB_helper_functions.R")
if (FIT_ZINB){
  All_fit_values <- list()
  for (t in Good_types) {
    print(c("Type:", t))
    new_df <- df[df$Type == t ,]
    Fit_values <- tapply(1:nrow(new_df), as.character(new_df$Type), function(x) sapply(colnames(new_df)[-1], function(g)Fit_model(new_df[x,], g, t)))
    All_fit_values <- c(All_fit_values, Fit_values)
  }
  save(All_fit_values, file=paste0(work.dir,"/All_fit_values.rda"))
} else {
load(paste0(work.dir,"/1000_fit_values.rda")) 
} 
df <- t(df)

##########################################################################################
### Setting RUN_ITER: ####################################################################
##########################################################################################

#run_iter <- 3000

##########################################################################################
### Find Beta for each pair of cell type using an optimizer: #############################
##########################################################################################

Err <- function(data, par) {
  with(data, sum((par * term1 +  term2)^2))
}

Beta_file_name = paste0(FACS.cells[run_iter],"_Beta.rda")
if (COMPUTE_BETA){
  Beta_list <- list()
  for (c in FACS.cells[run_iter]){
    B <- list()
    p1 <- list()
    p2 <- list()
    start_time = Sys.time()
    for (i in 1:dim(Good_pairs)[1]) {
      first_type <- as.character(Good_pairs[i,1])
      second_type <- as.character(Good_pairs[i,2])
      if (first_type == second_type) {
        B <- c(B, 1)
      } else {
        pair_name <- paste0(as.character(first_type), "_", as.character(second_type))
        list1 <- as.data.frame(t(All_fit_values[[first_type]]))
        list2 <- as.data.frame(t(All_fit_values[[second_type]]))
        y = df[genes, c]
        if(length(y) == 0) {print("ERROR: Something is wrong!!!!!!!")}
        Mu1 <- unlist(list1[genes, "Mu"])
        Mu2 <- unlist(list2[genes, "Mu"])
        Pi1 <- unlist(list1[genes, "Pi"])
        Pi2 <- unlist(list2[genes, "Pi"])
        term1 <- Mu1 - Mu1 * Pi1 - Mu2 + Mu2 * Pi2
        term2 <- Mu2 - Mu2 * Pi2 - y
        dat <- as.data.frame(cbind(term1, term2))
        results <- optim(par = runif(1, 0.09, 0.9),
                       fn = Err,
                       data = dat,
                       method = "Brent",
                       lower = 0.01,
                       upper = 0.99)
        B <- c(B, results$par)
      }
      p1 <- c(p1, first_type)
      p2 <- c(p2, second_type)
    }
    temp <- cbind.data.frame(unlist(B), unlist(p1), unlist(p2))
    colnames(temp) <- c("B", "Type1", "Type2")
    Beta_list[[c]] <- temp
    end_time = Sys.time()
    print(end_time - start_time)
  }
  print("Done!")
  save(Beta_list, file=paste0(work.dir, "/Beta_files/",Beta_file_name))
} else {load(paste0(work.dir, "/Beta_files/",Beta_file_name))}


##########################################################################################
### Find the probability of the observation given Beta: ##################################
##########################################################################################

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/ZINB_helper_functions.R")

pyg_file_name = paste0(FACS.cells[run_iter],"_pyg.rda")
count_threshold <- max(df)
if (COMPUTE_LIKELIHOOD){
  pyg <- list()
  for (c in FACS.cells[run_iter]){
    y <- df[genes,c]
    start_time = Sys.time()
    pp <- list()
    l <- list()
    for(j in Good_types){
      print(j)
      list <- as.data.frame(t(All_fit_values[[as.character(j)]]))
      pair_name <- paste0(as.character(j), "_", as.character(j))
      pyg[[pair_name]] <- lapply(names(y), function(x){Compute_pure_pyg(x, y[[x]], list, c, as.character(j))})
      pyg[[pair_name]][pyg[[pair_name]] == 0] <- 0.000001
      x <- sum(log(unlist(pyg[[pair_name]])))
      if (!is.infinite(x)) {
        l <- c(l,x)
        pp <- c(pp, j)
      }
    }
    Type_subset <- unlist(pp)[order(unlist(l), decreasing = TRUE)][1:5]
    Good_pairs_subset <- expand.grid(Type_subset, Good_types)
    for (i in 1:dim(Good_pairs_subset)[1]) {
      first_type <- as.character(Good_pairs_subset[i,1])
      second_type <- as.character(Good_pairs_subset[i,2])
      if (first_type != second_type){
        pair_name <- paste0(as.character(first_type), "_", as.character(second_type))
        print(pair_name)
        list1 <- as.data.frame(t(All_fit_values[[first_type]]))
        list2 <- as.data.frame(t(All_fit_values[[second_type]]))
        pyg[[pair_name]] <- lapply(names(y), function(x){Compute_pyg(x, y[[x]], list1, list2, c, first_type, second_type, count_threshold)})
      }
    }
    end_time = Sys.time()
    print(end_time - start_time)
  }
  save(pyg, file=paste0(work.dir, "/Probability_files/",pyg_file_name))
} else {load(paste0(work.dir, "/Probability_files/",pyg_file_name))}


##########################################################################################
### Find the probability of the observation given Beta: ##################################
##########################################################################################

cluster_lable_id <- as.data.frame(unique(FACS.anno[,c("cluster_label", "cluster_id")]))
cell <- FACS.cells[run_iter]

if (AGGREGATE_RESULTS){
  Finaldf <- list()
  py <- list()
  for (c in cell){
    Type1 <- list()
    Type2 <- list()
    cl1 <- list()
    cl2 <- list()
    logp <- list()
    for (i in 1:dim(Good_pairs_subset)[1]) {
      first_type <- as.character(Good_pairs_subset[i,1])
      second_type <- as.character(Good_pairs_subset[i,2])
      pair_name <- paste0(as.character(first_type), "_", as.character(second_type))
      pyg[[pair_name]][pyg[[pair_name]] == 0] <- 0.000001
      logp <- c(logp , sum(log(unlist(pyg[[pair_name]]))))
      Type1 <- c(Type1, first_type)
      Type2 <- c(Type2, second_type)
      cl1 <- c(cl1, cluster_lable_id[cluster_lable_id$cluster_id == first_type, "cluster_label"])
      cl2 <- c(cl2, cluster_lable_id[cluster_lable_id$cluster_id == second_type, "cluster_label"])
    }
    py[[c]] <- cbind.data.frame(unlist(logp), unlist(Type1), unlist(Type2), unlist(cl1), unlist(cl2))
    colnames(py[[c]]) <- c("logp", "Type1", "Type2", "cl1", "cl2")
    Finaldf[[c]] <- merge(py[[c]], Beta_list[[c]], by=c("Type1","Type2"))
    Finaldf[[c]][["cluster_id"]] <- FACS.anno[[cell, c( "cluster_id")]]
    Finaldf[[c]][["cluster_label"]] <- FACS.anno[[cell, c( "cluster_label")]]
    Finaldf[[c]]["pair_identity"] <- ifelse(Finaldf[[c]][,"Type1"] == Finaldf[[c]][,"cluster_id"] |  Finaldf[[c]][,"Type2"] == Finaldf[[c]][,"cluster_id"] , "U_or_V", "None") 
  }
  aggregated_results_file_name = paste0(FACS.cells[run_iter],"_results.rda")
  save(Finaldf, file=paste0(work.dir, "/ZINB_results_files/",aggregated_results_file_name))
}


#select.cell <- FACS.cells[run_iter]
#ggplot(Finaldf[[select.cell]][!is.na(Finaldf[[select.cell]][,"logp"]),], aes(logp, B, colour = pair_identity)) + 
#  geom_point(alpha = 0.4) + xlab("Log(p)") + ylab("B")

#Finaldf[[select.cell]][order(Finaldf[[select.cell]][,"logp"], decreasing = TRUE),][1:10,]
