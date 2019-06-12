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

FACS.anno <- as.data.frame(FACS.anno) 

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



temp <- table(FACS.anno[FACS.anno$cluster_label %in% highQ_type,"cluster_id"]) >= 10 
Good_types <- as.numeric(names(temp[temp]))
Good_pairs <- t(combn(Good_types, 2))
#Adding pure types as 1_1, 2_2 and ...
for (t in Good_types) {
  Good_pairs <- rbind(Good_pairs, c(t, t))
}

load(paste0(work.dir, "Test_GOOD_cells.rda"))
load(paste0(work.dir, "Test_BAD_cells.rda"))
rownames(FACS.anno) <- FACS.anno$sample_name
cluster_lable_id <- as.data.frame(unique(FACS.anno[,c("cluster_label", "cluster_id")]))
test.cells <- union(GOOD.cells, BAD.cells)
FACS.anno <- FACS.anno[test.cells,]
colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ]  <- paste0("rename",colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ])
colnames(FACS.counts) <- gsub("-", "_", colnames(FACS.counts))
FACS.counts <- FACS.counts[test.cells, Long_markers_list]
sum(rownames(FACS.anno) == rownames(FACS.counts)) == length(test.cells)
df <- cbind(FACS.anno$cluster_id, as.data.frame.matrix(FACS.counts))
colnames(df) <- c("Type", colnames(FACS.counts))
dim(df)


##########################################################################################
### Some initialization: #################################################################
##########################################################################################

genes <- Long_markers_list
df <- t(df)
load(paste0(work.dir, "/All_fit_values_1000_include47genes.rda"))
##########################################################################################
### Setting RUN_ITER: ####################################################################
##########################################################################################
#run_iter <- 5000
Beta_list <- list()
Finaldf <- list()

for (run_iter in 1:1){
  ##########################################################################################
  ### Find Beta for each pair of cell type using an optimizer: #############################
  ##########################################################################################
  
  Err <- function(data, par) {
    with(data, sum((par * term1 +  term2)^2))
  }
  
  for (c in test.cells[run_iter]){
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
  
  #Beta_file_name = paste0(FACS.cells[run_iter],"_Beta.rda")
  #save(Beta_list, file=paste0(work.dir, "/Beta_files/",Beta_file_name))
  
  ##########################################################################################
  ### Find the probability of the observation given Beta: ##################################
  ##########################################################################################
  
  source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/git_workspace/Rstudios/Rstudios/ZINB_helper_functions.R")
  
  count_threshold <- max(df)
  pyg <- list()
  for (c in test.cells[run_iter]){
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
  
  #pp[which.max(unlist(l))]
  #FACS.anno[[test.cells[run_iter], c( "cluster_id")]]
  #unlist(pp)[order(unlist(l))]
  ##########################################################################################
  ### Find the probability of the observation given Beta: ##################################
  ##########################################################################################
  
  cell <- test.cells[run_iter]

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
    Finaldf[[c]]["pair_identity"] <- ifelse(Finaldf[[c]][,"Type1"] == Finaldf[[c]][,"cluster_id"] |
                                              Finaldf[[c]][,"Type2"] == Finaldf[[c]][,"cluster_id"] , "U_or_V", "None")
  }
} 

Finaldf[[cell]][order(Finaldf[[cell]][,"logp"], decreasing = TRUE),][1:10,]


#is U and V contamination symmetric?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##########################################################################################
### Find the probability of the observation given Beta: ##################################
##########################################################################################
GOOD <- GOOD.FACS.cells
GOOD.FACS.cells <- GOOD.FACS.cells[1:1000] #I did the analysis for these cells
BAD.FACS.cells # I also did analysis for these cells

output <- list()
for (run_iter in 1:589){
   #if (!run_iter %in% c(3,16, 47, 78, 290, 291, 292, 293, 501, 503, 505, 507, 518, 522, 523, 529, 531, 534, 537, 540, 546, 553, 556, 560, 562, 564, 568, 572,  574, 578, 581, 586)){
  #if (!run_iter %in% c(70, 338, 410, 425, 521, 616, 676, 948)){ #Good cells 250
  #if (!run_iter %in% c(948)){ #Good cells 1000
  if (!run_iter %in% c(68, 220)){ #bad cells 47
  #if (!run_iter %in% c(70, 331, 338, 408, 409, 410, 411, 412, 413, 414, 415, 416, 418, 422, 424, 425, 426, 427, 428, 429, 430, 431, 432, 434, 435, 436, 437, 438, 439, 444, 
  #                     447, 448, 449, 450, 452, 521, 616, 676, 690, 929, 933, 948, 955)){ #good cells 47
    print(run_iter)
    select.cell <- BAD.FACS.cells[run_iter]
    load(paste0(work.dir,   "ZINB_results_files/", BAD.FACS.cells[run_iter],"_results.rda"))
    output[[select.cell]] <- Finaldf[[select.cell]]
  }
}

output.BAD.cells.47 <- output

study.cells <- FACS.anno[FACS.anno$cluster_label == "Low Quality VISp L5 PT Ctxn3 2", "sample_name"]
study.cells <- FACS.anno[FACS.anno$cluster_label == "Low Quality VISp L6 CT Ptprt_2", "sample_name"]
study.cells <- FACS.anno[FACS.anno$cluster_label == "Low Quality VISp L5 PT Ctxn3 1", "sample_name"]
study.cells <- FACS.anno[FACS.anno$cluster_label == "Low Quality ALM L6 CT Cpa6", "sample_name"]
study.cells <- FACS.anno[FACS.anno$cluster_label == "Low Quality Meis2 Adamts19", "sample_name"]
study.cells <- FACS.anno[FACS.anno$cluster_label == "Low Quality VISp L6 CT Ptprt_1", "sample_name"]
study.cells <- FACS.anno[FACS.anno$cluster_label == "Low Quality Astro Aqp4", "sample_name"]
study.cells <- FACS.anno[FACS.anno$cluster_label == "Low Quality Sst Chodl", "sample_name"]

output <- output.BAD.cells.47
select.cell <- study.cells[[1]][3]
output[[select.cell]][order(output[[select.cell]][,"logp"], decreasing = TRUE),][1:10,]

output <- output.BAD.cells.250
select.cell <- study.cells[[1]][3]
output[[select.cell]][order(output[[select.cell]][,"logp"], decreasing = TRUE),][1:10,]

output <- output.BAD.cells.1000
output[[select.cell]][order(output[[select.cell]][,"logp"], decreasing = TRUE),][1:10,]


study.cells <- GOOD.FACS.cells
select.cell <- "LS-15003_S52_E1-50"
output <- output.GOOD.cells.47
output[[select.cell]][order(output[[select.cell]][,"logp"], decreasing = TRUE),][1:10,]

output <- output.GOOD.cells.250
output[[select.cell]][order(output[[select.cell]][,"logp"], decreasing = TRUE),][1:10,]

output <- output.GOOD.cells.1000
output[[select.cell]][order(output[[select.cell]][,"logp"], decreasing = TRUE),][1:10,]

#select.cell <- BAD.FACS.cells[run_iter]
#ggplot(output[[select.cell]][!is.na(output[[select.cell]][,"logp"]),], aes(logp, B, colour = pair_identity)) + 
#  geom_point(alpha = 0.4) + xlab("Log(p)") + ylab("B")

plot(c(0.72, 0.91, 0.95), type="o", col="blue", ylim = c(0,1), axes=FALSE, xlab = "", ylab="% Correct predictions for GOOD cells")
axis(1, at=1:3, lab=c("47genes","250genes","1000genes"))
axis(2, las=1)
text(x=1.3, y=0.71,labels = "0.72")
axis(2, las=1, at=4*0:g_range[2])


correct <- c()
output <- output.GOOD.cells.47
cells <- names(output)

for (c in cells){
  output[[c]] <- output[[c]][order(output[[c]][,"logp"], decreasing = TRUE),]
  correct <- c(correct, sum(output[[c]][1,"cluster_label"] == output[[c]][1,"cl1"] & 
        output[[c]][1,"cluster_label"] == output[[c]][1,"cl2"]))}

sum(correct)/length(correct)

