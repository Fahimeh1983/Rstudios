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

contaminated_types <- c("Doublet SMC and Glutamatergic", "Doublet Astro Aqp4 Ex",
                        "Doublet Endo and Peri_1", "Doublet VISp L5 NP and L6 CT", "Doublet Endo Peri SMC")

highQ_type <- setdiff(unique(FACS.anno$cluster_label), LowQ_types)
cluster_label_id <- as.data.frame(unique(FACS.anno[,c("cluster_label", "cluster_id")]))
FACS.anno <- as.data.frame(FACS.anno) 
rownames(FACS.anno) <- FACS.anno$sample_name
load(paste0(work.dir, "/Test_GOOD_cells.rda"))
load(paste0(work.dir, "/Test_BAD_cells.rda"))
load(paste0(work.dir, "/validation_cells.rda"))
sum(validation.cells %in% c(test.GOOD.cells, test.BAD.cells)) == 0
sum(test.BAD.cells %in% c(test.GOOD.cells, validation.cells)) == 0
sum(test.GOOD.cells %in% c(validation.cells, test.BAD.cells)) == 0
test_and_validation_cells <- c(test.BAD.cells, test.GOOD.cells, validation.cells)

load(paste0(work.dir, "/Good_trained_types.rda"))
load(paste0(work.dir, "/Good_trained_pairs.rda"))

#FACS.anno <- FACS.anno[test_and_validation_cells,]
#colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ]  <- paste0("rename",colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ])
#colnames(FACS.counts) <- gsub("-", "_", colnames(FACS.counts))
#FACS.counts <- FACS.counts[test_and_validation_cells, Long_markers_list]
#sum(rownames(FACS.anno) == rownames(FACS.counts)) == length(test_and_validation_cells)
#df <- cbind(FACS.anno$cluster_id, as.data.frame.matrix(FACS.counts))
#colnames(df) <- c("Type", colnames(FACS.counts))
#dim(df)

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
  
  for (c in test_and_validation_cells[run_iter]){
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
  
  count_threshold <- 200000
  pyg <- list()
  for (c in test_and_validation_cells[run_iter]){
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
  
  #unlist(pp)[order(unlist(l))]
  ##########################################################################################
  ### Find the probability of the observation given Beta: ##################################
  ##########################################################################################
  
  cell <- test_and_validation_cells[run_iter]

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
      cl1 <- c(cl1, cluster_label_id[cluster_label_id$cluster_id == first_type, "cluster_label"])
      cl2 <- c(cl2, cluster_label_id[cluster_label_id$cluster_id == second_type, "cluster_label"])
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
### Looking at the results: ##############################################################
##########################################################################################
all_test_cells <- union(test.BAD.cells, test.GOOD.cells)
validation.cells <- validation.cells
contaminated.cells <- FACS.anno[FACS.anno$cluster_label %in% contaminated_types, "sample_name"]

test_output <- list()
validation_output <- list()
done_test_cells <- c()
done_validation_cells <- c()
for (run_iter in 1:5089){
  cell_name <- test_and_validation_cells[run_iter]
  file_name <- paste0(work.dir,   "ZINB_results_files/", cell_name,"_results.rda")
  if (file.exists(file_name)) {
    if (cell_name %in% all_test_cells){
      done_test_cells <- c(done_test_cells, cell_name)
      load(file_name)
      test_output[[cell_name]] <- Finaldf[[cell_name]]
      }
    else if (cell_name %in% validation.cells) { 
      done_validation_cells <- c(done_validation_cells, cell_name)
      load(file_name)
      validation_output[[cell_name]] <- Finaldf[[cell_name]]
      }
    else{print("Something is wrong")}
  }
}

# computing the accuracy on the good cells, we remove the poorQ cells for this part 
correct <- c()
combined_output <- c(test_output, validation_output)
for (c in setdiff(c(done_validation_cells, done_test_cells), test.BAD.cells)){
  combined_output[[c]] <- combined_output[[c]][order(combined_output[[c]][,"logp"], decreasing = TRUE),]
  tmp <- combined_output[[c]][1,"cluster_label"] == combined_output[[c]][1,"cl1"] & 
        combined_output[[c]][1,"cluster_label"] == combined_output[[c]][1,"cl2"]
  correct <- c(correct, tmp)
}

sum(correct)/length(correct) #80% accuarcy on the good cells

# Now we should compute the mean of logp for the validation set
validation_best_logp <- data_frame()
for (c in done_validation_cells) {
  validation_output[[c]] <- validation_output[[c]][order(validation_output[[c]][,"logp"], decreasing = TRUE),]
  if (validation_output[[c]][, "cluster_id"] == validation_output[[c]][,"Type1"] & 
      validation_output[[c]][, "cluster_id"] == validation_output[[c]][,"Type2"]) {
    validation_best_logp <- rbind(validation_best_logp, validation_output[[c]][1,])
  }
}


validation_mean_logp <- as.data.frame(validation_best_logp %>% 
                                        group_by(cluster_id) %>% 
                                        summarise(mean_logp = mean(logp), min_logp = min(logp), max_logp=max(logp)))


cell_identity <- function(cell, test.BAD.cells, test.GOOD.cells, contaminated.cells){
  if(cell %in% test.GOOD.cells) { return(c("PURE"))}
  if(cell %in% setdiff(test.BAD.cells, contaminated.cells)) {return(c("POORQ"))}
  if(cell %in% contaminated.cells){return(c("CONTAMIN"))}
}

prediction_identity <- function(output_1line, cellQ, validation_mean_logp){
  if (cellQ == "PURE"){
    if (output_1line$Type1 == output_1line$Type2 & output_1line$cl1 == output_1line$cluster_label){return(c("PURE_PURE_CORRECT"))}
    if (output_1line$Type1 == output_1line$Type2 & output_1line$cl1 != output_1line$cluster_label){return(c("PURE_PURE_WRONG"))}
    if (output_1line$Type1 != output_1line$Type2){return(c("PURE_CONTAMIN"))}
  } else if (cellQ == "POORQ") {
    if (output_1line$Type1 == output_1line$Type2){return(check_prediction_quality(cellQ, output_1line, validation_mean_logp))}
    if (output_1line$Type1 != output_1line$Type2){return(c("POORQ_CONTAMIN"))}
  } else if (cellQ == "CONTAMIN") {
    if (output_1line$Type1 == output_1line$Type2){return(check_prediction_quality(cellQ, output_1line, validation_mean_logp))}
    if (output_1line$Type1 != output_1line$Type2){return(c("CONTAMIN_CONTAMIN"))}
  }
}

check_prediction_quality <- function(cellQ, output_1line, validation_mean_logp){
  if (cellQ == "POORQ"){
    if (!is_empty(validation_mean_logp[validation_mean_logp$cluster_id == output_1line$Type1, "min_logp"])){
      if (output_1line$logp >= validation_mean_logp[validation_mean_logp$cluster_id == output_1line$Type1, "min_logp"]){
        return(c("POORQ_PURE"))
      } else {return(c("POORQ_POORQ"))}
    }else{return(NaN)}
  }
  if (cellQ == "CONTAMIN"){
    if (!is_empty(validation_mean_logp[validation_mean_logp$cluster_id == output_1line$Type1, "min_logp"])){
      if (output_1line$logp >= validation_mean_logp[validation_mean_logp$cluster_id == output_1line$Type1, "min_logp"]){
        return(c("CONTAMIN_PURE"))
      } else {return(c("CONTAMIN_POORQ"))}
    } else{return(NaN)}
  }
}

prediction <- c()
for (c in done_test_cells){
  cellQ <- cell_identity(c, test.BAD.cells, test.GOOD.cells, contaminated.cells)
  tmp <- test_output[[c]][order(test_output[[c]][,"logp"], decreasing = TRUE),][1,]
  prediction <- c(prediction, prediction_identity(tmp, cellQ, validation_mean_logp))
}

sum(prediction == "PURE_PURE_CORRECT") / sum(prediction == "PURE_PURE_CORRECT" | prediction == "PURE_PURE_WRONG" | prediction == "PURE_PURE_CONTAMIN")
sum(prediction == "PURE_PURE_WRONG") / sum(prediction == "PURE_PURE_CORRECT" | prediction == "PURE_PURE_WRONG" | prediction == "PURE_PURE_CONTAMIN")
sum(prediction == "POORQ_POORQ") / sum(prediction == "POORQ_POORQ" | prediction == "POORQ_PURE" | prediction == "POORQ_CONTAMIN")
sum(prediction == "CONTAMIN_CONTAMIN") / sum(prediction == "CONTAMIN_PURE" | prediction == "CONTAMIN_POORQ" | prediction == "CONTAMIN_CONTAMIN")

temp <- c()
for (c in done_test_cells[prediction == "POORQ_PURE"]) {
  temp <- rbind(temp, test_output[[c]][order(test_output[[c]][,"logp"], decreasing = TRUE),][1,])
}
temp 

ggplot(melt(cbind(Good_test_cells, Bad_test_cells, contaminated_test_cells)), aes(value  , fill = Var2)) + 
  geom_density(alpha = 0.5) + xlab("Min(logp) - predicted(logp)") + ylab("Density")

