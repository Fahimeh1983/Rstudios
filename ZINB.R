##########################################################################################
### Reading the counts: ##################################################################
##########################################################################################
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/scrattch.hicat/R/dend.markers.R")
#devtools::install_github("AllenInstitute/scrattch.io", ref = "dev")
require(ggplot2)
require(pscl)
require(MASS)
require(boot)
library(dplyr)
library(tidyr)
library(feather)
library(scrattch.io)
library(scrattch.hicat)

work.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/Patchseq_vs_FACs_cre_analysis/mouse_patchseq_VISp_20181220_collapsed40_cpm/"
#tome <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/tomes/facs/mouse_V1_ALM_20180520/transcrip.tome"
cpm <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/data.feather"
#markers.path <- paste0(work.dir, "/1000_marker_genes_include47.rda")
markers.path <- paste0(work.dir, "/select.markers.rda")

# Sample annotations
# Will throw Warnings, but these are OK - just because of how NAs are stored in HDF5 files.
cpm_counts <- read_feather(cpm)
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
#load(paste0(work.dir, "/Counts/FACS_counts.rda"))
load(markers.path)

FACS.anno <- as.data.frame(FACS.anno) 
Long_markers_list <- select.markers
cpm_counts <- cpm_counts[,c("sample_id",Long_markers_list)]
cpm_counts <- as.data.frame(cpm_counts)
rownames(cpm_counts) <- cpm_counts$sample_id
cpm_counts <- cpm_counts[,setdiff(colnames(cpm_counts), "sample_id")]
FACS.counts <- cpm_counts

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
load(paste0(work.dir, "/train_cells.rda"))
sum(validation.cells %in% c(test.GOOD.cells, test.BAD.cells)) == 0
sum(test.BAD.cells %in% c(test.GOOD.cells, validation.cells)) == 0
sum(test.GOOD.cells %in% c(validation.cells, test.BAD.cells)) == 0
test_and_validation_cells <- c(test.BAD.cells, test.GOOD.cells, validation.cells)
#FACS.counts <- FACS.counts[test_and_validation_cells, Long_markers_list]

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

FACS.anno <- FACS.anno[test_and_validation_cells,]
colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ]  <- paste0("rename",colnames(FACS.counts)[grepl( "Rik" , colnames(FACS.counts)) ])
colnames(FACS.counts) <- gsub("-", "_", colnames(FACS.counts))
Long_markers_list[grepl( "Rik" , Long_markers_list)] <-  paste0("rename", Long_markers_list[grepl( "Rik" , Long_markers_list) ])
Long_markers_list <- gsub("-", "_", Long_markers_list)
FACS.counts <- FACS.counts[test_and_validation_cells, Long_markers_list]
sum(rownames(FACS.anno) == rownames(FACS.counts)) == length(test_and_validation_cells)
df <- cbind(FACS.anno$cluster_id, as.data.frame.matrix(FACS.counts))
colnames(df) <- c("Type", colnames(FACS.counts))
tmp<- rownames(df)
df <- cbind.data.frame(lapply(df, as.integer))
rownames(df) <- tmp
dim(df)

##########################################################################################
### Some initialization: #################################################################
##########################################################################################

genes <- Long_markers_list
df <- t(df)
load(paste0(work.dir, "/All_fit_values_4020.rda"))

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
#contaminated.cells <- contaminated.cells$sample_name

test_output <- list()
validation_output <- list()
done_test_cells <- c()
done_validation_cells <- c()
for (run_iter in 1:5100){
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

prediction_identity <- function(output_1line, cell, cellQ, validation_mean_logp){
  Type1 <- as.character(output_1line$Type1) 
  Type2 <- as.character(output_1line$Type2) 
  cluster1 <- as.character(output_1line$cl1)
  cluster2 <- as.character(output_1line$cl2)
  label <- as.character(output_1line$cluster_label)
  B <- output_1line$B
  if (cellQ == "PURE"){
    if (Type1 == Type2 & cluster1 %in% highQ_type){
      if(cluster1 == label){
        return(c(cell, Type1, Type2, cluster1, cluster2, B, label, "PURE_PURE_CORRECT"))
      }else{
        return(c(cell, Type1, Type2, cluster1, cluster2, B, label, "PURE_PURE_WRONG"))
      }
    }
    if (Type1 != Type2){return(c(cell, Type1, Type2, cluster1, cluster2, B, label, "PURE_CONTAMIN"))}
  } else if (cellQ == "POORQ") {
    if (Type1 == Type2 & cluster1 %in% highQ_type){
      if (length(validation_mean_logp[validation_mean_logp$cluster_id == Type1, "min_logp"])!=0){
              if (output_1line$logp >= validation_mean_logp[validation_mean_logp$cluster_id == Type1, "min_logp"]){
                return(c(cell, Type1, Type2, cluster1, cluster2, B, label,"POORQ_PURE_CORRECT"))
              } else {return(c(cell, Type1, Type2, cluster1, cluster2, B, label, "POORQ_PURE_WRONG"))}
            }else{return(c(cell, Type1, Type2, cluster1, cluster2, B, label,NaN))}
    }
    if (Type1 != Type2){return(c(cell, Type1, Type2, cluster1, cluster2, B, label,"POORQ_CONTAMIN"))}
  } else if (cellQ == "CONTAMIN") {
    if (Type1 == Type2 & cluster1 %in% highQ_type){
      if (length(validation_mean_logp[validation_mean_logp$cluster_id == Type1, "min_logp"])!=0){
        if (output_1line$logp >= validation_mean_logp[validation_mean_logp$cluster_id == Type1, "min_logp"]){
          return(c(cell, Type1, Type2, cluster1, cluster2, B, label,"CONTAMIN_PURE_CORRECT"))
        } else {return(c(cell, Type1, Type2, cluster1, cluster2, B, label,"CONTAMIN_PURE_WRONG"))}
      }else{return(c(cell, Type1, Type2, cluster1, cluster2, B, label,NaN))}  
    }
    if (Type1 != Type2){return(c(cell, Type1, Type2, cluster1, cluster2, B, label,"CONTAMIN_CONTAMIN"))}
  }
}


prediction <- c()
for (c in done_test_cells){
  cellQ <- cell_identity(c, test.BAD.cells, test.GOOD.cells, contaminated.cells)
  tmp <- test_output[[c]][order(test_output[[c]][,"logp"], decreasing = TRUE),][1,]
  prediction <- rbind(prediction, prediction_identity(tmp,c , cellQ, validation_mean_logp))
}
prediction <- as.data.frame(prediction)
colnames(prediction) <- c( "sample_id", "Type1", "Type2","cl1", "cl2", "B", "cluster_label", "identity")

sum(prediction$identity == "PURE_PURE_CORRECT") 
sum(prediction$identity == "PURE_PURE_WRONG") 
sum(prediction$identity == "PURE_CONTMAIN") 
sum(prediction$identity == "POORQ_PURE_CORRECT")
sum(prediction$identity == "POORQ_PURE_WRONG")
sum(prediction$identity == "POORQ_CONTAMIN")
sum(prediction$identity == "CONTAMIN_PURE_CORRECT")
sum(prediction$identity == "CONTAMIN_PURE_WRONG")
sum(prediction$identity == "CONTAMIN_CONTAMIN")

study.cells <- done_test_cells[prediction$identity == "CONTAMIN_CONTAMIN"]
for (c in study.cells){
  tmp <- test_output[[c]][order(test_output[[c]][,"logp"], decreasing = TRUE),][1,]
  print(c(as.character(tmp$cl1), as.character(tmp$cl2), tmp$B, tmp$cluster_label))
}

##########################################################################################
### Reading in mapping results of FACS onto V1+ALM types: ################################
##########################################################################################
ref.data.rda.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2/")
load(paste0(ref.data.rda.path, "norm.dat.rda"))
dim(norm.dat)
norm.dat <- norm.dat[,rownames(FACS.anno)]
load(paste0(work.dir, "V1_ALM_dend_markers_attached.rda") )
load(paste0(work.dir, "/FACS_ALM_V1_TREE.rda"))
load(paste0(work.dir, "/FACS_ALM_V1_TREE_mapping.df.rda"))

cl <- setNames(FACS.anno$cl, FACS.anno$sample_name)
cl <- factor(cl)
cl.df <- unique(as.data.frame(FACS.anno[,c("cluster_id", "cluster_label", "cluster_color", "subclass_id", "subclass_label", "class_id", "class_label", "cl")]))
rownames(cl.df) <- cl.df$cl
cl.df <- cl.df[order(cl.df$cluster_id),]
cl.df <- cl.df[cl.df$class_label!="Low Quality",]
cl <- cl[cl%in% cl.df$cl]
cl <- droplevels(cl)


##########################################################################################
### Comparing mapping and clustering results: ############################################
##########################################################################################
#FACS anno lables come from clustering
#FACS.mapping.df is from mapping

rownames(prediction) <- prediction$sample_id
sum(rownames(FACS.anno) %in% rownames(FACS_Tree_mapping.df))
FACS_Tree_memb <- FACS_Tree_memb[rownames(FACS_Tree_mapping.df),]
sorted_cl <- t(apply(FACS_Tree_memb[ , labels(dend)], 1, function(x) names(x)[order(x, decreasing = TRUE)]))
sorted_bt <- t(apply(FACS_Tree_memb[ , labels(dend)], 1, function(x) x[order(x, decreasing = TRUE)]))

first_cl <- Renew_list(ls = sorted_cl[,1], ref.df = cl.df, label = "cl", new.label = "cluster_label")
second_cl <- Renew_list(ls = sorted_cl[,2], ref.df = cl.df, label = "cl", new.label = "cluster_label")

first_cl <- first_cl[rownames(FACS_Tree_mapping.df)] 
second_cl <- second_cl[rownames(FACS_Tree_mapping.df)]

first_bt <- sorted_bt[rownames(FACS_Tree_mapping.df),1]
second_bt <- sorted_bt[rownames(FACS_Tree_mapping.df),2]
  
FACS_Tree_mapping.df <- cbind(FACS_Tree_mapping.df, first_cl, second_cl, first_bt, second_bt)

study.cells <- done_test_cells[prediction$identity == "PURE_PURE_CORRECT" | prediction$identity == "PURE_PURE_WRONG" | prediction$identity == "PURE_PURE_CONTAMIN"]
#Accuracy on the good cell, ZINB model accuarcy was 86%
sum(FACS.anno[study.cells, "cluster_label"] == FACS_Tree_mapping.df[study.cells, "cluster_label"]) / length(study.cells)
sum(as.character(FACS.anno[study.cells, "cluster_label"]) == FACS_Tree_mapping.df[study.cells, "first_cl"]) / length(study.cells)

study.cells <- done_test_cells[prediction$identity == "POORQ_PURE_CORRECT" | prediction$identity == "POORQ_PURE_WRONG" | prediction$identity == "POORQ_CONTAMIN"]
tmp <- cbind(FACS_Tree_mapping.df[study.cells, c("cluster_label", "first_cl", "second_cl", "first_bt", "second_bt") ], prediction[study.cells, c("cl1", "cl2", "B", "identity")])
sum(FACS_Tree_mapping.df[study.cells, "cluster_label"] %in% highQ_type)

study.cells <- done_test_cells[prediction$identity == "CONTAMIN_PURE_CORRECT" | prediction$identity == "CONTAMIN_PURE_WRONG" | prediction$identity == "CONTAMIN_CONTAMIN"]
tmp <- cbind(FACS_Tree_mapping.df[study.cells, c("cluster_label", "first_cl", "second_cl", "first_bt", "second_bt") ], prediction[study.cells, c("cl1", "cl2", "B", "identity")])
tmp[tmp$identity == "CONTAMIN_CONTAMIN",]

##########################################################################################
### D-contamination and new mapping: #####################################################
##########################################################################################

dim(prediction)
prediction$Type1 <- as.numeric(prediction$Type1)
prediction$Type2 <- as.numeric(prediction$Type2)
prediction$B <- as.numeric(prediction$B)
prediction <- prediction %>% mutate(Final_cl = ifelse(B > 0.5, cl1, cl2)) 

genes <- Long_markers_list

new_df <- c()
for (c in (prediction$sample_id)) {
  BB <- prediction[c, "B"]
  contaminant <- ifelse( BB < 0.5, prediction[c, "Type1"], prediction[c, "Type2"])
  y = df[genes, c]
  new_y <- c()
  for (g in genes) {
    PP <- All_fit_values[[as.character(contaminant)]][,g][["Pi"]]
    MM <- All_fit_values[[as.character(contaminant)]][,g][["Mu"]]
    new_y[g] <- (1./ BB) * ( y[g] - (1 - BB) * (1 - PP) * MM) 
  }
  new_y[new_y<0] <- 0
  new_df <- cbind(new_df, new_y)
}

colnames(new_df) <- as.character(prediction$sample_id)
new_norm.dat <- log2(new_df+1)
norm.dat <- norm.dat[,names(cl)] 
norm.dat <- as.matrix(norm.dat)

all_test_cells <- all_test_cells[all_test_cells %in% colnames(new_norm.dat)]
#new_cl <- cl[all_test_cells]
class(new_norm.dat)

rownames(new_norm.dat)[grepl( "rename" , rownames(new_norm.dat)) ]  <- gsub("rename","",rownames(new_norm.dat)[grepl( "rename" , rownames(new_norm.dat)) ])
rownames(new_norm.dat) <- gsub("_", "-", rownames(new_norm.dat))
new_FACS_Tree_memb <- map_dend_membership(dend, cl= cl, norm.dat, new_norm.dat, colnames(new_norm.dat), bs.num=100, p=0.7, low.th=0.15)

##########################################################################################
### Comparing the mapping of contaminated and d-contaminated data: #######################
##########################################################################################
library(reshape2)
study.cells <- contaminated.cells
study.cells <- done_test_cells[prediction$identity == "CONTAMIN_CONTAMIN"]

old_sorted_bt <- t(apply(FACS_Tree_memb[study.cells, labels(dend)], 1, function(x) x[order(x, decreasing = TRUE)]))
old_first_bt <- old_sorted_bt[study.cells,1]
old_second_bt <- old_sorted_bt[study.cells, 2]

new_sorted_bt <- t(apply(new_FACS_Tree_memb[study.cells, labels(dend)], 1, function(x) x[order(x, decreasing = TRUE)]))
new_first_bt <- new_sorted_bt[study.cells, 1]
new_second_bt <- new_sorted_bt[study.cells, 2]

old_sorted_cl <- t(apply(FACS_Tree_memb[study.cells, labels(dend)], 1, function(x) names(x)[order(x, decreasing = TRUE)]))
old_first_cl <- old_sorted_cl[study.cells, 1]
old_second_cl <- old_sorted_cl[study.cells, 2]
old_first_cluster_label <- Renew_list(old_first_cl, cl.df, label = "cl", new.label = "cluster_label")

new_sorted_cl <- t(apply(new_FACS_Tree_memb[study.cells, labels(dend)], 1, function(x) names(x)[order(x, decreasing = TRUE)]))
new_first_cl <- new_sorted_cl[study.cells, 1]
new_second_cl <- new_sorted_cl[study.cells, 2]
new_first_cluster_label <- Renew_list(new_first_cl, cl.df, label = "cl", new.label = "cluster_label")

tmp <- cbind.data.frame(old_first_bt, new_first_bt, old_second_bt, new_second_bt)
tt <- melt(tmp)
ggplot(tt[tt$variable %in% c("old_first_bt", "new_first_bt"),], aes (value, fill = variable)) +
  geom_density(alpha=0.5)
ggplot(tt[tt$variable %in% c("old_second_bt", "new_second_bt"),], aes (value, fill = variable)) +
  geom_density(alpha=0.5)

prediction_Type1 <- Renew_list(set_names(prediction$Type1, rownames(prediction)), ref.df = cl.df, label = "cluster_id", new.label = "cl")
prediction_Type1<- prediction_Type1[study.cells]
prediction_cluster_label1 <- set_names(prediction$Final_cl, prediction$sample_id)
prediction_cluster_label1 <- prediction_cluster_label1[study.cells]

tmp <- cbind.data.frame(old_first_cl, new_first_cl, prediction_Type1 )
sum(tmp$old_first_cl == tmp$prediction_Type1)

tmp <- cbind.data.frame(old_first_cluster_label, new_first_cluster_label, prediction_cluster_label1 )
sum(tmp$old_first_cluster_label == tmp$prediction_cluster_label1 )
