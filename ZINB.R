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
#robject.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/R_Object/"
#markers_path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/select.markers.rda"

tome <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/tomes/facs/mouse_V1_ALM_20180520/transcrip.tome"
# Sample annotations
# Will throw Warnings, but these are OK - just because of how NAs are stored in HDF5 files.
FACS.anno <- read_tome_anno(tome)
# Read all counts as sparse matrix
# These are stored in samples (rows) x genes (columns) format
exon_counts <- read_tome_dgCMatrix(tome = tome, target = "/data/exon")
intron_counts <- read_tome_dgCMatrix(tome = tome, target = "/data/intron")
FACS.counts <- exon_counts + intron_counts
colnames(FACS.counts) <- read_tome_gene_names(tome)
rownames(FACS.counts) <- read_tome_sample_names(tome)
# See everything stored in the tome
#h5ls(tome)
print("Done!")

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

#counts <- intron[short_markers_list, FACS.cells] + exon[short_markers_list, FACS.cells]
FACS.cells <- rownames(FACS.counts)
FACS.counts <- FACS.counts[FACS.cells,short_markers_list]
rownames(FACS.anno) <- FACS.anno$sample_name
FACS.anno <- FACS.anno[FACS.cells,]
df <- cbind(FACS.anno$cluster_id, as.data.frame.matrix(FACS.counts))
colnames(df) <- c("Type", colnames(FACS.counts))


#load("/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/20190227_BT014-RSC-195_mouse_patchseq_star2.0_exon.Rdata")
#load("/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/20190227_BT014-RSC-195_mouse_patchseq_star2.0_intron.Rdata")
#load("/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/20190227_BT014-RSC-195_mouse_patchseq_star2.0_samp.dat.Rdata")
#patchseq.cells <- colnames(exon)
#patchseq_counts <- intron[short_markers_list, patchseq.cells] + exon[short_markers_list, patchseq.cells] 
#patchseq_counts <- t(patchseq_counts)
#patchseq_counts <- as.data.frame.matrix(patchseq_counts)
#names(patchseq_counts)[grepl( "Rik" , names( patchseq_counts )) ]  <- paste0("rename",names(patchseq_counts)[grep("Rik", colnames(patchseq_counts))])
#colnames(patchseq_counts) <- gsub("-", "_", colnames(patchseq_counts))
#sum(colnames(df) %in% colnames(patchseq_counts)) == length(short_markers_list)
#sample_id <- setNames(samp.dat$patched_cell_container, samp.dat$exp_component_name)
#sample_id <- sample_id[patchseq.cells]
#patchseq_counts <- cbind.data.frame(patchseq_counts[patchseq.cells, ], sample_id)
#rownames(patchseq_counts) <- patchseq_counts$sample_id
#patchseq_counts <- patchseq_counts[,short_markers_list]
#patchseq.cells <- rownames(patchseq_counts)
#patchseq_anno <- as.data.frame(read_feather(patchseq.anno.path))
#rownames(patchseq_anno) <- patchseq_anno$sample_id
#dim(patchseq_anno)
#print("Done")

##########################################################################################
### Some initialization: #################################################################
##########################################################################################
temp <- table(FACS.anno$cluster_id) >= 10 
Good_types <- as.numeric(names(temp[temp]))
Good_pairs <- t(combn(Good_types, 2))
genes <- short_markers_list
#df <- df[, 2:length(colnames(df))]
#df <- t(df)

##########################################################################################
### Functions: ###########################################################################
##########################################################################################

logit_fit <- function(model_name) {
  Pi = 1
  Alpha = NA
  Mu = 1
  result = list("Method" = "Logit" , "Mu" = Mu, "Alpha" = Alpha, "Pi" = Pi)
  #save(m1, file = paste0(work.dir,"/Zeroinfl_models/", model_name, ".rda"))
  return(result)
}

NB_fit <- function(data.df, equation, model_name){
  m1 <- glm.nb(equation, data = data.df)
  Pi = 0
  Alpha = 1. / m1$theta
  Mu = predict(m1, type = "response")[[1]]
  result = list("Method" = "NB" , "Mu" = Mu, "Alpha" = Alpha, "Pi" = Pi)
  #save(m1, file = paste0(work.dir,"/Zeroinfl_models/", model_name, ".rda"))
  return(result)
}

ZINB_fit <- function(data.df, equation, model_name){
  m1 <- zeroinfl(equation, data = data.df, dist = "negbin")
  Pi <- predict(m1, type = "zero")[[1]]
  Alpha = 1. / m1$theta
  Mu = predict(m1, type = "count")[[1]]
  #prob_of_zero = predict(m1, type = "prob")[1,1]
  result = list("Method" = "ZINB" , "Mu" = Mu, "Alpha" = Alpha, "Pi" = Pi)
  #save(m1, file = paste0(work.dir,"/Zeroinfl_models/", model_name, ".rda"))
  return(result)
}

Fit_model <- function(data.df, gene, type_id){
  model_name = paste0(type_id,"_", gene)
  sub_df <- data.df[, gene]
  if (length(sub_df) >= 10) {
    if (all(sub_df == 0)){
      result <- logit_fit(model_name)
    } else if (all(sub_df != 0)){
      eq = paste0(as.character(gene), " ~ 1")
      #print(eq)
      result <- NB_fit(data.df, equation = as.formula(eq), model_name)
    } else {
      eq = paste0(as.character(gene), " ~ 1")
      #print(eq)
      result <- try(ZINB_fit(data.df, equation = as.formula(eq), model_name))
      if("try-error" %in% class(result)) {
        print(c("Was not able to fit ZINB for gene:", gene, "Zero%", sum(sub_df == 0) / length(sub_df) ))
        result <- logit_fit(model_name)
        }
    }
    return(result)
  } else {
    result = list("Method" = "Not enough data" , "Mu" = NA, "Alpha" = NA, "Pi" = NA)
    return(result)
  }
}

##########################################################################################
### fit ZINB or NB or logit per gene per Type: ###########################################
##########################################################################################

All_fit_values <- list()
for (t in Good_types) {
  print(c("Type:", t))
  new_df <- df[df$Type == t ,]
  Fit_values <- tapply(1:nrow(new_df), as.character(new_df$Type), function(x) sapply(colnames(new_df)[-1], function(g)Fit_model(new_df[x,], g, t)))
  All_fit_values <- c(All_fit_values, Fit_values)
}
setwd("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/Patchseq_vs_FACs_cre_analysis/mouse_patchseq_VISp_20181220_collapsed40_cpm/")
save(All_fit_values, file="All_fit_values.rda")
load("All_fit_values.rda")

df <- df[,2:length(colnames(df))]
##########################################################################################
### Find Beta for each pair of cell type using an optimizer: #############################
##########################################################################################
df <- t(df)

Err <- function(data, par) {
  with(data, sum((par[1] * term1 - term2)^2))
}

Beta_list <- list()
for (c in FACS.cells[1:500]){
  Beta <- list()
  start_time = Sys.time()
  for (i in 1:dim(Good_pairs)[1]) {
    first_type <- as.character(Good_pairs[i,1])
    second_type <- as.character(Good_pairs[i,2])
    pair_name <- paste0(as.character(first_type), "_", as.character(second_type))
    list1 <- as.data.frame(t(All_fit_values[[first_type]]))
    list2 <- as.data.frame(t(All_fit_values[[second_type]]))
    #y = patchseq_counts[genes, c] 
    y = df[genes, c]
    if(length(y) == 0) {print("ERROR: Something is wrong!!!!!!!")}
    Mu1 <- unlist(list1[genes, "Mu"])
    Mu2 <- unlist(list2[genes, "Mu"])
    Pi1 <- unlist(list1[genes, "Pi"])
    Pi2 <- unlist(list2[genes, "Pi"])
    term1 <- Mu1 - Mu1 * Pi1 - Mu2 + Mu2 * Pi2
    term2 <- y - Mu2 + Mu2 * Pi2
    dat <- as.data.frame(cbind(term1, term2))
    results <- optim(par = runif(1, 0.51, 0.99), 
                     fn = Err, 
                     data = dat, 
                     method = "Brent", 
                     lower = 0, 
                     upper = 1000)
    Beta[pair_name] <-  results$par
  } 
  Beta_list[[c]] <- Beta
  end_time = Sys.time()
  print(end_time - start_time)
}
print("Done!")
save(Beta_list, file = "Beta_list.rda")
#load("Beta_list.rda")

##########################################################################################
### Find the probability of the observation given Beta: ##################################
##########################################################################################

Compute_likelihood_efficient <- function(contam_list, Method, Mu, Pi, Alpha) {
  if(Method == "NB"){
    n <- 1/Alpha
    prob <- n /(n + Mu)
    L <- lapply(contam_list, function(x){dnbinom(x, n, prob, log = FALSE)})
    #print("I am here1")
    
  } else if (Method == "ZINB") {
      LNB0 <- Pi + (1 - Pi) * (1 + Alpha * Mu) ^ (-1 / Alpha)
      n <- 1/Alpha
      prob <- n /(n + Mu)  
      LNB <- lapply(setdiff(contam_list, 0), function(x){dnbinom(x, n, prob, log = FALSE)})
      LNB <- lapply(LNB, function(x){(1-Pi) * x})
      L <- c(LNB0 , LNB)
      #print("I am here2")
      
  } else if (Method == "Logit") {
    L <- c(1, rep(0, length(contam_list) - 1))
    #print("I am here3")
    
  }
  L
}

Get_contamination_factor <- function(cell_name, pair_name){
  Beta_list[[c]][[pair_name]]
}
  
Get_contamination_counts <- function(beta, yg){
  contamin_count1 <- list()
  contamin_count2 <- list()
  for (k in 0:yg){
    contamin_count1 <- c(contamin_count1, round(k / beta))
    contamin_count2 <- c(contamin_count2, round(( yg - k) / (1 - beta)))
  }
  result <- list("contamin_count1" = contamin_count1 , "contamin_count2"= contamin_count2)
  return(result)
}

Compute_all_likelihood <- function(gene, yg, list1, list2, cell_name, pair_name) {
  beta <- Get_contamination_factor(cell_name, pair_name)
  con <- Get_contamination_counts(beta, yg)
  con1 <- unlist(con$contamin_count1) 
  con2 <- unlist(con$contamin_count2) 
  #print(con1)
  #print(con2)
  term1 <- unlist(Compute_likelihood_efficient(con1, Method = list1[[gene,"Method"]], Pi = list1[[gene,"Pi"]], Mu = list1[[gene,"Mu"]], Alpha = list1[[gene,"Alpha"]]))
  term2 <- unlist(Compute_likelihood_efficient(con2 , Method = list2[[gene,"Method"]], Pi = list2[[gene,"Pi"]], Mu = list2[[gene,"Mu"]], Alpha = list2[[gene,"Alpha"]]))
  sum(term1 * rev(term2))
}

All_cell_pyg <- list()
for (c in FACS.cells[1:500]){
  start_time = Sys.time()
  pyg <- list()
  for (i in 1:dim(Good_pairs)[1]) {
    #print(i)
    first_type <- as.character(Good_pairs[i,1])
    second_type <- as.character(Good_pairs[i,2])
    pair_name <- paste0(as.character(first_type), "_", as.character(second_type))
    list1 <- as.data.frame(t(All_fit_values[[first_type]]))
    list2 <- as.data.frame(t(All_fit_values[[second_type]]))
    #y <- patchseq_counts[, c]
    y <- df[,c]
    pyg[[pair_name]] <- lapply(names(y), function(x){Compute_all_likelihood(x, y[[x]], list1, list2, c, pair_name)})
  }
  All_cell_pyg[[c]] <- pyg
  end_time = Sys.time()
  print(end_time - start_time)
}

save(All_cell_pyg, file = "All_cell_pyg.rda")

##########################################################################################
### Find the probability of the observation given Beta: ##################################
##########################################################################################

results <- list()
for (c in FACS.cells[7]){
  for (i in 1:dim(Good_pairs)[1]) {
    first_type <- as.character(Good_pairs[i,1])
    second_type <- as.character(Good_pairs[i,2])
    pair_name <- paste0(as.character(first_type), "_", as.character(second_type))
    print(pair_name)
    results <- c(results , sum(log(unlist(All_cell_pyg[[c]][[pair_name]]))))
  }
}

results[is.infinite(unlist(results))] <- -1000 
max(unlist(results))
Good_pairs[which.max(unlist(results)),]
VIS.cldf[VIS.cldf$cluster_id == Good_pairs[which.max(unlist(results)),][1],c("cluster_label")]
VIS.cldf[VIS.cldf$cluster_id == Good_pairs[which.max(unlist(results)),][2],c("cluster_label")]
VIS.cldf[VIS.cldf$cluster_id == VIS.cl[FACS.cells[7]],c("cluster_label")]
Ty
