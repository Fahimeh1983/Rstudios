##########################################################################################
### Setting up the libraries: ############################################################
##########################################################################################

library(dendextend)
library(matrixStats)
library(feather)
library(dplyr)
library(scrattch.io)
library(scrattch.vis)
library(feather)
library(dendextend)
library(tibble)
library(Matrix)
Where = "/allen"
mydir = paste0(Where,"/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/Patchseq_vs_FACs_cre_analysis/mouse_patchseq_VISp_20181220_collapsed40_cpm/")
source(paste0(Where, "/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/heatmap.R"))
source(paste0(Where,"/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/de.genes.R"))
source(paste0(Where,"/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/dendro.R"))
source(paste0(Where,"/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/patchseq.R"))
source(paste0(Where,"/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Plot_tree_functions.R"))
source(paste0(Where,"/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R"))
setwd(mydir)
#devtools::install_github("AllenInstitute/scrattch.io")
options(stringsAsFactors = F)


##########################################################################################
### Setting up the path dirs: ############################################################
##########################################################################################

batch_date="20190722_BT014-RSC-215"
patchseq.dir =  "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/"
facs.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"
robject.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
FACs.robject.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/R_Object/"
ref.data.rda.path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"
select.markers.path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/select.markers.rda"
latest.mapping.memb.path ="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20190723_collapsed40_cpm/mapping.memb.with.bp.40.rda"
latest.mapping.df.path ="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20190723_collapsed40_cpm/mapping.df.with.bp.40.rda"
work.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/"

##########################################################################################
### Reading REF FACS data, dend and markers list: ########################################
##########################################################################################

FACs_anno <- Read_anno_cellonly_region(paste0(facs.dir, "anno.feather"), region = "VISp")
load(paste0(ref.data.rda.path, "ref.data.rda"))
load(paste0(ref.data.rda.path, "V1.dend.with.bp.40.rda"))
load(paste0(ref.data.rda.path, "V1.dend.list.with.bp.40.rda"))
plot(dend)
cl <- Renew_list(ls = cl, ref.df = cl.df, label = "cluster_id", new.label = "cluster_label")
select.cl <- labels(dend)
load(select.markers.path)

######################################################################################################
### Loading the query data ###########################################################################
######################################################################################################

tmp<-load(paste0(robject.dir,batch_date,"_mouse_patchseq_star2.0_cpm.Rdata"))
query.dat = cpmR

# loading samp.dat object
tmp<-load(paste0(robject.dir,batch_date,"_mouse_patchseq_star2.0_samp.dat.Rdata"))

keepcells = which(samp.dat$Region=="VISp" & samp.dat$Type=="patch_seq")
samp.dat = samp.dat[c(keepcells, which(samp.dat$Region=="TCx"),which(samp.dat$Region=="FCx"),which(samp.dat$Region=="MOp"),which(samp.dat$Region=="TEa")   ),]   #FCx is for Brian.  Rat samples mapped in mouse

query.dat = query.dat[,as.character(samp.dat$exp_component_name)]
colnames(query.dat)=as.character(samp.dat$patched_cell_container)

query.dat.norm = log2(as.matrix(query.dat+1))
idx=match(rownames(norm.dat), rownames(query.dat.norm))
query.dat.norm=query.dat.norm[idx,]

patchseq_anno <- Read_patchseq_anno(paste0(patchseq.dir,"anno.feather"))
sum(select.markers %in% rownames(query.dat.norm))
 
#Patchseq Cells of interests
locked_cells_spec_id = rownames(read.csv(paste0(work.dir, "mouse_met_Jun_14.csv"), check.names=FALSE, row.names = 1 ))
locked_cells_sample_id = patchseq_anno[patchseq_anno$spec_id_label %in% locked_cells_spec_id, "sample_id"]
length(locked_cells_spec_id) == length(locked_cells_sample_id)
dim(query.dat.norm[select.markers, locked_cells_sample_id])
#write.csv(query.dat.norm, file=paste0(work.dir, "/query_dat_norm.csv"))

##########################################################################################
######################################## Mapping data using Tree #########################
##########################################################################################

set.seed(1983)
#Patchseq_Tree_memb = map_dend_membership(dend, cl, norm.dat, query.dat.norm, colnames(query.dat.norm), bs.num=100, p=0.7, low.th=0.15)
#Patchseq_Tree_mapping.df <- summarize_cl(dend, Patchseq_Tree_memb, query.dat.norm, conf.th=0.7, min.genes=1, min.genes.ratio=0.3)
#save(Patchseq_Tree_mapping.df, file=file.path(paste0("Final_matrices/Patchseq_Tree_mapping.df", ".rda"))
#save(Patchseq_Tree_memb, file=file.path(paste0("Final_matrices/Patchseq_Tree_mapping.memb",".rda"))
#Patchseq_Tree_memb <- Patchseq_Tree_memb[locked_cells_sample_id,]
#Patchseq_Tree_mapping.df <- Patchseq_Tree_mapping.df[locked_cells_sample_id,]
#save(Patchseq_Tree_mapping.df, file=file.path(paste0("Final_matrices_locked_data/Patchseq_Tree_mapping.df", ".rda"))
#save(Patchseq_Tree_memb, file=file.path(paste0("Final_matrices_locked_data/Patchseq_Tree_mapping.memb",".rda"))
#load("Final_matrices_locked_data/AIM1.1/mapping.memb.rda")
#load("Final_matrices_locked_data/AIM1.1/mapping.df.rda")
load(latest.mapping.memb.path)
Patchseq_Tree_memb <- memb
load(latest.mapping.df.path)
Patchseq_Tree_mapping.df <- mapping.df
rm(memb)
rm(mapping.df)
Patchseq_Tree_memb <- Patchseq_Tree_memb[locked_cells_sample_id,]
Patchseq_Tree_mapping.df <- Patchseq_Tree_mapping.df[locked_cells_sample_id,]

#Mapping FACS data
FACs.cells <- colnames(norm.dat)
set.seed(1983)
#FACS_Tree_memb <- map_dend_membership(dend, cl, norm.dat, norm.dat, FACs.cells, bs.num=100, p=0.7, low.th=0.15)
#FACS_Tree_mapping.df <- summarize_cl(dend, FACS_Tree_memb, norm.dat, conf.th=0.7, min.genes=1, min.genes.ratio=0.3)
#save(FACS_Tree_mapping.df, file=file.path(paste0("Final_matrices_locked_data/FACS_Tree_mapping.df", ".rda")))
#save(FACS_Tree_memb, file=file.path(paste0("Final_matrices_locked_data/FACS_Tree_mapping.memb",".rda")))
#load("Final_matrices_locked_data/AIT2.3.1/mapping.df.rda")
#load("Final_matrices_locked_data/AIT2.3.1/mapping.memb.rda")


##########################################################################################
### Read mapped NN results  ##############################################################
##########################################################################################

Patchseq_NN_memb <- read.csv(paste0(work.dir, "/NN_10000epochs_batch500_patchseq_membership1.csv"), check.names=FALSE, row.names = 1)
FACs_NN_memb <- read.csv(paste0(work.dir,"/NN_10000epochs_batch500_FACS_membership1.csv") , check.names=FALSE, row.names = 1)
#Patchseq_NN_memb <- Patchseq_NN_memb[locked_cells_sample_id,]
#save(Patchseq_NN_memb, file=file.path(paste0("Final_matrices_locked_data/Patchseq_NN_10_10000ephocs_500batch_mapping.memb",".rda")))

#NN_FACs_memb = NN_FACs_memb / rowSums(FACs_NN_memb)
#NN_memb = NN_memb / rowSums(Patchseq_NN_memb)
colnames(Patchseq_NN_memb) <- select.cl
colnames(FACs_NN_memb) <- select.cl

##########################################################################################
### Read mapped Seurat results  ##########################################################
##########################################################################################

#source(paste0(Where,"/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R"))

#Seurat_memb = read.csv(paste0("Seurat_memb.csv"), check.names=FALSE, row.names = 1)
#colnames(Seurat_memb) <- gsub("\\."," " ,gsub("prediction.score.","",colnames(Seurat_memb)))
#colnames(Seurat_memb) <- gsub("L2 3", "L2/3", colnames(Seurat_memb))
#Seurat_memb <- Seurat_memb[,select.cl]


##########################################################################################
### Building mapping probability matrix ##################################################
##########################################################################################

set.seed(1983)
#Tree_mapping_probability = compute_mapping_probability(memb = FACS_Tree_memb, select.cells = FACs.cells, 
#                                                       select.cl = select.cl, ref.cl= cl)

#save(Tree_mapping_probability, file=file.path(paste0("Final_matrices_locked_data/REF_Tree_mapping_probability.rda")))
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/REF_mapping_probability.rda")

NN_mapping_probability = compute_mapping_probability(memb = FACs_NN_memb, select.cells = FACs.cells, 
                                                     select.cl = select.cl, ref.cl = cl)
#save(NN_mapping_probability, file=file.path(paste0("Final_matrices_locked_data/REF_NN_mapping_probability.rda")))
#load("Final_matrices_locked_data/AIT2.3.2/REF_mapping_probability.rda")


###########################################################################################
### KLdiv #################################################################################
###########################################################################################
set.seed(1983)
#Tree_FACs_KLdiv = compute_KLdiv(select.cl = select.cl, 
#                           select.cells = FACs.cells, 
#                           mapping_probability = Tree_mapping_probability, 
#                           memb = FACS_Tree_memb)
#save(Tree_FACs_KLdiv, file=file.path(paste0("Final_matrices_locked_data/Tree_FACs_KLdiv", ".rda")))
#load("Final_matrices_locked_data/AIT2.3.1/KLdiv.rda")

Patchseq_Tree_KLdiv = compute_KLdiv(select.cl = select.cl, 
                                    select.cells = locked_cells_sample_id, 
                                    #select.cells = colnames(query.dat.norm),
                                    mapping_probability = Tree_mapping_probability, 
                                    memb = Patchseq_Tree_memb)

#save(Patchseq_Tree_KLdiv, file=file.path(paste0("Final_matrices/Tree_Patchseq_KLdiv", ".rda")))
#save(Patchseq_Tree_KLdiv, file=file.path(paste0("Final_matrices_locked_data/Tree_Patchseq_KLdiv", ".rda")))
#load("Final_matrices_locked_data/AIM1.1/KLdiv.rda")

NN_FACs_KLdiv = compute_KLdiv(select.cl = select.cl, 
                              select.cells = FACs.cells, 
                              mapping_probability = NN_mapping_probability, 
                              memb = as.matrix(FACs_NN_memb))
#save(NN_FACs_KLdiv, file=file.path(paste0("Final_matrices_locked_data/NN_FACs_KLdiv", ".rda")))
#load("Final_matrices_locked_data/AIT2.3.2/KLdiv.rda")

Patchseq_NN_KLdiv = compute_KLdiv(select.cl = select.cl, 
                                  select.cells = locked_cells_sample_id, 
                                  #select.cells = colnames(query.dat.norm), 
                                  mapping_probability = NN_mapping_probability, 
                                  memb = as.matrix(Patchseq_NN_memb))

#save(Patchseq_NN_KLdiv, file=file.path(paste0("Final_matrices/NN_Patchseq_KLdiv", ".rda")))
#load("Final_matrices_locked_data/AIM1.2/KLdiv.rda")

##########################################################################################
### Correlation ##########################################################################
##########################################################################################

FACs_anno <- FACs_anno[FACs.cells, ] #Ask zizhen why some cells are not in norm.dat
FACs.cl.med <- Compute_median_gene_expression(anno_file = FACs_anno, norm.dat = norm.dat , markers = select.markers)
#save(FACs.cl.med, file=file.path(paste0("Final_matrices/FACs.cl.med", ".rda")))
#load("Final_matrices/FACs.cl.med.rda")

Patchseq_FACs_cor <- Compute_correlation_mat(markers = select.markers, cells = colnames(query.dat.norm), #locked_cells_sample_id, 
                                             query.dat.norm = query.dat.norm, train.cl.med = FACs.cl.med)

FACs_FACs_cor <- Compute_correlation_mat(markers = select.markers, cells = FACs.cells, 
                                         query.dat.norm = norm.dat, train.cl.med = FACs.cl.med)

###########################################################################################
### Aggreagte results #####################################################################
###########################################################################################

source(paste0(Where,"/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R"))

#Patchseq_Tree_memb <- Patchseq_Tree_memb[colnames(query.dat.norm),select.cl]
#NN_memb <- Patchseq_NN_memb[colnames(query.dat.norm), select.cl]
Patchseq_Tree_memb <- Patchseq_Tree_memb[locked_cells_sample_id,select.cl]
Patchseq_NN_memb <- Patchseq_NN_memb[locked_cells_sample_id, select.cl]
#Seurat_memb <- Seurat_memb[locked_cells_sample_id, select.cl]

Tree_3_cls <- Get_3_best_cl(Patchseq_Tree_memb, select.cl)
colnames(Tree_3_cls) <- c("Tree_first_cl", "Tree_second_cl", "Tree_third_cl")

Tree_3_bts <- Get_3_best_bt(Patchseq_Tree_memb, select.cl)
colnames(Tree_3_bts) <- c("Tree_first_bt", "Tree_second_bt", "Tree_third_bt")

NN_3_cls <- Get_3_best_cl(Patchseq_NN_memb, select.cl)
colnames(NN_3_cls) <- c("NN_first_cl", "NN_second_cl", "NN_third_cl")

NN_3_bts <- Get_3_best_bt(Patchseq_NN_memb, select.cl)
colnames(NN_3_bts) <- c("NN_first_bt", "NN_second_bt", "NN_third_bt")

#Seurat_3_cls <- Get_3_best_cl(Seurat_memb, select.cl)
#colnames(Seurat_3_cls) <- c("Seurat_first_cl", "Seurat_second_cl", "Seurat_third_cl")

#Seurat_3_bts <- Get_3_best_bt(Seurat_memb, select.cl)
#colnames(Seurat_3_bts) <- c("Seurat_first_bt", "Seurat_second_bt", "Seurat_third_bt")

Tree_3_KL <- Get_3_best_KL(memb = Patchseq_Tree_memb, ref.cl = select.cl, KLdiv = Patchseq_Tree_KLdiv)
colnames(Tree_3_KL) <-c("Tree_first_KL", "Tree_second_KL", "Tree_third_KL")

NN_3_KL <- Get_3_best_KL(memb = Patchseq_NN_memb, ref.cl = select.cl, KLdiv = Patchseq_NN_KLdiv)
colnames(NN_3_KL) <- c("NN_first_KL", "NN_second_KL", "NN_third_KL")

Tree_3_cor <- Get_3_best_cor(memb = Patchseq_Tree_memb, ref.cl = select.cl, cor = Patchseq_FACs_cor)
colnames(Tree_3_cor) <- c("Tree_first_cor", "Tree_second_cor", "Tree_third_cor")

NN_3_cor <- Get_3_best_cor(memb = Patchseq_NN_memb, ref.cl = select.cl, cor = Patchseq_FACs_cor)
colnames(NN_3_cor) <- c("NN_first_cor", "NN_second_cor", "NN_third_cor")

#cells <- colnames(query.dat.norm)
cells <- locked_cells_sample_id
results <- cbind.data.frame(Tree_3_cls[cells,],
                            NN_3_cls[cells,],
                            Tree_3_bts[cells,],
                            NN_3_bts[cells,],
                            Tree_3_KL[cells,],
                            NN_3_KL[cells,],
                            Tree_3_cor[cells,],
                            NN_3_cor[cells,])#,
#Seurat_3_cls[locked_cells_sample_id,],
#Seurat_3_bts[locked_cells_sample_id,])


##### FACs data

#Tree_3_cls <- Get_3_best_cl(FACS_Tree_memb, select.cl)
#colnames(Tree_3_cls) <- c("Tree_first_cl", "Tree_second_cl", "Tree_third_cl")

#Tree_3_bts <- Get_3_best_bt(FACS_Tree_memb, select.cl)
#colnames(Tree_3_bts) <- c("Tree_first_bt", "Tree_second_bt", "Tree_third_bt")

#NN_3_cls <- Get_3_best_cl(FACs_NN_memb, select.cl)
#colnames(NN_3_cls) <- c("NN_first_cl", "NN_second_cl", "NN_third_cl")

#NN_3_bts <- Get_3_best_bt(FACs_NN_memb, select.cl)
#colnames(NN_3_bts) <- c("NN_first_bt", "NN_second_bt", "NN_third_bt")

#Tree_3_KL <- Get_3_best_KL(memb = FACS_Tree_memb, ref.cl = select.cl, KLdiv = Tree_FACs_KLdiv)
#colnames(Tree_3_KL) <-c("Tree_first_KL", "Tree_second_KL", "Tree_third_KL")

#NN_3_KL <- Get_3_best_KL(memb = FACs_NN_memb, ref.cl = select.cl, KLdiv = NN_FACs_KLdiv)
#colnames(NN_3_KL) <- c("NN_first_KL", "NN_second_KL", "NN_third_KL")

#Tree_3_cor <- Get_3_best_cor(memb = FACS_Tree_memb, ref.cl = select.cl, cor = FACs_FACs_cor)
#colnames(Tree_3_cor) <- c("Tree_first_cor", "Tree_second_cor", "Tree_third_cor")

#NN_3_cor <- Get_3_best_cor(memb = FACs_NN_memb, ref.cl = select.cl, cor = FACs_FACs_cor)
#colnames(NN_3_cor) <- c("NN_first_cor", "NN_second_cor", "NN_third_cor")

#results <- cbind.data.frame(Tree_3_cls[FACs.cells,],
#                            NN_3_cls[FACs.cells,],
#                            Tree_3_bts[FACs.cells,],
#                            NN_3_bts[FACs.cells,],
#                            Tree_3_KL[FACs.cells,],
#                            NN_3_KL[FACs.cells,],
#                            Tree_3_cor[FACs.cells,],
#                            NN_3_cor[FACs.cells,])

##########################################################################################
### Assign cell identities: ##############################################################
##########################################################################################
cells <- locked_cells_sample_id
results <- cbind.data.frame(Tree_3_cls[cells,],
                            NN_3_cls[cells,],
                            Tree_3_bts[cells,],
                            NN_3_bts[cells,],
                            Tree_3_KL[cells,],
                            NN_3_KL[cells,],
                            Tree_3_cor[cells,],
                            NN_3_cor[cells,],
                            patchseq_anno[cells, c("topLeaf_id", "topLeaf_label", "topLeaf_color")])#,
Original_cols <- colnames(results)

results <- results %>% 
  rownames_to_column("id") %>%
  mutate(Tree_not_finall_call = ifelse(Tree_first_cor > 0.5  & Tree_first_KL < 2, "Good", "PoorQ")) %>%
  mutate(Tree_call = case_when(Tree_not_finall_call == "Good" & Tree_first_bt >= 0.9 ~ "Core",
                               Tree_not_finall_call == "Good" & Tree_first_bt < 0.9 &
                                 Tree_first_bt + Tree_second_bt >= 0.7 &
                                 Tree_first_bt / Tree_second_bt >= 2 ~ "I1", 
                               Tree_not_finall_call == "Good" & Tree_first_bt < 0.9 &
                                 Tree_first_bt + Tree_second_bt >= 0.7 &
                                 Tree_first_bt / Tree_second_bt < 2 ~ "I2",
                               Tree_not_finall_call == "Good" & Tree_first_bt < 0.9 &
                                 Tree_first_bt + Tree_second_bt < 0.7 ~ "I3",
                               Tree_not_finall_call == "PoorQ" ~ "PoorQ",
                               TRUE ~ "Other")) %>%
  mutate(NN_not_finall_call = ifelse(NN_first_cor > 0.5 & NN_first_KL < 2, "Good", "PoorQ")) %>%
  mutate(NN_call = case_when(NN_not_finall_call == "Good" & NN_first_bt >= 0.85 ~ "Core",
                             NN_not_finall_call == "Good" & NN_first_bt < 0.85 &
                               NN_first_bt + NN_second_bt >= 0.6 &
                               NN_first_bt / NN_second_bt >= 2 ~ "I1", 
                             NN_not_finall_call == "Good" & NN_first_bt < 0.85 &
                               NN_first_bt + NN_second_bt >= 0.6 &
                               NN_first_bt / NN_second_bt < 2 ~ "I2",
                             NN_not_finall_call == "Good" & NN_first_bt < 0.85 &
                               NN_first_bt + NN_second_bt < 0.6 ~ "I3",
                             NN_not_finall_call == "PoorQ" ~ "PoorQ",
                             TRUE ~ "Other")) %>%
  mutate(Tree_id = case_when(Tree_call == "Core" ~ 1,
                             Tree_call == "I1" ~ 2,
                             Tree_call == "I2" ~ 3,
                             Tree_call == "I3" ~ 4,
                             Tree_call == "PoorQ" ~ 5)) %>%
  mutate(Tree_color = case_when(Tree_call == "Core" ~ "Green",
                                Tree_call == "I1" ~ "Blue",
                                Tree_call == "I2" ~ "red",
                                Tree_call == "I3" ~ "Orange",
                                Tree_call == "PoorQ" ~ "Purple")) %>%
  mutate(NN_id = case_when(NN_call == "Core" ~ 1,
                           NN_call == "I1" ~ 2,
                           NN_call == "I2" ~ 3,
                           NN_call == "I3" ~ 4,
                           NN_call == "PoorQ" ~ 5)) %>%
  mutate(NN_color = case_when(NN_call == "Core" ~ "Green",
                              NN_call == "I1" ~ "Blue",
                              NN_call == "I2" ~ "red",
                              NN_call == "I3" ~ "Orange",
                              NN_call == "PoorQ" ~ "Purple")) %>%
  column_to_rownames("id") 

results <- results[,c(Original_cols, "Tree_call", "NN_call", "Tree_color", "NN_color", "Tree_id", "NN_id")]
sum(results$Tree_first_cl == results$NN_first_cl & results$Tree_call=="Core"  & results$NN_call=="Core" )
sum(results$Tree_first_cl == results$NN_first_cl & (results$Tree_call=="Core" | results$Tree_call=="I1") & (results$NN_call=="Core" | results$NN_call=="I1"))
sum(results$Tree_first_cl == results$NN_first_cl & (results$Tree_call=="Core" | results$Tree_call=="I1" | results$Tree_call=="I2") & 
      (results$NN_call=="Core" | results$NN_call=="I1" | results$NN_call=="I2"))
sum(results$Tree_first_cl == results$NN_first_cl & (results$Tree_call=="Core" | results$Tree_call=="I1" | results$Tree_call=="I2" | results$Tree_call=="I3") & 
      (results$NN_call=="Core" | results$NN_call=="I1" | results$NN_call=="I2" | results$NN_call=="I3"))

sum(results$Tree_call!="PoorQ"  & results$NN_call!="PoorQ" )

#sum(results$Tree_first_cl == results$NN_first_cl & results$Tree_call=="I2" & results$Tree_call=="I2")
#sum(results$Tree_first_cl == results$NN_first_cl & results$Tree_call=="I3" & results$Tree_call=="I3")
#sum(results$Tree_call!="PoorQ" & results$NN_call!="PoorQ")

sum(results$Tree_call!="PoorQ" & results$NN_call!="PoorQ")
#write.csv(results,"Final_matrices_locked_data/FACs_results.csv")
#write.csv(results,"Final_matrices_locked_data/Patchseq_results.csv")

#Tree_cls <- c("Tree_first_cl", "Tree_second_cl", "Tree_third_cl",
#              "Tree_first_bt", "Tree_second_bt", "Tree_third_bt",
#              "Tree_first_KL", "Tree_second_KL", "Tree_third_KL", 
#              "Tree_first_cor", "Tree_second_cor", "Tree_third_cor",
#              "Tree_call", "Tree_color", "Tree_id")
#NN_cls <- colnames(results)[!colnames(results) %in% Tree_cls]
#write.csv(results[,Tree_cls], "Final_matrices_locked_data/Tree_Patchseq_results.csv")
#write.csv(results[,NN_cls], "Final_matrices_locked_data/NN_Patchseq_results.csv")
#write.csv(results[,Tree_cls], "Final_matrices_locked_data/Tree_FACS_results.csv")
#write.csv(results[,NN_cls], "Final_matrices_locked_data/NN_FACS_results.csv")

results[locked_cells_sample_id,"Old_cluster"] <- Patchseq_Tree_mapping.df[locked_cells_sample_id, "cl"]
Leaf_node_cells <- rownames(Patchseq_Tree_mapping.df)[Patchseq_Tree_mapping.df$resolution.index==1]
Internal_node_cells <- rownames(Patchseq_Tree_mapping.df)[(Patchseq_Tree_mapping.df$resolution.index < 1 &Patchseq_Tree_mapping.df$resolution.index > 0.7)]
PoorQ_cells <- setdiff(rownames(Patchseq_Tree_mapping.df), c(Leaf_node_cells, Internal_node_cells))
Leaf_node_cells <- intersect(Leaf_node_cells, locked_cells_sample_id)
Internal_node_cells <- intersect(Internal_node_cells, locked_cells_sample_id)
PoorQ_cells <- intersect(PoorQ_cells, locked_cells_sample_id)

results[Leaf_node_cells, "Old_call"] <- c("Leaf_node")
results[Internal_node_cells, "Old_call"] <- c("Internal_node")
results[PoorQ_cells, "Old_call"] <- c("PoorQ")

results[Leaf_node_cells, "Old_color"] <- c("Green")
results[Internal_node_cells, "Old_color"] <- c("Blue")
results[PoorQ_cells, "Old_color"] <- c("Purple")

results[Leaf_node_cells, "Old_id"] <- c(1)
results[Internal_node_cells, "Old_id"] <- c(2)
results[PoorQ_cells, "Old_id"] <- c(3)

ggplot(results, aes(Tree_first_bt , fill = Tree_call)) + 
  geom_density(alpha = 0.3) + xlim(c(0,1)) + 
  xlab("first cluster call confidence") + ylab("Density")

ggplot(results, aes(Tree_second_bt  , fill = Tree_call)) + 
  geom_density(alpha = 0.3) + xlim(c(0,1)) +
  xlab("second cluster call confidence") + ylab("Density")

ggplot(results, aes(NN_first_bt , fill = NN_call)) + 
  geom_density(alpha = 0.3) + xlim(c(0,1)) + 
  xlab("first cluster call confidence") + ylab("Density")

ggplot(results, aes(NN_second_bt  , fill = NN_call)) + 
  geom_density(alpha = 0.3) + xlim(c(0,1)) +
  xlab("second cluster call confidence") + ylab("Density")

##########################################################################################
### Core, I1, I2 size ####################################################################
##########################################################################################

source(paste0(Where,"/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R"))

core_I1_size <- Get_size_of_core_I1_cells(results)
I2_confusions <- Get_size_of_I2_confusions(results)

write.csv(core_I1_size, "/Shared_with_IVSCC/Core_I1_size.csv")
write.csv(I2_confusions,"/Shared_with_IVSCC/I2_type_confusion.csv")

##########################################################################################
### River plot ###########################################################################
##########################################################################################
results <- results %>% rownames_to_column("sample_id")
ref_color_label_id = unique(FACs_anno[,c("cluster_id","cluster_label", "cluster_color")])

#ref_id <- data.frame(cluster_lable = select.cl, cluster_id = seq(93))
#ref_color_label_id <- left_join(ref_color, ref_id)
colnames(ref_color_label_id) <- c("Tree_first_cl_id", "Tree_first_cl", "Tree_firt_cl_color")
dim(results)
results <-  left_join(results, ref_color_label_id)
colnames(ref_color_label_id) <- c("NN_first_cl_id", "NN_first_cl", "NN_firt_cl_color")
results <- left_join(results, ref_color_label_id)
dim(results)
rownames(results) <- results$sample_id

GABAcells <- locked_cells_sample_id[patchseq_anno[locked_cells_sample_id, "subclass_label"] %in% c("Vip", "Sst", "Pvalb", "Lamp5", "Sncg")]
core_I1_GABAcells <- GABAcells[results[GABAcells, "Tree_call"] %in% c("Core", "I1") & results[GABAcells, "NN_call"] %in% c("Core", "I1")]
  
source(paste0( Where, "/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Plot_utils.R"))
library(dplyr)

tmp <- results %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id ,Tree_id, Tree_call, Tree_color, NN_id, NN_call, NN_color) %>%
  `colnames<-` (c("sample_id", "map_cluster_id", "map_cluster_label", 
                  "map_cluster_color", "cluster_id", "cluster_label", "cluster_color"))

river_plot(tmp, min.cells=0, min.frac=0)

tmp <- results[core_I1_GABAcells,] %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id , Tree_first_cl_id, Tree_first_cl, Tree_firt_cl_color, NN_first_cl_id, NN_first_cl, NN_firt_cl_color) %>%
  `colnames<-` (c("sample_id", "map_cluster_id", "map_cluster_label", 
                  "map_cluster_color", "cluster_id", "cluster_label", "cluster_color"))

river_plot(tmp, min.cells=0, min.frac=0)

tmp <- results %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id ,Tree_id, Tree_call, Tree_color, NN_id, NN_call, NN_color) %>%
  `colnames<-` (c("sample_id", "map_cluster_id", "map_cluster_label", 
                  "map_cluster_color", "cluster_id", "cluster_label", "cluster_color"))

Tree_cl_label <- set_names(results$Tree_first_cl, rownames(results))
Tree_cl_color <- as.character(Renew_list(ls = Tree_cl_label, ref.df = cl.df,label = "cluster_label",new.label = "cluster_color"))
Tree_cl_id <- as.character(Renew_list(ls = Tree_cl_label, ref.df = cl.df,label = "cluster_label",new.label = "cluster_id"))
NN_cl_label <- set_names(results$NN_first_cl, rownames(results))
NN_cl_color <- as.character(Renew_list(ls = NN_cl_label, ref.df = cl.df,label = "cluster_label",new.label = "cluster_color"))
NN_cl_id <- as.character(Renew_list(ls = NN_cl_label, ref.df = cl.df,label = "cluster_label",new.label = "cluster_id"))
results <- cbind.data.frame(results, Tree_cl_color, Tree_cl_id, NN_cl_color, NN_cl_id)
results$Tree_cl_color <- as.character(results$Tree_cl_color)
results$NN_cl_color <- as.character(results$NN_cl_color)
results$Tree_cl_id <- as.integer(results$Tree_cl_id)
results$NN_cl_id <- as.integer(results$NN_cl_id)
rm(NN_cl_id,NN_cl_color, NN_cl_label, Tree_cl_color, Tree_cl_id, Tree_cl_label)

INH_cl <- cl.df[cl.df$class_label == "GABAergic", "cluster_label"]
EXC_cl <-  cl.df[cl.df$class_label == "Glutamatergic", "cluster_label"]


source(paste0(Where,"/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Plot_utils.R"))

tmp <- results[results$Tree_call != "PoorQ" & results$NN_call!="PoorQ" & results$Tree_first_cl %in% INH_cl,] %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id ,Tree_cl_id, Tree_first_cl, Tree_cl_color, NN_cl_id, NN_first_cl, NN_cl_color) %>%
  `colnames<-` (c("sample_id", "map_cluster_id", "map_cluster_label", 
                  "map_cluster_color", "cluster_id", "cluster_label", "cluster_color"))
river_plot(tmp, min.cells=0, min.frac=0)

go_63x_cells = patchseq_anno[patchseq_anno$go_no_go_63x_label == "63x go", "sample_id"]
go_63x_cells = go_63x_cells[go_63x_cells %in% locked_cells_sample_id]

sub <- results[go_63x_cells,]

tmp <- sub[sub$Tree_call != "PoorQ" & sub$NN_call!="PoorQ" & sub$Tree_first_cl %in% INH_cl,] %>% 
  rownames_to_column("id") %>% 
  dplyr::select(id ,Tree_cl_id, Tree_first_cl, Tree_cl_color, NN_cl_id, NN_first_cl, NN_cl_color) %>%
  `colnames<-` (c("sample_id", "map_cluster_id", "map_cluster_label", 
                  "map_cluster_color", "cluster_id", "cluster_label", "cluster_color"))
river_plot(tmp, min.cells=0, min.frac=0)


##########################################################################################
############################################### PLOTS ####################################
##########################################################################################

ggplot(data = melt(Tree_mapping_probability[select.cl, select.cl]), aes(x=Var1, y=Var2, fill=value)) +
  geom_tile()+ theme(axis.text = element_text(size=5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("clustering_cluster_label") + ylab("NN_mapping_cluster_label") +  
  scale_fill_gradient(low = "white", high = "red")


NN = Get_nth_predicted_bt(Patchseq_NN_memb, select.cl, 1)
Tree = Get_nth_predicted_bt(Patchseq_Tree_memb, select.cl, 1)
plot_geomdensity_categorized(Tree,  NN, category_names = c("Tree", "NN"), xlabel = "Confidence")

nth_predicted_cl = Get_nth_predicted_cl(memb = FACS_Tree_memb, ref.cl = select.cl, nth_cl = 1)
Tree = Get_matrix_member(nth_predicted_cl, Tree_FACs_KLdiv)
nth_predicted_cl = Get_nth_predicted_cl(memb = FACs_NN_memb, ref.cl = select.cl, nth_cl = 1)
NN = Get_matrix_member(nth_predicted_cl, NN_FACs_KLdiv)
plot_geomdensity_categorized(Tree,  NN, category_names = c("Tree", "NN"), xlabel = "Confidence")

nth_predicted_cl = Get_nth_predicted_cl(memb = Patchseq_Tree_memb, ref.cl = select.cl, nth_cl = 1)
Tree = Get_matrix_member(nth_predicted_cl, Patchseq_Tree_KLdiv)
nth_predicted_cl = Get_nth_predicted_cl(memb = NN_memb, ref.cl = select.cl, nth_cl = 1)
NN = Get_matrix_member(nth_predicted_cl, Patchseq_NN_KLdiv)
plot_geomdensity_categorized(Tree,  NN, category_names = c("Tree", "NN"), xlabel = "Confidence")


nth_predicted_cl = Get_nth_predicted_cl(memb =FACS_Tree_memb, ref.cl = select.cl, nth_cl = 1)
FACs = as.data.frame(Get_matrix_member(nth_predicted_cl, FACs_FACs_cor))
FACs['Method'] <- c("FACs_FACs")
colnames(FACs) <- c("cor", "Method")
nth_predicted_cl = Get_nth_predicted_cl(memb = Patchseq_Tree_memb, ref.cl = select.cl, nth_cl = 1)
Tree = as.data.frame(Get_matrix_member(nth_predicted_cl, Patchseq_FACs_cor))
Tree['Method'] <- c("Tree_Patchseq_FACs")
colnames(Tree) <- c("cor", "Method")
nth_predicted_cl = Get_nth_predicted_cl(memb = NN_memb, ref.cl = select.cl, nth_cl = 1)
NN = as.data.frame(Get_matrix_member(nth_predicted_cl, Patchseq_FACs_cor))
NN['Method'] <- c("NN_Patchseq_FACs")
colnames(NN) <- c("cor", "Method")
tmp <- rbind(NN, Tree,FACs)
ggplot(tmp, aes(cor, fill = Method)) + 
  geom_density(alpha = 0.5) +
  xlab("Correlation with the best cluster") + ylab("Density")



######################################################################################################
######################################### Playing with data ##########################################
######################################################################################################

#Are there any duplicate?
distinct_query <- as.data.frame(query.dat.norm) %>% distinct()
dim(distinct_query)
distinct_FACS <- as.data.frame(norm.dat) %>% distinct()
dim(distinct_FACS)

#Is the data imbalanced?
plot(sort(table(cl)))

#Are some local features enough?
temp <- norm.dat[select.markers,]
plot(temp[,1000], temp[,2000])

################################################################################################
################################ Gaussian mixture model ########################################
################################################################################################

library(mclust)
BIC <- mclustBIC(results$Tree_primary_bt)
plot(BIC)
summary(BIC)
mod1 <- Mclust(results$Tree_primary_bt, x = BIC)
summary(mod1, parameters = TRUE)
plot(mod1, what = "classification")

mod4 <- densityMclust(results$NN_primary_bt, G=10)
summary(mod4)
#plot(mod4, what = "BIC")
plot(mod4, what = "density", data = results$NN_primary_bt, breaks = 40)


################################################################################################
################################ Ephys information #############################################
#Exc_highernodes = c("n5", "n14", "n6", "n15", "n39", "n40", "n31", "n30", "n28", "n29", "n33", "n34", "n35",
#                    "n7", "n9", "n10", "n11", "n41", "n45", "n43", "n46", "n47", "n48", "n50", "n51", "n53",
#                    "n54", "n55", "n56")
good.cells <- intersect(Tree.good.cells, NN.good.cells)
# We have the class and subclass information based on T_type for all the cells
#Tclass info
Exc_cls <- unique(FACs_anno[FACs_anno$class_label == "Glutamatergic", "cluster_label"])
Inh_cls <- unique(FACs_anno[FACs_anno$class_label == "GABAergic", "cluster_label"])
results <- results %>% 
  rownames_to_column("sample_id") %>% 
  mutate(Tree_primary_class = ifelse( Tree_primary_cl %in% Exc_cls , "Exc", "Inh")) %>%
  column_to_rownames("sample_id")
results <- results %>% 
  rownames_to_column("sample_id") %>%
  mutate(NN_primary_class = ifelse( NN_primary_cl %in% Exc_cls , "Exc", "Inh")) %>%
  column_to_rownames("sample_id")
results <- results %>% 
  rownames_to_column("sample_id") %>%
  mutate(Seurat_class = ifelse( Seurat_cl %in% Exc_cls , "Exc", "Inh")) %>%
  column_to_rownames("sample_id")

#Tsubclass info
pvalb_cls <- unique(FACs_anno[FACs_anno$subclass_label == "Pvalb", "cluster_label"])
Sst_cls <- unique(FACs_anno[FACs_anno$subclass_label == "Sst", "cluster_label"])
Vip_cls <- unique(FACs_anno[FACs_anno$subclass_label == "Vip", "cluster_label"])
Lamp5_cls <- unique(FACs_anno[FACs_anno$subclass_label == "Lamp5", "cluster_label"])
Mies2_cls <- unique(FACs_anno[FACs_anno$subclass_label == "Meis2", "cluster_label"])
Sncg_cls <- unique(FACs_anno[FACs_anno$subclass_label == "Sncg", "cluster_label"])
Serpinf1_cls <- unique(FACs_anno[FACs_anno$subclass_label == "Serpinf1", "cluster_label"])

results <- results %>%  
  rownames_to_column("sample_id") %>% 
  mutate(Tree_primary_subclass = case_when(
    Tree_primary_cl %in% Exc_cls ~ "Exc",
    Tree_primary_cl %in% pvalb_cls ~ "Pvalb",
    Tree_primary_cl %in% Sst_cls ~ "Sst",
    Tree_primary_cl %in% Vip_cls ~ "Vip",
    Tree_primary_cl %in% Lamp5_cls ~ "Lamp5",
    TRUE ~ "Other")) %>%
  mutate(NN_primary_subclass = case_when(
    NN_primary_cl %in% Exc_cls ~ "Exc",
    NN_primary_cl %in% pvalb_cls ~ "Pvalb",
    NN_primary_cl %in% Sst_cls ~ "Sst",
    NN_primary_cl %in% Vip_cls ~ "Vip",
    NN_primary_cl %in% Lamp5_cls ~ "Lamp5",
    TRUE ~ "Other")) %>%
  mutate(Seurat_subclass = case_when(
    Seurat_cl %in% Exc_cls ~ "Exc",
    Seurat_cl %in% pvalb_cls ~ "Pvalb",
    Seurat_cl %in% Sst_cls ~ "Sst",
    Seurat_cl %in% Vip_cls ~ "Vip",
    Seurat_cl %in% Lamp5_cls ~ "Lamp5",
    TRUE ~ "Other")) %>%
  column_to_rownames("sample_id")

#For the Eclass and subclass, we have a subset of cells available

E_type <- read.csv("E_type_predictions.csv")
E_subclass_type <- read.csv("MET_subclass_predictions.csv")

MET_sample_id <- sapply(E_type$index, function(x) patchseq_anno[patchseq_anno$spec_id_label==x, "sample_id"])
rownames(E_type) <- MET_sample_id

MET_sample_id <- sapply(E_subclass_type$index, function(x) patchseq_anno[patchseq_anno$spec_id_label==x, "sample_id"])
rownames(E_subclass_type) <- MET_sample_id

E_type = E_type[MET_sample_id,]
E_subclass_type = E_subclass_type[MET_sample_id,]

#Assinging all the E class and subclass info
tmp <- cbind(results[MET_sample_id, ], E_type[MET_sample_id,], E_subclass_type[MET_sample_id, ]) %>% 
  rownames_to_column("sample_id") %>%
  mutate(E_class = case_when(
    Exc > 0.7 ~ "Exc",
    Inh > 0.7 ~ "Inh",
    TRUE ~ "Other")) %>%
  mutate(E_subclass = case_when(
    Excitatory > 0.7 ~ "Exc",
    Lamp5 > 0.7 ~ "Lamp5",
    Vip > 0.7 ~ "Vip",
    Sst > 0.7 ~ "Sst",
    Pvalb > 0.7 ~ "Pvalb",
    TRUE ~ "Other")) %>%
  column_to_rownames("sample_id") %>%
  select(c("E_class", "E_subclass"))

results[,"E_class"] <- NA
results[,"E_subclass"] <- NA
results[MET_sample_id, "E_class"] <- tmp[MET_sample_id,"E_class"]
results[MET_sample_id, "E_subclass"] <- tmp[MET_sample_id,"E_subclass"]

results <- results %>% 
  rownames_to_column("id") %>%
  mutate(NN_call= ifelse(id %in% NN.good.cells, "Good", "Bad")) %>%
  mutate(Tree_call= ifelse(id %in% Tree.good.cells, "Good", "Bad")) %>%
  mutate(Old_tree_call = case_when(
    id %in% Leaf_node_cells ~ "Leaf_node",
    id %in% Internal_node_cells ~ "Internal_node",
    id %in% PoorQ_cells ~ "PoorQ",
    TRUE ~ "Other")) %>%
  column_to_rownames("id")

intersect(MET_sample_id, tmp2)
tmp <- rownames(results %>% rownames_to_column("sample_id") %>% filter(E_class=="Exc") %>% column_to_rownames("sample_id"))
#tmp2 <- setdiff(union(Tree.good.cells, NN.good.cells), intersect(Tree.good.cells, NN.good.cells))
tmp2 <- intersect(Tree.good.cells, NN.good.cells)
study.cells <- intersect(tmp2, tmp)
tmp2 <- intersect(intersect(NN.good.cells, Tree.good.cells), tmp)
study.cells <- tmp2

NN <- as.data.frame(table(results[study.cells, "NN_primary_class"]))
Tree <- as.data.frame(table(results[study.cells, "Tree_primary_class"]))
Seurat <- as.data.frame(table(results[study.cells, "Seurat_class"]))
tmp <- merge(merge(NN, Tree,  by="Var1", all.x = TRUE), Seurat, by = "Var1", all.x = TRUE)
colnames(tmp) <- c("Var1", "NN", "Tree", "Seurat")
tmp <- melt(tmp)
colnames(tmp) <- c("Variable", "Method", "Freq")

ggplot(tmp,aes(x=Variable,y=Freq,fill=Method))+
  geom_bar(stat="identity",position="dodge") +
  xlab("Class")+ylab("Class size") 


NN <- as.data.frame(table(results[study.cells, "NN_primary_subclass"])) %>% 
  filter(Var1 %in% c("Exc_subclass", "Vip", "Pvalb", "Sst", "Lamp5"))
Tree <- as.data.frame(table(results[study.cells, "Tree_primary_subclass"])) %>% 
  filter(Var1 %in% c("Exc_subclass", "Vip", "Pvalb", "Sst", "Lamp5"))
Seurat <- as.data.frame(table(results[study.cells, "Seurat_subclass"])) %>% 
  filter(Var1 %in% c("Exc_subclass", "Vip", "Pvalb", "Sst", "Lamp5"))
Ephys <- as.data.frame(table(results[study.cells, "E_subclass"]))

tmp <- merge(merge(merge(NN, Tree, by="Var1"), Seurat, by = "Var1", all.x = TRUE), Ephys, by = "Var1", all.x = TRUE)
colnames(tmp) <- c("Var1", "NN", "Tree", "Seurat", "Ephys")
tmp <- melt(tmp)
colnames(tmp) <- c("Variable", "Method", "Freq")
ggplot(tmp,aes(x=Variable,y=Freq,fill=Method))+
  geom_bar(stat="identity",position="dodge") +
  xlab("Subclass")+ylab("Subclass size") 


################################################################################################
################################ Cre information ###############################################

Build_cre_info_table <- function(ref.anno, ref.cells, ref.cl, cl) {
  
  cre_info <- ref.anno[ref.cells,] %>% 
    as.data.frame()  %>%  
    group_by(genotype_label, cluster_label) %>% 
    summarise(Freq = n())
  
  ref_cluster_size <- table(cl[ref.cells]) / sum(table(cl))
  
  cre_info <- dcast(as.data.frame(cre_info),formula = genotype_label ~ cluster_label, drop = FALSE, value.var = "Freq")
  rownames(cre_info) <- cre_info$genotype_label
  cre_info <- cre_info[, ref.cl]
  cre_info <- cre_info/rowSums(cre_info, na.rm = TRUE)
  cre_info <- log2(cre_info[,ref.cl]/ ref_cluster_size[ref.cl])
  cre_info
}

Get_genotype_ref_probability <- function(ref.cre, cells, predicted_cl, genotype){
  genotype <- genotype[cells]
  predicted_cl <- predicted_cl[cells]
  cre_info <- sapply(1:length(cells), function(x)ref.cre[genotypes[x], predicted_cl[x]])
  names(cre_info) <- cells
  cre_info
}


FACs_cre_info <- Build_cre_info_table(FACs_anno, FACs.cells, select.cl, cl)

ggplot(data = melt(as.matrix(FACs_cre_info[tt,])), aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Cluster labels") + ylab("Genotype label") +  
  scale_fill_gradient(low = "white", high = "red")


genotypes = setNames(patchseq_anno$genotype_label, rownames(patchseq_anno))

cells <- query.dat.cells
cells <-union(Tree.good.cells, NN.good.cells)
predicted_cl = setNames(results[cells, "Tree_primary_cl"], cells)
Tree <- Get_genotype_ref_probability(FACs_cre_info, cells, predicted_cl, genotypes)

predicted_cl = set_names(results[cells, "NN_primary_cl"], cells)
NN <- Get_genotype_ref_probability(FACs_cre_info, cells, predicted_cl, genotypes)

predicted_cl = set_names(results[cells, "Seurat_cl"], cells)
Seurat <- Get_genotype_ref_probability(FACs_cre_info, cells, predicted_cl, genotypes)

df <- cbind.data.frame(Tree, NN, Seurat)
df[is.na(df)] <- -1000.
Method <- setNames(colnames(df)[apply(df,1,which.max)], rownames(df))
Best_value <- apply(df, 1, max)
df <- cbind(df, Method, Best_value)

df<- df %>% rownames_to_column("id") %>%
  mutate(Method = case_when(
    Best_value == -1000 ~ "missing",
    Best_value != -1000 & Tree == NN  & NN == Seurat ~ "common",
    Best_value != -1000 & Tree == NN  & Tree == Best_value ~ "common_tree_NN",
    Best_value != -1000 & Tree == Seurat  & Tree == Best_value ~ "common_tree_Seurat",
    Best_value != -1000 & Seurat == NN  & Seurat == Best_value ~ "common_NN_Seurat",
    #Best_value != -1000 & Seurat != NN & Seurat != Tree & Tree != NN ~ colnames(df)[apply(df,1,which.max)],
    TRUE ~ colnames(df)[apply(df,1,which.max)])) %>%
  column_to_rownames("id")

tmp <- table(df$Method)

tmp<- unlist(list( tmp["missing"],tmp["common"],
                   tmp["NN"] + tmp["common_tree_NN"]+ tmp["common_NN_Seurat"],
                   tmp["Seurat"] + tmp["common_NN_Seurat"]+tmp["common_tree_Seurat"],
                   tmp["Tree"]+tmp["common_tree_NN"]+ tmp["common_tree_Seurat"]))

barplot(tmp)

#Gave to Nathan
tmp <- rownames(results[results$Tree_primary_KL <1,])
df <- patchseq_anno[tmp, c("sample_id", "spec_id_label")]
df$spec_id_label <- as.character(df$spec_id_label)
write.csv(df, "Informative_patchseq_cells.csv")


################################################################################################
################################ Thresholding  #################################################



#This must HOLD
sum(as.character(Tree_mapping.df[Leaf_node_cells, "cl"]) == results[Leaf_node_cells, "Tree_primary_cl"]) == length(Leaf_node_cells)


dim(results)
results$Tree_primary_Zscore <- (results$Tree_primary_cor - mean(results$Tree_primary_cor)) / sd(results$Tree_primary_cor)
results$NN_primary_Zscore <- (results$NN_primary_cor - mean(results$NN_primary_cor)) / sd(results$NN_primary_cor)
results$Tree_likelihood = log(results$Tree_secondary_bt + 0.0000001/(results$Tree_primary_bt))
results$NN_likelihood = log(results$NN_secondary_bt + 0.0000001/(results$NN_primary_bt))
hist(results$Tree_likelihood)
hist(results$NN_likelihood)

Tree.core.cells
Tree.core.cells <- rownames(results)[(results$Tree_primary_KL < 2.5 & (results$Tree_primary_bt >0.7 ) & results$Tree_primary_Zscore > -2)]
#Tree.core.cells <- rownames(results)[(results$Tree_primary_KL < 2 & results$Tree_primary_Zscore > -2)]
#Tree.core.cells <- rownames(results)[(results$Tree_primary_KL < 2 & (results$Tree_primary_bt + results$Tree_secondary_bt > 0.9) & results$Tree_primary_Zscore > -2)]
Tree.intermediate.cells <- rownames(results)[results$Tree_primary_KL < 2]
Tree.intermediate.cells <- setdiff(Tree.intermediate.cells, Tree.core.cells)
Tree.poorQ.cells <- setdiff(rownames(results), c(Tree.core.cells, Tree.intermediate.cells))
#This must HOLD
(length(Tree.intermediate.cells) + length(Tree.poorQ.cells) + length(Tree.core.cells)) == length(rownames(results))
hist(results[Tree.core.cells, "Tree_primary_bt"])
hist(results[Tree.core.cells, "Tree_primary_KL"])
hist(results[Tree.core.cells, "Tree_primary_cor"])
sum(Leaf_node_cells %in% Tree.core.cells)


NN.core.cells <- rownames(results)[(results$NN_primary_KL < 2.5  & results$NN_primary_Zscore > -2 & results$NN_primary_bt > 0.7)]
NN.intermediate.cells <- rownames(results)[results$NN_primary_KL < 2]
NN.intermediate.cells <- setdiff(NN.intermediate.cells, NN.core.cells)
study.cells <- setdiff(NN.core.cells, Tree.core.cells)
hist(results[NN.core.cells, "Tree_primary_bt"])
hist(results[NN.core.cells, "NN_primary_bt"])
sum(NN.core.cells %in% Tree.core.cells) / length(NN.core.cells)
sum(Tree.core.cells %in% NN.core.cells) / length(Tree.core.cells)
length(Tree.core.cells) 
length(NN.core.cells) 


hist(results[study.cells, "Tree_primary_bt"])
hist(results[study.cells, "NN_primary_bt"])

#NN.intermediate.cells <- rownames(results)[results$Tree_primary_KL < 2]
#NN.intermediate.cells <- setdiff(Tree.intermediate.cells, Tree.core.cells)
#NN.poorQ.cells <- setdiff(rownames(results), c(Tree.core.cells, Tree.intermediate.cells))
#This must HOLD
#(length(Tree.intermediate.cells) + length(Tree.poorQ.cells) + length(Tree.core.cells)) == length(rownames(results))




sum(results[Tree.core.cells, "Tree_primary_cl"] == results[Tree.core.cells, "NN_primary_cl"]) / length(Tree.core.cells)
max(results[results$Tree_primary_bt>0.7, "Tree_primary_KL"])



plot_cl_cor_cat <- function(core.cells, intermediate.cells, lowQ.cells, cell_cat){
  
  cormat1= return_cormat(markers = select.markers, cells = core.cells, 
                         query.dat.norm = query.dat.norm, train.cl.med = train.cl.med)
  best.cl= return_corcl(cormat1)
  best.cl.score1 = return_corescore(cormat1)
  
  
  cormat2= return_cormat(markers = select.markers, cells = intermediate.cells, 
                         query.dat.norm = query.dat.norm, train.cl.med = train.cl.med)
  best.cl2= return_corcl(cormat2)
  best.cl.score2 = return_corescore(cormat2)
  
  cormat3= return_cormat(markers = select.markers, cells = lowQ.cells, 
                         query.dat.norm = query.dat.norm, train.cl.med = train.cl.med)
  best.cl3= return_corcl(cormat3)
  best.cl.score3 = return_corescore(cormat3)
  
  tmp1 <- best.cl.score1 %>% as.data.frame() %>% rownames_to_column("sample_id") %>% mutate(cat=cell_cat[1]) 
  tmp2 <- best.cl.score2 %>% as.data.frame() %>% rownames_to_column("sample_id") %>% mutate(cat=cell_cat[2]) 
  tmp3 <- best.cl.score3 %>% as.data.frame() %>% rownames_to_column("sample_id") %>% mutate(cat=cell_cat[3]) 
  
  tmp <- rbind(tmp2,tmp1, tmp3)
  colnames(tmp) <- c("sample_id","best.cl.score", "cat")
  
  ggplot(tmp, aes(x=best.cl.score, fill = cat)) +
    geom_histogram(alpha = 0.5, position = 'identity', binwidth = 0.02) +
    scale_x_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1, 1.1)) 
  
}



hist(results$NN_secondary_bt/(results$NN_primary_bt + 0.0000001) , breaks = 50)
hist(results$Tree_secondary_bt/(results$Tree_primary_bt + 0.0000001) , breaks = 50)

#tmp1 <- rownames(NN_sum_stat)[NN_sum_stat$primary_bt>.7]
#tmp1 <- rownames(NN_sum_stat)[log(NN_sum_stat$primary_bt/(NN_sum_stat$seconday_bt + 0.0000001)) >0.4]
tmp1 <- rownames(NN_sum_stat)[(NN_sum_stat$primary_bt + NN_sum_stat$seconday_bt) > 0.9]
tmp1<- tmp1[NN_sum_stat[tmp1,"primary_KL"]<2.5]
tmp1 <- tmp1[NN_sum_stat[tmp1,"primary_cor"]>0.6]

tmp2 <- rownames(Tree_sum_stat)[(Tree_sum_stat$primary_bt + Tree_sum_stat$seconday_bt) > 0.9]
tmp2 <- tmp2[Tree_sum_stat[tmp2,"primary_KL"]<2.5]
tmp2 <- tmp2[Tree_sum_stat[tmp2,"primary_cor"]>0.6]

sum(tmp2 %in% tmp1) / length(tmp2)
sum(tmp1 %in% tmp2) / length(tmp1)
length(tmp1)
length(tmp2)


gr1 <- rownames(NN_sum_stat[(NN_sum_stat$primary_bt >= 0.7 ), ])
gr1 <- gr1[NN_sum_stat[gr1, "primary_KL"] <2]
gr2 <- rownames(NN_sum_stat[(NN_sum_stat$primary_KL < 2 ), ])
gr2 <- setdiff(gr2, gr1)

plot(NN_sum_stat[gr1, "primary_cor"])
plot(NN_sum_stat[gr2, "primary_cor"])

hist(Tree_sum_stat[study.cells1, "primary_cor"])
hist(Tree_sum_stat[study.cells1, "primary_KL"])

study.cells1 <- rownames(NN_sum_stat[(NN_sum_stat$primary_bt > 0.7 ), ])
hist(NN_sum_stat[study.cells1, "primary_KL"])
hist(NN_sum_stat[study.cells1, "primary_cor"])

plot( Tree_sum_stat$primary_KL, Tree_sum_stat$primary_bt)
plot( NN_sum_stat$primary_KL, NN_sum_stat$primary_bt)

plot( Tree_sum_stat$primary_cor, Tree_sum_stat$primary_bt)
plot( NN_sum_stat$primary_cor, NN_sum_stat$primary_bt)

plot( Tree_sum_stat$primary_KL, Tree_sum_stat$primary_cor)
plot( NN_sum_stat$primary_KL, NN_sum_stat$primary_cor)

Leaf_node_cells <- rownames(Tree_mapping.df)[Tree_mapping.df$resolution.index==1]
Internal_node_cells <- rownames(Tree_mapping.df)[(Tree_mapping.df$resolution.index<1 & Tree_mapping.df$resolution.index>0.7)]
PoorQ_cells <- rownames(Tree_mapping.df)[(Tree_mapping.df$resolution.index<=0.7)]

hist(Tree_sum_stat[Leaf_node_cells, "primary_bt"])
hist(Tree_sum_stat[Internal_node_cells, "primary_bt"])
hist(Tree_sum_stat[PoorQ_cells, "primary_bt"])

plot_cl_cor_cat(core.cells = Leaf_node_cells, 
                intermediate.cells = Internal_node_cells, 
                lowQ.cells  = PoorQ_cells, 
                cell_cat = c("Leaf node cells", "Internal_node_cells", "PoorQ_cells"))


plot_cl_cor_cat(core.cells = Tree.Core.cells, 
                intermediate.cells = Tree.intermediate.cells, 
                lowQ.cells  = Tree.poorQ.cells, 
                cell_cat = c("Core_cells", "Intermediate_cells", "PoorQ_cells"))


plot_cl_cor_cat(core.cells = NN.Core.cells, 
                intermediate.cells = NN.intermediate.cells, 
                lowQ.cells  = NN.poorQ.cells, 
                cell_cat = c("Core_cells", "Intermediate_cells", "PoorQ_cells"))


study.cells <- Leaf_node_cells[!Leaf_node_cells %in% Tree.Core.cells]
hist(Tree_KL_bt[study.cells])


#Study only the Core and transition new cells to remove bad cells from there
study.cells <- c(Tree.Core.cells, Tree.intermediate.cells)
cormat= return_cormat(markers = select.markers, cells = study.cells, 
                      query.dat.norm = query.dat.norm, train.cl.med = train.cl.med)
best.cl= return_corcl(cormat)
best.cl.score = return_corescore(cormat)

df <- best.cl %>% as.data.frame() %>%
  mutate(best.cl.score)
colnames(df) <- c("cl_label", "best.cl.score")
rownames(df) <- names(best.cl)

best.Z = tapply(rownames(df), df$cl_label, function(x)(df[x,"best.cl.score"] - mean(df[x,"best.cl.score"]))/sd(df[x,"best.cl.score"]))
remove_cells <- unlist(sapply(best.Z, function(x)names(x)[x < -2]), use.names = FALSE)
new_Core_cells <- setdiff(Tree.Core.cells, remove_cells)
new_Transition_cells <- setdiff(Tree.intermediate.cells, remove_cells)
new_LowQ_cells <- c(Tree.poorQ.cells, remove_cells)
length(c(new_Core_cells, new_Transition_cells, new_LowQ_cells))

sum(patchseq_anno[select.cells, "cluster_label"] == mapping.df[select.cells, "cl"])
sst_nodes <- c("n105", "n91", "n108", "n95", "n96", "n94", "n107", "n99", "n98", "n100", "n101", "n90", select.cl[grep("Sst", select.cl)][-1])
old_sst_cells <- rownames(mapping.df)[as.character(mapping.df$cl) %in% sst_nodes]
sum(old_sst_cells %in% old_Transition_cells)
sum(old_sst_cells %in% old_Core_cells)

study.cells <- old_LowQ_cells[old_LowQ_cells%in% new_Transition_cells]















































#Align this data with FACs using Seurat
library(Seurat)
library("tibble")
library("tidyr")
library("dplyr")
library(ggplot2)
library(cowplot)
library(umap)
#install.packages("umap")


#Merge FACs and patchseq data
sum(rownames(norm.dat[select.markers,]) == rownames(query.dat.norm[select.markers,])) == length(select.markers)
all_data <- cbind(norm.dat[select.markers,], query.dat.norm[select.markers,])
metadata <- data.frame(row.names = c(colnames(norm.dat[select.markers,]), colnames(query.dat.norm[select.markers,])))
metadata[colnames(norm.dat[select.markers,]), "tech"] = "FACs"
metadata[colnames(query.dat.norm[select.markers,]), "tech"] = "Patchseq"

seuratobj <- CreateSeuratObject(counts = all_data, meta.data = metadata)
cell.list <- SplitObject(object = seuratobj, split.by = "tech")

for (i in 1:length(x = cell.list)) {
  #cell.list[[i]] <- NormalizeData(object = cell.list[[i]], verbose = FALSE)
  VariableFeatures(object = cell.list[[i]]) <- select.markers
}

reference.list <- cell.list[c("Patchseq", "FACs")]
anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, anchor.features = 4020)

#Old data
#cortex.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
#DefaultAssay(object = cortex.integrated) <- "integrated"
#features <- VariableFeatures(object = cortex.integrated) 
#DefaultAssay(object = cortex.integrated) <- "RNA"
#VariableFeatures(object = cortex.integrated)  <- features
#cortex.integrated <- ScaleData(object = cortex.integrated, verbose = FALSE)
#cortex.integrated <- RunPCA(object = cortex.integrated, npcs = 30, verbose = FALSE)
#cortex.integrated <- RunUMAP(object = cortex.integrated, reduction = "pca", dims = 1:30)
#p1 <- DimPlot(object = cortex.integrated, reduction = "umap", group.by = "tech")
#plot(p1)

#New intergrated data
cortex.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(object = cortex.integrated) <- "integrated"
cortex.integrated <- ScaleData(object = cortex.integrated, verbose = FALSE)
cortex.integrated <- RunPCA(object = cortex.integrated, npcs = 30, verbose = FALSE)
cortex.integrated <- RunUMAP(object = cortex.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(object = cortex.integrated, reduction = "umap", group.by = "tech")
plot(p1)

aligned_data = as.matrix(GetAssayData(cortex.integrated, assay="integrated"))
non_aligned_data = as.matrix(GetAssayData(cortex.integrated, assay="RNA"))
write.csv(aligned_data[,colnames(query.dat.norm)], file = "aligned_patchseq.csv")
write.csv(aligned_data[,colnames(norm.dat)], file = "aligned_FACs.csv")

#Preparing the M1 data for analysis
#/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/Miniatlas/transcriptome_joint/

############################## The whole thing in the block is just for reading the gene expression data ######################
###############################################################################################################################
setwd("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/Miniatlas/transcriptome_joint/")
library(matrixStats)
library(dplyr)
library(scrattch.hicat)
library(RANN)
library(Matrix)
library(ggplot2)
source('/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/scrattch.hicat/R/joint_analysis.R')
source('/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/scrattch.hicat/R/merge_cl.R')
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
blue.red <-colorRampPalette(c("blue", "white", "red"))
sets = grep("cells_|nuclei_", dir(), value=T)
save(sets, file="sets.rda")
cols = c("cluster_id", "cluster_label", "cluster_color", "class_label")
cl.result =  sapply(sets, function(x){
  print(x)
  load(file.path(x, "cl.final.rda"))
  cl.df$cluster_id = 1:nrow(cl.df)
  if(!is.factor(cl)){
    cl = as.factor(cl)
  }
  cl.clean = droplevels(cl[cl %in% row.names(cl.df)[!cl.df$class_label %in% c("Low Quality","Noise")]])
  cl.df = cl.df[levels(cl.clean), cols]
  return(list(cl=cl.clean, cl.df=cl.df))
},simplify=FALSE)
cl.list = sapply(cl.result, function(x)x[[1]], simplify=F)
cl.df.list = sapply(cl.result, function(x)x[[2]], simplify=F)
dat.list = sapply(sets, function(x){
  load(file.path( x, "norm.dat.rda"))
  norm.dat[, names(cl.list[[x]])]
  #dat.list[[x]][,names(cl.list[[x]])]
},simplify=FALSE)
ref.sets = c("10X_cells_AIBS","10X_cells_V3_AIBS", "10X_nuclei_V3_AIBS", "SmartSeq_cells_AIBS", "SmartSeq_nuclei_AIBS")
de.param.list = list(de_param(q1.th=0.4, q.diff.th=0.7, de.score.th=150, min.cells=15),
                     de_param(q1.th=0.4, q.diff.th=0.7, de.score.th=150, min.cells=10),
                     de_param(q1.th=0.4, q.diff.th=0.7, de.score.th=150, min.cells=10),
                     de_param(q1.th=0.5, q.diff.th=0.7, de.score.th=150, min.cells=4),
                     de_param(q1.th=0.4, q.diff.th=0.7, de.score.th=150, min.cells=4))
names(de.param.list) = ref.sets
comb.dat = prepare_joint(dat.list, de.param.list = de.param.list, cl.list = cl.list, cl.df.list = cl.df.list)
load("V4/comb.markers.rda")
SmartSeq_cells_AIBS.norm.dat = comb.dat$dat.list$SmartSeq_cells_AIBS[comb.markers,]
TenX_cells_AIBS.norm.dat = comb.dat$dat.list$`10X_cells_AIBS`[comb.markers,]


###############################################################################################################################
###############################################################################################################################

TenX_cells_AIBS = colnames(TenX_cells_AIBS.norm.dat)
SmartSeq_cells_AIBS = colnames(SmartSeq_cells_AIBS.norm.dat)

M1.cl = read.csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/Miniatlas/transcriptome_joint/V4/cl.csv")
M1.cl.df = read.csv("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/Miniatlas/transcriptome_joint/V4/cl.df.csv")

SmartSeq_cells_AIBS.labels <- M1.cl[M1.cl$X %in% SmartSeq_cells_AIBS,] %>% mutate(cluster_label = M1.cl.df[x,"cluster_label"])
colnames(SmartSeq_cells_AIBS.labels) <- c("sample_id", "factor_cl", "cl")

TenX_cells_AIBS.labels <- M1.cl[M1.cl$X %in% TenX_cells_AIBS,] %>% mutate(cluster_label = M1.cl.df[x,"cluster_label"])
colnames(TenX_cells_AIBS.labels) <- c("sample_id", "factor_cl", "cl")

TenX_cells = as.character(TenX_cells_AIBS.labels$sample_id)
SmartSeq_cells = as.character(SmartSeq_cells_AIBS.labels$sample_id)

path = paste0(mydir , "TenX_cells.csv")
write.csv(TenX_cells, file = path)

path = paste0(mydir , "SmartSeq_cells.csv")
write.csv(SmartSeq_cells, file = path)

path = paste0(mydir , "comb_markers.csv")
write.csv(comb.markers, file = path)

path = paste0(mydir , "SmartSeq_cells_AIBS.norm.csv")
writeMM(SmartSeq_cells_AIBS.norm.dat[comb.markers, SmartSeq_cells], file = path)

path = paste0(mydir , "TenX_cells_AIBS.norm.csv")
writeMM(TenX_cells_AIBS.norm.dat[comb.markers, TenX_cells], file = path)
TenX_cells_AIBS.norm.dat = readMM("TenX_cells_AIBS.norm.csv")

path = paste0(mydir , "SmartSeq_cells_AIBS_labels.csv")
write.csv(SmartSeq_cells_AIBS.labels, file = path)
path = paste0(mydir , "TenX_cells_AIBS_labels.csv")
write.csv(TenX_cells_AIBS.labels, file = path)

M1_select_cl <- M1.cl.df$cluster_label
path = paste0(mydir , "M1_select_cl.csv")
write.csv(M1_select_cl, file = path)

TenX_NN_results = read.csv(paste0(mydir, "TenX_NN_results.csv"))
predicted_NN_cl <- setNames(  factor(TenX_NN_results$predicted_cl+1, levels = length(M1_select_cl)), TenX_NN_results$sample_id)
ref_cl <- setNames(  factor(TenX_NN_results$cl, levels = M1_select_cl), TenX_NN_results$sample_id)

predicted_NN_cl <- setNames(factor(TenX_NN_results$predicted_cl+1, levels =seq(81)), TenX_NN_results$sample_id)
ref_cl <- setNames(factor(TenX_NN_results$factor_cl+1, levels =seq(81)), TenX_NN_results$sample_id)

confusion_matrix = TenX_NN_results[,c("cl", "NN_cl")] %>% table()
library(reshape2)
ggplot(data = melt(confusion_matrix), aes(x=cl, y=NN_cl, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Clustering lable") + ylab("NN Mapping lables") +  
  scale_fill_gradient(low = "white", high = "red")

rownames(TenX_NN_results) <- TenX_NN_results$X
rownames(TenX_cells_AIBS.labels) <- TenX_cells_AIBS.labels$sample_id

NN <- factor(TenX_NN_results[TenX_cells, "NN_cl"], levels = M1_select_cl)
Clustering <- factor(TenX_cells_AIBS.labels[TenX_cells, "cl"])
NN <- table(NN)
Clustering <- table(Clustering)
Clustering <- Clustering[M1_select_cl]
NN <- NN[M1_select_cl]
tmp <- melt(rbind(Clustering, NN))
colnames(tmp) <- c("Method", "cl", "value")

ggplot(tmp ,aes(x=cl,y=value,fill=Method))+
  geom_bar(stat="identity",position="dodge") +
  xlab("")+ylab("Cluster size") + theme(axis.text.x = element_text(angle = 90, hjust = 1))

rownames(TenX_NN_results)[TenX_NN_results$NN_cl != TenX_NN_results$cl]
mat <- build_confusion_matrix(mat1 = TenX_NN_results, mat2 = TenX_NN_results, 
                              colname1 = "NN_cl", colname2 = "cl", 
                              cells = TenX_cells)
plot_confusion_matrix(confusionmatrix = mat, cells = TenX_cells, ref.cl = M1_select_cl)
