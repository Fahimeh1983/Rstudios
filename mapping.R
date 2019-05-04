library(dendextend)
source("patchseq/heatmap.R")
source("patchseq/de.genes.R")
source("patchseq/dendro.R")
source("patchseq/patchseq.R")
library(matrixStats)
library(feather)

#just update batch_date and source it
#batch_date="20190227_BT014-RSC-195"

blue.red <-colorRampPalette(c("blue", "white", "red"))
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

######################################################################################################
### Setting up some paths ############################################################################
######################################################################################################

ref.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520"
res.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"
query.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"

######################################################################################################
### Reading ref data #################################################################################
######################################################################################################
#?????????????????????????? Where to put this
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/REF_mapping_probability.rda")

bp.collapse.th = 40
bp.name.add = NULL
if (!is.null(bp.collapse.th)) {
  bp.name.add = paste0(".with.bp.", bp.collapse.th)
}

###load reference data and tree
# make sure to load the files created by build_reference_dend
tmp.load1 = load(file=file.path(res.dir, "ref.data.rda")) # should include cl, cl.df, norm.dat. # The loaded cl is not used because it only includes cluster ids but not cluster labels 
tmp.load2 = load(file.path(file=res.dir, file=paste0("V1.dend", bp.name.add,".rda"))) # should include the pruned V1 tree
tmp.load3 = load(file.path(res.dir, file=paste0("V1.dend.list", bp.name.add,".rda"))) # should include dend.list

plot(dend)

#anno has ALM, so we can't use it
#ref.samp.dat = as.data.frame(read_feather(file.path(ref.dir, "anno.feather")))
#rownames(ref.samp.dat) <- ref.samp.dat$sample_id
#ref.samp.dat <- subset(ref.samp.dat, select = -c(sample_id))
#select.cl <- labels(dend)
#ref.samp.dat <- ref.samp.dat[ref.samp.dat$cluster_label %in% select.cl,]
#cl = setNames(as.character(ref.samp.dat$cluster_label), rownames(ref.samp.dat)) 
#cl = setNames(factor(cl, levels=unique(ref.samp.dat$cluster_label)), names(cl))

rownames(cl.df)=cl.df$cluster_id
cltmp=cl.df[as.character(cl),"cluster_label"]
names(cltmp)=names(cl)
cl=factor(cltmp)

######################################################################################################
### Loading the query data ###########################################################################
######################################################################################################

#### load query data (either load from R objects or from an already existed shiny folder or use the FACS data to map to itself)
#### option 1 : load Patchseq data

tmp<-load(paste0(query.dir,batch_date,"_mouse_patchseq_star2.0_cpm.Rdata"))
query.dat = cpmR

# loading samp.dat object
tmp<-load(paste0(query.dir,batch_date,"_mouse_patchseq_star2.0_samp.dat.Rdata"))

keepcells = which(samp.dat$Region=="VISp" & samp.dat$Type=="patch_seq")
samp.dat = samp.dat[c(keepcells, which(samp.dat$Region=="TCx"),which(samp.dat$Region=="FCx"),which(samp.dat$Region=="MOp"),which(samp.dat$Region=="TEa")   ),]   #FCx is for Brian.  Rat samples mapped in mouse

query.dat = query.dat[,as.character(samp.dat$exp_component_name)]
colnames(query.dat)=as.character(samp.dat$patched_cell_container)

query.dat.norm = log2(as.matrix(query.dat+1))
idx=match(rownames(norm.dat), rownames(query.dat.norm))
query.dat.norm=query.dat.norm[idx,]

######################################################################################################
### Some initialization ##############################################################################
######################################################################################################

set.seed(1983)
query.dat.cells = colnames(query.dat.norm)
FACS.cells <- colnames(norm.dat)

######################################################################################################
### Mapping data using Tree ##########################################################################
######################################################################################################

set.seed(1983)
# The key algorithm : Run 100 iterations, in each one sample 70% of the cells.
# low.th is the minimum differnce in Pearson correlation required to decide on which branch to map to. If the difference 
# is lower than this threshold, a random branch is chosen.
# The resulted memb object is a matrix where every row is a cell, and columns are the nodes of the dendrogram. The values are the bootstrap support for that node.
memb = map_dend_membership(dend, cl, norm.dat, query.dat.norm, query.dat.cells, bs.num=100, p=0.7, low.th=0.15)
# analyze the results to generate the final output, i.e. for every cell the deepest node with th=0.8 confidence.
# Also apply constraints on the minimum number of genes (or genes ratio)
#mapping.df = summarize_cl(dend, memb, query.dat.norm, conf.th=0.7, min.genes=1, min.genes.ratio=0.3)

######################################################################################################
### KLdiv ############################################################################################
######################################################################################################

set.seed(1983)
select.cl = labels(dend)
Patchseq_Tree_KLdiv = compute_KLdiv(select.cl = select.cl, 
                                    select.cells = query.dat.cells, 
                                    mapping_probability = Tree_mapping_probability, 
                                    memb = memb)


######################################################################################################
### Correlation ######################################################################################
######################################################################################################


FACs_anno <- Read_anno_cellonly_region(paste0(ref.dir, "anno.feather"), region = "VISp")
FACs_anno <- FACs_anno[FACS.cells, ] 
FACs.cl.med <- Compute_median_gene_expression(anno_file = FACs_anno, norm.dat = norm.dat , markers = select.markers)
save(FACs.cl.med, file=file.path(paste0("Final_matrices/FACs.cl.med", ".rda")))
load("Final_matrices/FACs.cl.med.rda")

Patchseq_FACs_cor <- Compute_correlation_mat(markers = select.markers, cells = colnames(query.dat.norm), #locked_cells_sample_id, 
                                             query.dat.norm = query.dat.norm, train.cl.med = FACs.cl.med)


######################################################################################################
### Aggreagte results ################################################################################
######################################################################################################

Patchseq_Tree_memb <- Patchseq_Tree_memb[query.dat.cells, select.cl]

Tree_3_cls <- Get_3_best_cl(Patchseq_Tree_memb, select.cl)
colnames(Tree_3_cls) <- c("Tree_first_cl", "Tree_second_cl", "Tree_third_cl")

Tree_3_bts <- Get_3_best_bt(Patchseq_Tree_memb, select.cl)
colnames(Tree_3_bts) <- c("Tree_first_bt", "Tree_second_bt", "Tree_third_bt")

Tree_3_KL <- Get_3_best_KL(memb = Patchseq_Tree_memb, ref.cl = select.cl, KLdiv = Patchseq_Tree_KLdiv)
colnames(Tree_3_KL) <-c("Tree_first_KL", "Tree_second_KL", "Tree_third_KL")

Tree_3_cor <- Get_3_best_cor(memb = Patchseq_Tree_memb, ref.cl = select.cl, cor = Patchseq_FACs_cor)
colnames(Tree_3_cor) <- c("Tree_first_cor", "Tree_second_cor", "Tree_third_cor")

cells <- query.dat.cells
cells <- locked_cells_sample_id
results <- cbind.data.frame(Tree_3_cls[cells,],
                            Tree_3_bts[cells,],
                            Tree_3_KL[cells,],
                            Tree_3_cor[cells,])


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
  column_to_rownames("id") 

results <- results[,c(Original_cols, "Tree_call", "Tree_color", "Tree_id")]


######################################################################################################
### Write outputs ####################################################################################
######################################################################################################

write.csv(results,"Final_matrices_locked_data/Patchseq_results.csv")
# combine with reference annotations
#rownames(samp.dat)=samp.dat$patched_cell_container
#mapping.df = cbind(mapping.df, samp.dat[row.names(mapping.df),])
#save(mapping.df, file=file.path(paste0("mapping.df", bp.name.add, ".rda")))
#save(memb, file=file.path(paste0("mapping.memb",bp.name.add,".rda")))
#write.csv(mapping.df, file=file.path(paste0("mapping.df", bp.name.add, ".csv")))
#write.csv(memb, file=file.path( paste0("mapping.memb",bp.name.add, ".csv")))

######################################################################################################
### Functions used ###################################################################################
######################################################################################################

compute_KLdiv <- function(select.cl, select.cells, memb, mapping_probability){
  library("LaplacesDemon")
  memb <- memb[select.cells, select.cl]
  mapping_probability <- mapping_probability[select.cl, select.cl]
  KLdiv <- matrix(nrow = length(select.cells), ncol = length(select.cl))
  mapping_probability[mapping_probability==0] <- 0.0000001
  mapping_probability <- mapping_probability/rowSums(mapping_probability)
  memb[memb==0] <- 0.0000001
  memb <- memb/rowSums(memb)
  
  for (cell in 1:length(select.cells)) {
    for (cl in 1:length(select.cl)){
      P <- memb[cell,select.cl]
      Q <- mapping_probability[cl, select.cl]
      KLdiv[cell, cl] <- KLD(px= P,py =Q)$sum.KLD.px.py
    }
  }
  
  rownames(KLdiv) <- rownames(memb)
  colnames(KLdiv) <- as.character(select.cl)
  KLdiv
}


Read_anno_cellonly_region <- function(annopath, region){
  anno <- read_feather(annopath)
  region = region
  anno <- anno %>% filter(class_label %in% c("GABAergic", "Glutamatergic") & region_label ==region)
  anno <- as.data.frame(anno)
  rownames(anno) <- anno$sample_id
  anno
}

Compute_median_gene_expression <- function(anno_file, norm.dat, markers) {
  train.cl.med <- do.call("cbind", tapply(anno_file$sample_id, anno_file$cluster_label, function(x) rowMedians(norm.dat[markers, x, drop=F])))
  rownames(train.cl.med) <- markers
  train.cl.med
}

Compute_correlation_mat <- function(markers, cells, query.dat.norm, train.cl.med){
  cormat = cor(as.matrix(query.dat.norm[markers,cells]), train.cl.med[markers,])
  cormat
}


Get_nth_predicted_cl <- function(memb, ref.cl, nth_cl){
  sorted = t(apply(memb[ , ref.cl], 1, function(x) names(x)[order(x, decreasing = TRUE)]))
  sorted[, nth_cl]
}

Get_nth_predicted_bt <- function(memb, ref.cl, nth_bt){
  sorted = t(apply(memb[ ,ref.cl], 1, function(x) x[order(x, decreasing = TRUE)]))
  sorted[, nth_bt]
}

Get_3_best_cor <- function(memb, ref.cl, cor){
  first = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 1)
  first = Get_matrix_member(first, cor)
  second = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 2)
  second = Get_matrix_member(second, cor)
  third = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 3)
  third = Get_matrix_member(third, cor)
  cbind(first, second, third)
}

Get_3_best_KL <- function(memb, ref.cl, KLdiv){
  first = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 1)
  first = Get_matrix_member(first, KLdiv)
  second = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 2)
  second = Get_matrix_member(second, KLdiv)
  third = Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 3)
  third = Get_matrix_member(third, KLdiv)
  cbind(first, second, third)
}

Get_3_best_cl <- function(memb, ref.cl){
  first <- Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 1)
  second <- Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 2)
  third <- Get_nth_predicted_cl(memb = memb, ref.cl = ref.cl, nth_cl = 3)
  cbind(first, second, third)
}

Get_3_best_bt <- function(memb, ref.cl){
  first <- Get_nth_predicted_bt(memb = memb, ref.cl = ref.cl, nth_bt = 1)
  second <- Get_nth_predicted_bt(memb = memb, ref.cl = ref.cl, nth_bt = 2)
  third <- Get_nth_predicted_bt(memb = memb, ref.cl = ref.cl, nth_bt = 3)
  cbind(first, second, third)
}


