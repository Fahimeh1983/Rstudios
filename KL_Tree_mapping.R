library(dendextend)
library(matrixStats)
library(feather)
library(tibble)
source("patchseq/heatmap.R")
source("patchseq/de.genes.R")
source("patchseq/dendro.R")
source("patchseq/patchseq.R")
source("patchseq/Mapping_helper_functions.R")

### TODO: The batch date should be updated everytime
#just update batch_date and source it
#batch_date="20190506_BT014-RSC-206"

######################################################################################################
### Setting up some paths ############################################################################
######################################################################################################
### TODO: Change these if you changed them in the mapping.R
ref.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"
res.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"
query.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"

######################################################################################################
### Reading ref data #################################################################################
######################################################################################################
### TODO: the following two lines should be modified if the number of FACS cells has changed  
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/REF_mapping_probability.rda")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/select.markers.rda")

bp.collapse.th = 40
bp.name.add = NULL
if (!is.null(bp.collapse.th)) {
  bp.name.add = paste0(".with.bp.", bp.collapse.th)
}

###load reference data and tree
tmp.load1 = load(file=file.path(res.dir, "ref.data.rda")) # should include cl, cl.df, norm.dat. # The loaded cl is not used because it only includes cluster ids but not cluster labels 
tmp.load2 = load(file.path(file=res.dir, file=paste0("V1.dend", bp.name.add,".rda"))) # should include the pruned V1 tree
tmp.load3 = load(file.path(res.dir, file=paste0("V1.dend.list", bp.name.add,".rda"))) # should include dend.list

plot(dend)

rownames(cl.df)=cl.df$cluster_id
cltmp=cl.df[as.character(cl),"cluster_label"]
names(cltmp)=names(cl)
cl=factor(cltmp)

######################################################################################################
### Loading the query data ###########################################################################
######################################################################################################

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
#Patchseq_Tree_memb = map_dend_membership(dend, cl, norm.dat, query.dat.norm, query.dat.cells, bs.num=100, p=0.7, low.th=0.15)
load(paste0("mapping.memb",bp.name.add,".rda"))
Patchseq_Tree_memb <- memb
rm(memb)

######################################################################################################
### KLdiv ############################################################################################
######################################################################################################

set.seed(1983)
select.cl = labels(dend)
Patchseq_Tree_KLdiv = compute_KLdiv(select.cl = select.cl, 
                                    select.cells = query.dat.cells, 
                                    mapping_probability = Tree_mapping_probability, 
                                    memb = Patchseq_Tree_memb)


######################################################################################################
### Correlation ######################################################################################
######################################################################################################

FACs_anno <- Read_anno_cellonly_region(paste0(ref.dir, "anno.feather"), region = "VISp")
FACs_anno <- FACs_anno[FACS.cells, ] 
FACs.cl.med <- Compute_median_gene_expression(anno_file = FACs_anno, norm.dat = norm.dat , markers = select.markers)

Patchseq_FACs_cor <- Compute_correlation_mat(markers = select.markers, cells = query.dat.cells,  
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
  column_to_rownames("id") 

results <- results[,c(Original_cols, "Tree_call")]


######################################################################################################
### Write outputs ####################################################################################
######################################################################################################

# combine with reference annotations
write.csv(results, file=file.path(paste0("KL_tree_results", bp.name.add, ".csv")))
KL_tree_results <- results
rm(results)
