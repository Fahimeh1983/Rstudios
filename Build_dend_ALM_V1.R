library(scrattch.hicat)
library(feather)
library(dplyr)
source(paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R"))
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/scrattch.hicat/R/dend.markers.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/scrattch.hicat/R/markers.R")
source("/allen//programs/celltypes/workgroups/rnaseqanalysis/yzizhen/My_R/patchseq/patchseq.R")


##########################################################################################
### Setting up the path dirs: ############################################################
##########################################################################################

ref.data.rda.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2/")
work.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/Patchseq_vs_FACs_cre_analysis/mouse_patchseq_VISp_20181220_collapsed40_cpm/"
anno.path = "/data/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather"

##########################################################################################
### Reading REF FACS data, dend and markers list: ########################################
##########################################################################################

#FACs_anno <- Read_anno_cellonly_region(paste0(facs.dir, "anno.feather"), region = "VISp")
load(paste0(ref.data.rda.path, "norm.dat.rda"))
load(paste0(ref.data.rda.path, "Tasic_2018_dend.rda"))
FACS.anno <- read_feather(anno.path)

# For making the tree for VISp cells 
FACS.anno <- as.data.frame(FACS.anno) %>% 
  filter(region_label == "VISp") %>% 
  filter(!grepl("ALM", cluster_label)) %>% 
  filter(class_label !="Low Quality")
dim(FACS.anno)

rownames(FACS.anno) <- FACS.anno$sample_id
cl <- setNames(FACS.anno$cl, FACS.anno$sample_id)
cl <- factor(cl)
cl.df <- unique(as.data.frame(FACS.anno[,c("cluster_id", "cluster_label", "cluster_color", "subclass_id", "subclass_label", "class_id", "class_label", "cl")]))
rownames(cl.df) <- cl.df$cl
cl.df <- cl.df[order(cl.df$cluster_id),]
cl.df <- cl.df[cl.df$class_label!="Low Quality",]
cl <- cl[cl%in% cl.df$cl]
cl <- droplevels(cl)

dend
dim(cl.df)
plot(dend)
norm.dat <- norm.dat[, names(cl)]

##########################################################################################
### Building dend that has markers attached:  ############################################
##########################################################################################

de.param = de_param(low.th=0.6, padj.th=0.01, lfc.th=1, q1.th=0.5, q2.th=NULL,q.diff.th=0.7, de.score.th=80, min.cells=4)
display_result = display_cl(cl, norm.dat, n.markers=50, de.param = de.param, prefix = "test")
#save(display_result, file =paste0(work.dir, "display_results.rda") )
#load(paste0(work.dir, "display_results.rda"))
#de.genes= display_result$de.genes

#tmp=get_gene_score(de.genes, top.n=100, max.num=4000, bin.th=4) 
#up.gene.score=tmp$up.gene.score
#down.gene.score=tmp$down.gene.score

###Find markers of dendrogram.

de.genes=display_result$de.genes
tmp=get_gene_score(de.genes, top.n=100, max.num=2000, bin.th=4) 
up.gene.score=tmp$up.gene.score
down.gene.score=tmp$down.gene.score
select.genes = row.names(up.gene.score)
#load(paste0(workdir, "V1_ALM_dend_markers_attached.rda"))
n.markers = 30 # use up to 30 markers per node

labels(dend)[48] <- "L6b Col8a1 Rprm"
missing_labels <- labels(dend)[!labels(dend) %in% cl.df$cluster_label]
dend_pruned <- prune_dend(dend = dend, rm.labels = missing_labels)
#labels(dend)[41] <- "L6 CT Nxph2 Sla"

#labels(dend) = row.names(cl.df)[match(labels(dend),cl.df$cluster_label)]

dend = select_dend_markers(dend_pruned, norm.dat=norm.dat, cl=cl, de.genes=de.genes,
                           up.gene.score=up.gene.score[select.genes,], 
                           down.gene.score=down.gene.score[select.genes,], n.markers=n.markers)
dend = select_pos_dend_markers(dend= dend, norm.dat = norm.dat, cl = cl)

#save(dend, file =paste0(work.dir, "V1_ALM_dend_markers_attached.rda"))
#load(paste0(work.dir, "V1_ALM_dend_markers_attached.rda") )

##########################################################################################
### Mapping FACS on to FACS data:  #######################################################
##########################################################################################
# We have built the dend with the norm.dat which does not have the bad cells
# Now we are going to map the norm.dat with all the bad cells onto the tree

source(paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/patchseq.R"))
library(scrattch.hicat)

Good_norm.dat <- norm.dat
load(paste0(ref.data.rda.path, "norm.dat.rda"))
#norm.dat <- as.matrix(norm.dat)
FACs.cells <- rownames(FACS.anno)
norm.dat <- norm.dat[, FACs.cells]
norm.dat <- as.matrix(norm.dat)
Good_norm.dat <- as.matrix(Good_norm.dat)

FACS_Tree_memb = map_dend_membership(dend, cl=cl, Good_norm.dat, norm.dat, FACs.cells, bs.num=100, p=0.7, low.th=0.15)
#save(FACS_Tree_memb, file= paste0(work.dir, "/FACS_ALM_V1_TREE.rda"))
#load(paste0(work.dir, "/FACS_ALM_V1_TREE.rda"))

FACS_Tree_mapping.df <- summarize_cl(dend, FACS_Tree_memb, norm.dat, conf.th=0.7, min.genes=1, min.genes.ratio=0.3)
FACS_Tree_mapping.df$sample_id <- rownames(FACS_Tree_mapping.df)
cl.df$cl <- as.factor(cl.df$cl)
FACS_Tree_mapping.df <- left_join(FACS_Tree_mapping.df, cl.df)
FACS_Tree_mapping.df <- FACS_Tree_mapping.df %>% mutate(cluster_label = ifelse(is.na(cluster_label), cl, cluster_label))
rownames(FACS_Tree_mapping.df) <- FACS_Tree_mapping.df$sample_id
FACS_Tree_mapping.df <- FACS_Tree_mapping.df[,setdiff(colnames(FACS_Tree_mapping.df), "sample_id")]
save(FACS_Tree_mapping.df, file= paste0(work.dir, "/FACS_ALM_V1_TREE_mapping.df.rda"))
#load(paste0(work.dir, "/FACS_ALM_V1_TREE_mapping.df.rda"))



