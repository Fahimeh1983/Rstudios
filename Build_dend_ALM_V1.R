library(scrattch.hicat)
source(paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R"))


##########################################################################################
### Setting up the path dirs: ############################################################
##########################################################################################

ref.data.rda.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2/")
workdir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/Patchseq_vs_FACs_cre_analysis/mouse_patchseq_VISp_20181220_collapsed40_cpm/"

##########################################################################################
### Reading REF FACS data, dend and markers list: ########################################
##########################################################################################

#FACs_anno <- Read_anno_cellonly_region(paste0(facs.dir, "anno.feather"), region = "VISp")
load(paste0(ref.data.rda.path, "norm.dat.rda"))
load(paste0(ref.data.rda.path, "Tasic_2018_dend.rda"))
load(paste0(ref.data.rda.path, "Tasic_2018_cldf.rda"))
cl.df <- as.data.frame(cl.df)
library(tidyr)
library(dplyr)
cl.df <- cl.df[, c("cluster_id", "cluster_label", "cluster_color", "subclass_id", "subclass_label", "class_id", "class_label")]
cl.df <- cl.df[cl.df$class_label!="Low Quality",]
cl.df$cl <- rownames(cl.df)
plot(dend)
#cl <- Renew_list(ls = cl, ref.df = cl.df, label = "cluster_id", new.label = "cluster_label")
norm.dat <- norm.dat[, names(cl)]

##########################################################################################
### Building dend that has markers attached:  ############################################
##########################################################################################

de.param = de_param(low.th=0.6, padj.th=0.01, lfc.th=1, q1.th=0.5, q2.th=NULL,q.diff.th=0.7, de.score.th=80, min.cells=4)
display_result = display_cl(cl, norm.dat, n.markers=50, de.param = de.param, prefix = "test")
save(display_result, file =paste0(workdir, "display_results.rda") )
de.genes= display_result$de.genes

tmp=get_gene_score(de.genes, top.n=100, max.num=2000, bin.th=4) 
up.gene.score=tmp$up.gene.score
down.gene.score=tmp$down.gene.score

###Find markers of dendrogram.

de.genes=display_result$de.genes
tmp=get_gene_score(de.genes, top.n=100, max.num=2000, bin.th=4) 
up.gene.score=tmp$up.gene.score
down.gene.score=tmp$down.gene.score
select.genes = row.names(up.gene.score)
load(paste0(workdir, "V1_ALM_dend_markers_attached.rda"))
n.markers = 30 # use up to 30 markers per node
labels(dend) =row.names(cl.df)[match(labels(dend),cl.df$cluster_label)]
dend = select_dend_markers(dend, norm.dat=norm.dat, cl=cl, de.genes=de.genes,
                           up.gene.score=up.gene.score[select.genes,], 
                           down.gene.score=down.gene.score[select.genes,], n.markers=n.markers)
dend = select_pos_dend_markers(dend= dend, norm.dat = norm.dat, cl = cl)

#save(dend, file =paste0(workdir, "V1_ALM_dend_markers_attached.rda") )
load(paste0(workdir, "V1_ALM_dend_markers_attached.rda") )

##########################################################################################
### Mapping FACS on to FACS data:  #######################################################
##########################################################################################

norm.dat <- as.matrix(norm.dat)
FACs.cells <- colnames(norm.dat)
FACS_Tree_memb = map_dend_membership(dend, cl=droplevels(cl[cl %in% labels(dend)]), norm.dat, norm.dat, FACs.cells, bs.num=100, p=0.7, low.th=0.15)
FACS_Tree_mapping.df <- summarize_cl(dend, FACS_Tree_memb, norm.dat, conf.th=0.7, min.genes=1, min.genes.ratio=0.3)
#save(FACS_Tree_memb, file= paste0(workdir, "/FACS_ALM_V1_TREE.rda"))
#save(FACS_Tree_mapping.df, file= paste0(workdir, "/FACS_ALM_V1_TREE_mapping.df.rda"))
#load(paste0(workdir, "/FACS_ALM_V1_TREE.rda"))
#load(paste0(workdir, "/FACS_ALM_V1_TREE_mapping.df.rda"))
colnames(FACS_Tree_memb)[colnames(FACS_Tree_memb) %in% cl.df$cl] <- cl.df[cl.df$cl %in% colnames(FACS_Tree_memb), "cluster_label"]
FACS_Tree_mapping.df <- left_join(FACS_Tree_mapping.df, cl.df)
FACS_Tree_mapping.df <- FACS_Tree_mapping.df %>% mutate(cluster_label = ifelse(is.na(cluster_label), cl, cluster_label))
rownames(FACS_Tree_mapping.df) <- colnames(norm.dat)
rownames(FACS_Tree_memb) <- colnames(norm.dat)



