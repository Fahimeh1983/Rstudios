# Question: Why for some specific cres, we have many data on FACs assigned to an specific cluster but there are not present in Patchseq
# For this, I need to read the anno file 

setwd("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/Patchseq_vs_FACs_cre_analysis/mouse_patchseq_VISp_20181220_collapsed40_cpm")
library(dendextend)
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/heatmap.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/de.genes.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/dendro.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20181220_collapsed40_cpm/patchseq/patchseq.R")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Plot_tree_functions.R")
library(matrixStats)
library(feather)

#just update batch_date and source it
batch_date="20181220_BT014-RSC-188"

blue.red <-colorRampPalette(c("blue", "white", "red"))
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

ref.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520"
res.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"

bp.collapse.th = 40
bp.name.add = NULL
if (!is.null(bp.collapse.th)) {
  bp.name.add = paste0(".with.bp.", bp.collapse.th)
}

###################################################
##########  load reference data and tree ##########
# make sure to load the files created by build_reference_dend
tmp.load1 = load(file=file.path(res.dir, "ref.data.rda")) # should include cl, cl.df, norm.dat. # The loaded cl is not used because it only includes cluster ids but not cluster labels 
tmp.load2 = load(file.path(file=res.dir, file=paste0("V1.dend", bp.name.add,".rda"))) # should include the pruned V1 tree
tmp.load3 = load(file.path(res.dir, file=paste0("V1.dend.list", bp.name.add,".rda"))) # should include dend.list

plot(dend)

rownames(cl.df)=cl.df$cluster_id
cltmp=cl.df[as.character(cl),"cluster_label"]
names(cltmp)=names(cl)
cl=factor(cltmp)

###################################################
################# load query data #################
#(either load from R objects or from an already existed shiny folder or use the FACS data to map to itself)

query.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
tmp<-load(paste0(query.dir,batch_date,"_mouse_patchseq_star2.0_cpm.Rdata"))
query.dat = cpmR
#patchseq_anno <- read_feather("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/anno.feather")
#patchseq_anno <- as.data.frame(patchseq_anno)
#rownames(patchseq_anno) <- patchseq_anno$sample_id

# loading patchseq samp.dat object
tmp<-load(paste0(query.dir,batch_date,"_mouse_patchseq_star2.0_samp.dat.Rdata"))
keepcells = which(samp.dat$Region=="VISp" & samp.dat$Type=="patch_seq")
samp.dat = samp.dat[c(keepcells, which(samp.dat$Region=="TCx"),which(samp.dat$Region=="FCx")   ),]   #FCx is for Brian.  Rat samples mapped in mouse
query.dat = query.dat[,as.character(samp.dat$exp_component_name)]
colnames(query.dat)=as.character(samp.dat$patched_cell_container)
query.dat.norm = log2(as.matrix(query.dat+1))
idx=match(rownames(norm.dat), rownames(query.dat.norm))
query.dat.norm=query.dat.norm[idx,]

###################################################
######### loading FACs samp.dat object ############

FACs.anno.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/Mm_VISp_14236_20180912/"
FACs.samp.dat <- read_feather(paste0(FACs.anno.dir,"anno.feather"))
FACs.cells = colnames(norm.dat)
FACs.samp.dat <- as.data.frame(FACs.samp.dat)
rownames(FACs.samp.dat) <- FACs.samp.dat$sample_id
FACs.samp.dat <- FACs.samp.dat[FACs.cells,]

###################################################
########### Some more initializations #############

query.dat.cells = colnames(query.dat.norm)
FACs.cells = colnames(norm.dat)

###################################################
######### Map patchseq data on ref cl #############
memb = map_dend_membership(dend, cl, norm.dat, query.dat.norm, query.dat.cells, bs.num=100, p=0.7, low.th=0.15)
mapping.df = summarize_cl(dend, memb, query.dat.norm, conf.th=0.7, min.genes=0, min.genes.ratio=0)
rownames(samp.dat)=samp.dat$patched_cell_container
mapping.df = cbind(mapping.df, samp.dat[row.names(mapping.df),])
save(mapping.df, file=file.path(paste0("mapping.df", bp.name.add, ".rda")))
save(memb, file=file.path(paste0("mapping.memb",bp.name.add,".rda")))
write.csv(mapping.df, file=file.path(paste0("mapping.df", bp.name.add, ".csv")))
write.csv(memb, file=file.path( paste0("mapping.memb",bp.name.add, ".csv")))
load("mapping.df.with.bp.40.rda")
load("mapping.memb.with.bp.40.rda")

###################################################
########### Map FACseq data on ref cl #############
FACs.memb = map_dend_membership(dend, cl, norm.dat, norm.dat, FACs.cells, bs.num=100, p=0.7, low.th=0.15)
FACs.mapping.df = summarize_cl(dend, FACs.memb, norm.dat, conf.th=0.7, min.genes=0, min.genes.ratio=0.)
FACs.mapping.df = cbind(FACs.mapping.df, FACs.samp.dat[row.names(FACs.mapping.df),])
save(FACs.mapping.df, file=file.path(paste0("FACs.mapping.df", bp.name.add, ".rda")))
save(FACs.memb, file=file.path(paste0("FACs.mapping.memb",bp.name.add,".rda")))
write.csv(FACs.mapping.df, file=file.path(paste0("FACs.mapping.df", bp.name.add, ".csv")))
write.csv(FACs.memb, file=file.path( paste0("FACs.mapping.memb",bp.name.add, ".csv")))
load("FACs.mapping.df.with.bp.40.rda")
load("FACs.mapping.memb.with.bp.40.rda")


###################################################
### Comparing the clustering and mapping for FACs##

#cl and cl.df have all the information from the clustering, even in the anno file, there are two column, if it is "cl" then in means that comes from
#clustering but if it is "cluster_label" then that comes from mapping 
FACs_mapping_vs_clustering_df <- cbind.data.frame(cl[FACs.cells],FACs.mapping.df[FACs.cells,"cl"])
FACs_mapping_vs_clustering_df <- as.data.frame(FACs_mapping_vs_clustering_df)
colnames(FACs_mapping_vs_clustering_df) <- c("clustering_cluster_label", "mapping_cluster_label")
rownames(FACs_mapping_vs_clustering_df) <- FACs.cells

tree_levels <- levels(FACs.mapping.df$cl)
FACs_mapping_vs_clustering_df$mapping_cluster_label <- factor(FACs_mapping_vs_clustering_df$mapping_cluster_label, levels = tree_levels)
FACs_mapping_vs_clustering_df$clustering_cluster_label <- factor(FACs_mapping_vs_clustering_df$clustering_cluster_label, levels = tree_levels)

### Looking and bad and good mapped cells
bad_mapped_cells <- FACs_mapping_vs_clustering_df$clustering_cluster_label != FACs_mapping_vs_clustering_df$mapping_cluster_label
bad_mapped_cells <- rownames(FACs_mapping_vs_clustering_df)[bad_mapped_cells]
good_mapped_cells <- setdiff(FACs.cells, bad_mapped_cells)
cat("Acuraccy of mapping:",length(good_mapped_cells)/length(FACs.cells)*100, "%")

freq_table<- FACs_mapping_vs_clustering_df[bad_mapped_cells,] %>% 
  select(mapping_cluster_label) %>% 
  table() %>% 
  as.data.frame() %>% 
  `colnames<-`(c("label","Freq")) 

dendro_data <- as.ggdend(dend)
dend_nodes_table <- dendro_data$nodes %>% mutate(label = get_nodes_attr(dend, "label"))

table <- FACs_mapping_vs_clustering_df %>% 
  filter(as.character(clustering_cluster_label)!=as.character(mapping_cluster_label)) %>% 
  droplevels() %>% 
  table()
## Plot the heatmap of clustering and mapping comparison for FACs
n_bad_cells_mapped <- list(scores=freq_table[freq_table$label == dend_nodes_table$label,"Freq"])
Plot_info_on_dend_nodes(dend, n_bad_cells_mapped, FALSE, "Number of FACs cells that mapped wrongly at each node")

## Further checking
select.cl <- as.character(unique(cl))
df <- FACs.memb[FACs.cells, select.cl]
sorted_mapped_cl <- t(apply(df, 1, function(x) names(x)[order(x, decreasing = TRUE)])) %>% `colnames<-`(seq(ncol(df)))
sorted_bt_values <- t(apply(df, 1, function(x) x[order(x, decreasing = TRUE)]))  %>% `colnames<-`(seq(ncol(df)))

if (length(good_mapped_cells) == sum(sorted_mapped_cl[good_mapped_cells,"1"] == cl[good_mapped_cells])){
  print("All good cells are mapped correctly to the leafs")
}
print("All the values in the histogram must be more than 0.7")
hist(sorted_bt_values[good_mapped_cells,"1"], xlab = "bt-support", ylab = "number oc cells", main="Histogram of bt-support of good cells")

#If all the bad-mapped cells are not mapped only because their bs-support is lower than 0.7, then we should see in the below:
#There are some cells among these bad_mapped_cells that can be re-mapped if we decrease the bt-support
#These cells are mapped correclt to their cls but since the bt is small, they are rejected
recovery_bad_mapped_cells <- sorted_mapped_cl[bad_mapped_cells,"1"] == cl[bad_mapped_cells]
recovery_bad_mapped_cells <- bad_mapped_cells[recovery_bad_mapped_cells]
hist(sorted_bt_values[recovery_bad_mapped_cells,"1"], xlab = "bt-support", ylab = "number oc cells", main="Histogram of bt-support of bad cells1")

remaining_bad_mapped_cells <- setdiff(bad_mapped_cells, recovery_bad_mapped_cells)
leaf_bad_mapped_cells <- sorted_bt_values[remaining_bad_mapped_cells,"1"]>0.7
leaf_bad_mapped_cells <- remaining_bad_mapped_cells[leaf_bad_mapped_cells]
hist(sorted_bt_values[leaf_bad_mapped_cells,"1"],xlab = "bt-support", ylab = "number oc cells", main="Histogram of bt-support of bad cells2")
#All these selected cells among bad-mapped-cells have a high bt value for a leaf node but they are not the real cl for that cell
# This means that if our hypothesis for facs.mapping.df is correct then the cl in that data frame is simply the sorted mapped cl
# This must be TRUE
length(leaf_bad_mapped_cells) == sum(FACs_mapping_vs_clustering_df[leaf_bad_mapped_cells,"mapping_cluster_label"] == sorted_mapped_cl[leaf_bad_mapped_cells,"1"])
#Also this list should have identical number of cells:length(leaf_bad_mapped_cells)
FACs_mapping_vs_clustering_df[bad_mapped_cells,] %>% filter(mapping_cluster_label %in% select.cl) %>% select(mapping_cluster_label)

higher_node_bad_mapped_cells <- setdiff(remaining_bad_mapped_cells,leaf_bad_mapped_cells)
hist(sorted_bt_values[higher_node_bad_mapped_cells,"1"], xlab = "bt-support", ylab = "number oc cells", main="Histogram of bt-support of bad cells3")

#VERY IMPORTANT: This amount of cells are confused between two leafs
sum(sorted_mapped_cl[higher_node_bad_mapped_cells,"2"]==FACs_mapping_vs_clustering_df[higher_node_bad_mapped_cells,"clustering_cluster_label"])
tmp<-sorted_mapped_cl[higher_node_bad_mapped_cells,"2"]==FACs_mapping_vs_clustering_df[higher_node_bad_mapped_cells,"clustering_cluster_label"]
tmp<-higher_node_bad_mapped_cells[tmp]
FACs_mapping_vs_clustering_df[tmp,"mapping_cluster_label"]
sort(table(FACs_mapping_vs_clustering_df[tmp,"mapping_cluster_label"]))
confused_cl_pairs <- apply(sorted_mapped_cl[tmp,c("1","2")],1,paste,collapse="-")
sort(table(confused_cl_pairs))

confused_cl_combined_bt <- apply(sorted_bt_values[tmp,c("1","2")],1,sum)

library(reshape2)
Freq_table <- FACs_mapping_vs_clustering_df %>% 
  group_by(clustering_cluster_label, mapping_cluster_label) %>% 
  summarise(Freq = n())
mapping_levels <- get_nodes_attr(dend,"label")
clustering_levels <- levels(cl)
mapping_levels <- c(clustering_levels, mapping_levels[!(mapping_levels %in% clustering_levels)])
node_levels <- mapping_levels[!(mapping_levels %in% clustering_levels)]
tmp <- expand.grid(clustering_levels, mapping_levels)
colnames(tmp) <- c("clustering_levels", "mapping_levels")
colnames(Freq_table) <- c("clustering_levels", "mapping_levels", "Freq")
expanded <- merge(tmp, Freq_table, by= c("clustering_levels", "mapping_levels"),all = T)
expanded[is.na(expanded)] <- 0
Freq_matrix <- matrix(expanded$Freq, nrow = length(clustering_levels), ncol = length(mapping_levels))
colnames(Freq_matrix) <- as.character(mapping_levels)
rownames(Freq_matrix) <- as.character(clustering_levels)
for (c in colnames(Freq_matrix)){
  for (r in rownames(Freq_matrix)) {
    Freq_matrix[r,c] <- expanded[(expanded$clustering_levels== r & expanded$mapping_levels== c),"Freq"]
  }
}
Freq_matrix <- Freq_matrix/rowSums(Freq_matrix)

ggplot(data = melt(Freq_matrix[clustering_levels,node_levels]), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Clustering lable") + ylab("Mapping lables") +  
  scale_fill_gradient(low = "white", high = "red")

ggplot(data = melt(Freq_matrix[clustering_levels, clustering_levels]), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=5)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Clustering lable") + ylab("Mapping lables") +  
  scale_fill_gradient(low = "white", high = "red")

###################################################
##### Building the mapping_probability matrix #####

library(data.table)
FACs.memb <- as.data.frame.matrix(FACs.memb) 
FACs.memb$sample_id <- rownames(FACs.memb)
FACs_mapping_vs_clustering_df$sample_id <- rownames(FACs_mapping_vs_clustering_df)
membership_freq_table <- as.data.frame(left_join(FACs.memb[,c(select.cl,"sample_id")], FACs_mapping_vs_clustering_df,"sample_id"))
rownames(membership_freq_table) <- membership_freq_table$sample_id
membership_freq_table <- membership_freq_table[FACs.cells,]

tmp <- do.call("rbind",tapply(row.names(membership_freq_table),membership_freq_table$clustering_cluster_label, function(x)colSums(membership_freq_table[x,select.cl])))
membership_freq_table <- tmp
membership_freq_table <- membership_freq_table/rowSums(membership_freq_table)
missing_rows <- colnames(membership_freq_table)[!colnames(membership_freq_table) %in% rownames(membership_freq_table)]

ggplot(data = melt(data.matrix(membership_freq_table[clustering_levels, clustering_levels])), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Clustering lable") + ylab("Mapping lables") +  
  scale_fill_gradient(low = "white", high = "red")

mapping_probability <- membership_freq_table
save(mapping_probability, file=file.path(paste0("mapping_probability",bp.name.add,".rda")))
write.csv(mapping_probability, file=file.path(paste0("mapping_probability", bp.name.add, ".csv")))
load("mapping_probability.with.bp.40.rda")


###################################################
####### Looking into one suspicious cluster #######
suspicious_cluster <- rownames(mapping_probability)[diag(mapping_probability)<0.5]
select.cells <- names(cl)[cl %in% suspicious_cluster]
# df <- FACs.memb[select.cells,select.cl]
# sorted_mapped_cl <- t(apply(df, 1, function(x) names(x)[order(x, decreasing = TRUE)])) %>% `colnames<-`(paste0("cl",seq(ncol(df))))
# sorted_bt_values <- t(apply(df, 1, function(x) x[order(x, decreasing = TRUE)]))  %>% `colnames<-`(paste0("cl",seq(ncol(df))))

dend.list = dend_list(dend)
attr(dend.list[["n35"]],"markers")
attr(dend.list[["n82"]],"markers.byCl")
IEG <-read.table("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/IEG_genes_analysis/V1_ALM/process_24411/Early_late_genes.csv", sep = ",", header = TRUE)
early_genes <- IEG$Gene_name[IEG$category=="early"]
#These genes are among early genes and we think this cluster needs to be combined with the other one
names(attr(dend.list[["n35"]],"markers"))[names(attr(dend.list[["n35"]],"markers")) %in% early_genes]



###################################################
## Using mapping_probability matrix on Patchseq ###

select.cl <- as.character(unique(cl))
select.cells <- query.dat.cells
mapping_probability <- mapping_probability[select.cl, select.cl]
df <- memb[select.cells,select.cl]
sorted_mapped_cl <- t(apply(df, 1, function(x) names(x)[order(x, decreasing = TRUE)])) %>% `colnames<-`(paste0("cl",seq(ncol(df))))
sorted_bt_values <- t(apply(df, 1, function(x) x[order(x, decreasing = TRUE)]))  %>% `colnames<-`(paste0("cl",seq(ncol(df))))
sorted_bt_values <- as.data.frame(sorted_bt_values)
sorted_mapped_cl <- as.data.frame.matrix(sorted_mapped_cl)
sorted_bt_values <- sorted_bt_values[select.cells,]
sorted_mapped_cl <- sorted_mapped_cl[select.cells,]



df <- as.data.frame(apply(sorted_bt_values, 1, function(x){colnames(sorted_bt_values)[cumsum(x) >0.7][1]}))
colnames(df) <- c("stop_cl")
df <- as.numeric(gsub('cl([0-99]+)*','\\1',df$stop_cl))
df <- as.data.frame(df)
colnames(df) <- c("stop_cl")
rownames(df) <- rownames(sorted_bt_values)
sorted_bt_values <- cbind(sorted_bt_values, df)

for (row in rownames(sorted_bt_values)){
    cl_stop = sorted_bt_values[row,"stop_cl"]
    cl_list <- sorted_mapped_cl[row,1:cl_stop]
    prop_list <- list()
    if (length(cl_list)==1){
       df[row, "mapping_quality"] <- mapping_probability[as.character(cl_list[[1]]),as.character(cl_list[[1]])]
    }else{
      for (i in 2:length(cl_list)){
           prop_list <- c(prop_list, max(mapping_probability[as.character(cl_list[[1]]),as.character(cl_list[[i]])], mapping_probability[as.character(cl_list[[i]]), as.character(cl_list[[1]])]))
      }
      df[row,"mapping_quality"] <- min(unlist(prop_list))
    }
}

mapping.df[select.cells,"mapping_quality"] <- df[select.cells, "mapping_quality"]
mapping.df[select.cells,"stop_cl"] <- df[select.cells, "stop_cl"]
mapping.df$mapping_quality <- as.numeric(mapping.df$mapping_quality)
hist(mapping.df$mapping_quality, breaks = 50, xlab = "mapping_quality", main = "Histogram of mapping_quality")
hist(mapping.df$resolution.index, breaks = 50, xlab = "mapping_quality", main = "Histogram of resolution index")
low_quality_cells <- rownames(mapping.df)[mapping.df$mapping_quality==0]
high_quality_cells <- setdiff(query.dat.cells, low_quality_cells)
low_resolution_cells <- rownames(mapping.df)[mapping.df$resolution.index!=1]
high_resolution_cells <- setdiff(query.dat.cells, low_resolution_cells)
mapping.df[low_quality_cells, "mapping_quality_cat"] <- c("low_quality_cells") 
mapping.df[high_quality_cells, "mapping_quality_cat"] <- c("high_quality_cells") 
mapping.df[low_resolution_cells, "resolution_quality_cat"] <- c("low_resolution_cells") 
mapping.df[high_resolution_cells, "resolution_quality_cat"] <- c("high_resolution_cells") 


ggplot(mapping.df, aes(resolution.index, fill = mapping_quality_cat)) +
  geom_histogram(alpha = 0.5, position = 'identity', binwidth = 0.02) +
  xlab("resolution index") +
  ylab("Number of cells")

ggplot(mapping.df, aes(mapping_quality, fill = resolution_quality_cat)) +
  geom_histogram(alpha = 0.5, position = 'identity', binwidth = 0.02) +
  xlab("Mapping quality") +
  ylab("Number of cells")

ggplot(mapping.df, aes(mapping_quality, fill = mapping_quality_cat)) +
  geom_histogram(alpha = 0.5, position = 'identity', binwidth = 0.02) +
  xlab("Mapping quality") +
  ylab("Number of cells")


#None of the cells that have a mapping_quality==0 has a resolution index of 1
select.cells <- rownames(mapping.df)[mapping.df$mapping_quality==0]
hist(mapping.df[select.cells,"resolution.index"])

#one group of those with mapping_quality==0 are confused between two clusters
# 27 cells confused between "Sst Crh 4930553C11Rik ", "Sst Crhr2 Efemp1"0.9106318, 0.3173333
select.cells <- rownames(mapping.df)[(mapping.df$mapping_quality>0.3 & mapping.df$resolution.index<1)]
mapping.df[select.cells,"resolution.index"]
mapping.df[select.cells,"mapping_quality"]
get_pair_matrix(mapping_probability, sorted_mapped_cl[select.cells,1], sorted_mapped_cl[select.cells,2])

#17 cell confused between L5 IT VISp Hsd11b1 Endou     L5 IT VISp Whrn Tox2 0.9144093, 0.2086957
#83 between Sst Calb2 Pdlim5         Sst Calb2 Necab1 0.8441000 0.2131667
select.cells <- rownames(mapping.df)[(mapping.df$mapping_quality>0.2  & mapping.df$mapping_quality < 0.3)]
mapping.df[select.cells,"resolution.index"]
mapping.df[select.cells,"mapping_quality"]


#24 cells confused 
select.cells <- rownames(mapping.df)[(mapping.df$mapping_quality<0.01 & mapping.df$resolution.index >0.9)]
mapping.df[select.cells, "cl"]#!!!!!!!!!!!!!!!!!!!# Ask Zizhen
hist(mapping.df[select.cells,"resolution.index"], breaks = 50)
mapping.df[select.cells,"mapping_quality"]




select.cells <- rownames(mapping.df)[mapping.df$mapping_quality == 0]
hist(mapping.df[select.cells,"resolution.index"],breaks = 50, xlab = "resolution.index", main = "histogram of resolution index")
sum(mapping.df[select.cells,"resolution.index"]==1)


select.cells <- rownames(mapping.df)[(mapping.df$mapping_quality == 0 & mapping.df$resolution.index >0.8)]
hist(mapping.df[select.cells,"resolution.index"],breaks = 50, xlab = "resolution.index", main = "histogram of resolution index")
sum(mapping.df[select.cells,"resolution.index"]==1)


first_cl <- apply(memb[,select.cl], 1, function(x) names(x)[order(x, decreasing = TRUE)[1]])
first_cl <- data.frame(first_cl)
memb <- as.data.frame.matrix(memb)
rownames(first_cl) <- rownames(memb)
tmp <- cbind(first_cl[select.cells,], memb[select.cells,select.cl])
colnames(tmp) <- c("first_cl",select.cl)
df<-do.call("rbind",tapply(row.names(tmp),tmp$first_cl, function(x)colMeans(tmp[x,select.cl])))
non_missing_rows <- rownames(df)
missing_rows <- setdiff(colnames(df), rownames(df))
df <- rbind(df,t(as.data.frame(sapply(missing_rows, function(x)rep(0,length(select.cl)), simplify = F))))
rownames(df) <- c(non_missing_rows, missing_rows)
patchseq_mapping_probability <- df 

ggplot(data = melt(data.matrix(patchseq_mapping_probability[select.cl,select.cl])), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Mapped cluster") + ylab("all clusters") +  
  scale_fill_gradient(low = "white", high = "red")

save(patchseq_mapping_probability, file=file.path(paste0("patchseq_mapping_probability",bp.name.add,".rda")))
write.csv(patchseq_mapping_probability, file=file.path(paste0("patchseq_mapping_probability", bp.name.add, ".csv")))
load("patchseq_mapping_probability.with.bp.40.rda")

sub <- mapping_probability[select.cl, select.cl] - patchseq_mapping_probability[select.cl, select.cl]
ggplot(data = melt(data.matrix(sub[select.cl,select.cl])), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  #xlab("First_cluster") + ylab("Mapping") +  
  scale_fill_gradient(low = "white", high = "red")


cl.name.mat = sapply(sorted_mapped_cl, as.character)
row.names(cl.name.mat) = row.names(sorted_mapped_cl)
val.mat = as.matrix(sorted_bt_values)
sapply(select.cells, function(x)setNames(val.mat[x,], cl.name.mat[x,]), simplify=F)


# low_quality_cells <- rownames(df)[df[,"map_prop"]==0]
# other_cells <- query.dat.cells[!query.dat.cells %in% c(low_quality_cells)]
# 
# 
# mapping.df[low_quality_cells, "mapping_quality"] <- "low_quality_cells"
# mapping.df[other_cells, "mapping_quality"] <- "other_cells"

# ggplot(mapping.df, aes(resolution.index, fill = resolution.index)) + 
#   geom_histogram(alpha = 0.5, position = 'identity', binwidth = 0.05) +
#   xlab("res.index") +
#   ylab("Number of cells")
# 
# ggplot(mapping.df, aes(resolution.index, fill = mapping_quality)) + 
#   geom_histogram(alpha = 0.5, position = 'identity', binwidth = 0.05) +
#   xlab("res.index") +
#   ylab("Number of cells")


#Cells that have high res.index but they are still low quality based on our analysis
select.cells <- low_quality_cells[mapping.df[low_quality_cells,"resolution.index"]==1]
head(sort(memb[select.cells[1],select.cl], decreasing = TRUE))
mapping.df[select.cells[1],"cl"]













