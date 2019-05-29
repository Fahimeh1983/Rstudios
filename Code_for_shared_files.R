######################################################################################################
### Code used to share data with other Gabe ##########################################################
######################################################################################################

last_shiny <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20190506_test/"
Fahimehdir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/patchseq-work-dir/Patchseq_vs_FACs_cre_analysis/mouse_patchseq_VISp_20181220_collapsed40_cpm/"
KL_tree_results <- read.csv(paste0(last_shiny, "KL_tree_results.with.bp.40.csv"), row.names = 1, check.names = FALSE)
load(paste0(last_shiny, "mapping.df.with.bp.40.rda"))

#Patchseq Cells of interests
locked_cells_spec_id = rownames(read.csv(paste0(Fahimehdir, "2018_mouse_met_dataset.csv"), check.names=FALSE, row.names = 1 ))
locked_cells_sample_id = rownames(mapping.df[mapping.df$cell_id %in% locked_cells_spec_id,])

Tree_t_type_layer_info <- cbind(KL_tree_results[locked_cells_sample_id, 
                                                 c("Tree_first_cl", "Tree_second_cl", "Tree_third_cl", 
                                                   "Tree_first_bt", "Tree_second_bt", "Tree_third_bt", 
                                                   "Tree_call")], mapping.df[locked_cells_sample_id, "structure"])
colnames(Tree_t_type_layer_info) <- c("Tree_first_cl", "Tree_second_cl", "Tree_third_cl", 
                                      "Tree_first_bt", "Tree_second_bt", "Tree_third_bt", 
                                      "Tree_call", "structure")
write.csv(Tree_t_type_layer_info, paste0(Fahimehdir,"Shared_with_IVSCC/Tree_t_type_layer_info_05082019.csv"))
yours <- read.csv(paste0(last_shiny, "mapping.df.with.bp.40.lastmap.csv"), check.names = FALSE, row.names = 1)
