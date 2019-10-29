get_node_dend <- function(x, match_attr, match_value) {
  
  output <- NULL
  
  for(i in seq_len(length(x))) {
    
    if(attr(x[[i]], match_attr) == match_value) {
      
      output <- x[[i]]
      
    } else {
      if(is.list(x[[i]])) {
        nest <- get_node_dend(x[[i]], match_attr, match_value)
        if(!is.null(nest)) {
          output <- nest
        }
      }
    }
    
  }
  return(output)
  
}

build_layer_plot <- function(anno,
                             dend,
                             dendcluster_ids,
                             seed_val = 42,
                             right_pad = 10,
                             modify_layer_label = FALSE,
                             patch= FALSE) {

  
  cluster_anno <- anno %>%
      select(dendcluster_id, cluster_id, cluster_label, cluster_color) %>%
      unique()
  
  if (patch) {
    keep_layers <- c("1", "2/3", "4", "5","6a" ,"6b")
  } else {
    keep_layers <- c("L1","L2/3","L4","L5","L6")
    #keep_layers <- c("L2/3","L4","L5","L6","L6b")
  }
  
  
  xpad <- 0.1
  ypad <- 0.05
  

  if (patch) {
    if(modify_layer_label){
      anno$layer_label <- anno$Revisited_layer_label
      anno$layer_id <- anno$Revisited_layer_id
      anno$layer_color <- anno$Revisited_layer_color
      #keep_layers <- c("VISp1", "VIS1_ho", "VISp2/3", "VIS23_ho", "VISp4", "VIS4_ho", 
      #                 "VISp5", "VIS5_ho", "VISp6a", "VIS6a_ho", "VISp6b", "VIS6b_ho")
      keep_layers <- c("VIS1_ho", "VIS23_ho", "VIS4_ho", "VIS5_ho", "VIS6a_ho", "VIS6b_ho")
      #keep_layers <- c("VISp1", "VISp2/3", "VISp4",  "VISp5",  "VISp6a",  "VISp6b")
    }
  }
  
  layer_ranges <- data.frame(layer_label = rev(keep_layers),
                             ymin = (1:6 - 1) + ypad,
                             ymax = (1:6) - ypad)
  
  filtered_anno <- anno %>%
      filter(dendcluster_id %in% dendcluster_ids) %>%
      filter(Revisited_layer_label %in% keep_layers)
  

  cluster_ranges <- filtered_anno %>%
    select(cluster_id, cluster_color, cluster_label, dendcluster_id) %>%
    unique() %>%
    arrange(dendcluster_id) %>%
    mutate(xmin = 1:n() - 1 + xpad,
           xmax = 1:n()     - xpad,
           xmid = 1:n() - 0.5)
  
  set.seed(seed_val)
  

  plot_data <- filtered_anno %>%
      select(sample_id,
             dendcluster_id, cluster_color, cluster_label,
             layer_id, layer_color, layer_label) %>%
      left_join(layer_ranges) %>%
      left_join(cluster_ranges) %>%
      group_by(dendcluster_id, layer_id) %>%
      mutate(x = runif(n(),xmin + xpad, xmax - xpad),
             y = runif(n(),ymin + ypad, ymax - ypad),
             fill_color = cluster_color)

  
  # Layer background rectangles
  layer_rects <- layer_ranges %>%
    mutate(xmin = min(cluster_ranges$xmin) - xpad, xmax = max(cluster_ranges$xmax) + xpad) %>%
    mutate(fill = c( "#ECE09C","#F7F2DA","#A7D7DF","#C1E5E7","#D4EDED", "#FDE4DF" ))
    #mutate(fill = c("#FDE4DF","#ECE09C","#F7F2DA","#A7D7DF","#C1E5E7"))
    
  # Cluster color highlights at bottom of the plot
  cluster_rects <- cluster_ranges %>%
    mutate(ymin = -ypad, ymax = ypad)
  
  # Filter the dendrogram
  prune_dend_labels <- labels(dend)[!labels(dend) %in% filtered_anno$cluster_label]
  if (patch) {prune_dend_labels <- labels(dend)[!labels(dend) %in% filtered_anno$topLeaf_label]}
  filtered_dend <- dend %>%
    prune.dendrogram(prune_dend_labels)
  dend_seg <- as.ggdend(filtered_dend)$segments %>%
    mutate(y = (y/max(y))*3 + max(layer_rects$ymax) + ypad,
           yend = (yend/max(yend))*3 + max(layer_rects$ymax) + ypad,
           x = x - 0.5,
           xend = xend - 0.5)
  
  pad_rect <- data.frame(ymin = min(layer_rects$ymin),
                         ymax = max(layer_rects$ymax),
                         xmin = max(layer_rects$xmax),
                         xmax = max(layer_rects$xmax) + max(layer_rects$xmax)*(right_pad/100)/(1 - right_pad/100))
  
  p <- ggplot() +
    # right side padding for alignment
    geom_rect(data = pad_rect,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = "#FFFFFF",
                  color = "#FFFFFF")) +
    # dendrogram segments
    geom_segment(data = dend_seg,
                 aes(x = x, xend = xend,
                     y = y, yend = yend,
                     size = lwd,
                     color = col),
                 lineend = "square") +
    # layer background rectangles
    geom_rect(data = layer_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fill)) +
    # cluster label rectangles
    geom_rect(data = cluster_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = cluster_color)) +
    geom_rect(data = cluster_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = -2 - ypad, ymax = -2,
                  fill = cluster_color)) +
    # jittered cell points
    geom_point(data = plot_data,
               aes(x = x,
                   y = y,
                   color = cluster_color),
               size = 0.1) +
    # cluster name labels
    geom_rect(data = cluster_ranges,
              aes(xmin = xmid - 0.5 + xpad/2,
                  xmax = xmid + 0.5 - xpad/2,
                  ymax = 0 - ypad,
                  ymin = -2),
              fill = "#CAD7D7")+
    geom_text(data = cluster_ranges,
              aes(x = xmid,
                  y = -2 + ypad,
                  label = cluster_label),
              angle = 90,
              vjust = 0.3,
              hjust = 0,
              size = 2.5) +
    # Plot settings
    scale_color_identity() +
    scale_size(range = c(0.5,1), guide = FALSE) +
    scale_fill_identity() +
    scale_y_continuous(limits = c(-2.1, 8)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_void(base_size = 7) +
    theme(text = element_text(size = 6),
          legend.box.spacing = unit(0,"pt"))
  
  p
}

group_violin_plot2 <- function (genes = c("Hspa8", "Snap25", "Gad2", "Vip"), group_by = "final", 
                                clusters = 1:10, data_source = "internal", sort = F, logscale = F, 
                                showcounts = T, rotatecounts = F, fontsize = 7, labelheight = 25, 
                                max_width = 10, anno.feather = NULL, data.feather = NULL, dend = dend) 
{
  library(dplyr)
  library(ggplot2)
  genes <- rev(genes)

  if (is.null(anno.feather) & is.null(data.feather)) {
    data_file <- paste0(data_source, "/data.feather")
    anno_file <- paste0(data_source, "/anno.feather")
    data <- feather::feather(data_file)
    anno <- feather::read_feather(anno_file) %>% dplyr::mutate_if(is.factor, 
                                                                  as.character)
    anno <- anno[anno$cluster_label %in% labels(dend) & anno$region_label == "VISp",]
    data <- get_feather_data(anno, data, genes, group_by, 
                               clusters)
  }else {
    anno <- anno.feather %>% dplyr::mutate_if(is.factor, as.character)
    data <- get_feather_data(anno = anno, data = data.feather, genes = genes, group_by = group_by,
                             group_ids =  clusters) 
                             
  }
  genes <- sub("-", ".", genes)
  genes <- genes[genes %in% names(data)]
  data <- data %>% select(-xpos) %>% mutate(xpos = plot_id)
  genes[grepl("^[0-9]", genes)] <- paste0("X", genes[grepl("^[0-9]", 
                                                           genes)])
  names(data)[grepl("^[0-9]", genes)] <- paste0("X", names(data)[grepl("^[0-9]", 
                                                                       genes)])
  
  ngenes <- length(genes)
  nclust <- length(unique(data$plot_id))
  max_vals <- data %>% select(one_of(genes)) %>% summarise_each(funs(max)) %>% 
    unlist()
  data[genes] <- data[genes] + runif(nrow(data), 0, 1e-05)
  for (i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- max_vals[i]
    if (logscale) {
      data[gene] <- log10(data[gene] + 1)/log10(gene_max + 
                                                  1) * 0.9 + i
    }
    else {
      data[gene] <- data[gene]/gene_max * 0.9 + i
    }
  }
  print(names(data))
  header_labels <- build_header_labels(data = data, grouping = "plot", 
                                       ymin = ngenes + 1,
                                       label_height = labelheight, 
                                       label_type = "simple")
  max_labels <- data.frame(x = (nclust + 0.5) * 1.01, y = 1:ngenes + 
                             0.5, label = sci_label(max_vals))
  max_header <- data.frame(x = (nclust + 0.5) * 1.01, y = ngenes + 
                             1, label = "Max value")
  max_width <- nclust * (max_width/100)/(1 - max_width/100)
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  cluster_data <- data %>% group_by(plot_label, plot_color, 
                                    plot_id) %>% summarise(cn = n()) %>% as.data.frame(stringsAsFactors = F) %>% 
    arrange(plot_id) %>% mutate(labely = ngenes + label_y_size * 
                                  0.05, cny = max(header_labels$ymax) - 0.1 * label_y_size, 
                                xpos = plot_id)
  background_rects <- data.frame(ymin = 1:length(genes),
                                 ymax = 1:length(genes) + 1,
                                 xmin = 0.5,
                                 xmax = nclust + 0.5) %>%
    mutate(fill = ifelse(row_number() %% 2 == 0, "#CAD7D7", "#FFFFFF"))
  
  p <- ggplot() + 
    scale_fill_identity() + 
    scale_y_continuous("", 
                       breaks = 1:length(genes) + 0.45, 
                       labels = genes, 
                       expand = c(0,0)) + 
    scale_x_continuous("", expand = c(0, 0)) + 
    theme_classic(fontsize) + 
    theme(axis.text = element_text(size = rel(1)), 
          axis.text.x = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.y = element_text(face = "italic", color = "#000000"),
          legend.position = "none",
          axis.line = element_line(size = 0.1)) +
    geom_rect(data = background_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fill))
  
  for (i in 1:length(genes)) {
    
    p <- p + geom_violin(data = data, 
                         aes_string(x = "xpos", 
                                    y = genes[i], 
                                    #fill = "plot_color"),
                                    fill = as.factor(data$xpos)),
                         fill = "#000000",
                         scale = "width", 
                         size = 0.1,
                         adjust = 2) +
      stat_summary(data = data, aes_string(x = "xpos", 
                                           y = genes[i]), 
                   fun.y = "median", 
                   fun.ymin = "median", 
                   fun.ymax = "median", 
                   geom = "point", 
                   size = 0.1)
    
    
  }
  p <- p + geom_rect(data = header_labels, 
                     aes(xmin = xmin, 
                         ymin = ymin, 
                         xmax = xmax, 
                         ymax = ymax, 
                         fill = color)) + 
    geom_rect(aes(xmin = nclust +  0.5, 
                  xmax = nclust + 0.5 + max_width, 
                  ymin = 1, 
                  ymax = max(header_labels$ymax)), 
              fill = "white") + 
    geom_text(data = max_labels, 
              aes(x = x, 
                  y = y, 
                  label = label), 
              hjust = 0, 
              vjust = 0.35, 
              size = pt2mm(fontsize), 
              parse = TRUE)
  
  if (showcounts) {
    if (rotatecounts) {
      p <- p + geom_text(data = cluster_data, aes(y = cny, 
                                                  x = xpos, label = cn), angle = 90, vjust = 0.35, 
                         hjust = 1, size = pt2mm(fontsize))
    }
    else {
      p <- p + geom_text(data = cluster_data, aes(y = cny, 
                                                  x = xpos, label = cn), size = pt2mm(fontsize))
    }
  }
  p
}



get_feather_data <- function (anno, data, genes, group_by, group_ids)#, dend=NULL) 
{
  id_cols <- names(anno)[grepl("_id$", names(anno)) & names(anno) != 
                           "sample_id"]
  anno[id_cols] <- lapply(anno[id_cols], as.numeric)
  data_names <- names(data)
  if (sum(genes %in% data_names) != length(genes)) {
    not_found <- genes[!toupper(genes) %in% toupper(data_names)]
    warning(paste(paste0(not_found, collapse = ", "), "not found in feather data!"))
    genes <- data_names[toupper(data_names) %in% toupper(genes)]
  }
  data_cols <- which(data_names %in% c("sample_id", genes))
  gene_data <- data[, data_cols]
  colnames(gene_data) <- gsub("-", ".", colnames(gene_data))
  genes <- gsub("-", ".", genes)
  all_anno <- anno %>% dplyr::rename_(plot_id = paste0(group_by, 
                                                       "_id"), plot_label = paste0(group_by, "_label"), plot_color = paste0(group_by, 
                                                                                                                            "_color"))
  cluster_order <- data.frame(group_ids = group_ids) %>% dplyr::mutate(cluster_x = 1:n())
  data <- dplyr::left_join(all_anno, gene_data, by = "sample_id") %>% 
    dplyr::filter(plot_id %in% group_ids) %>% dplyr::left_join(cluster_order, 
                                                               by = c(plot_id = "group_ids")) %>% dplyr::arrange(cluster_x) %>% 
    dplyr::mutate(xpos = 1:n()) %>% dplyr::select(-plot_id) %>% 
    dplyr::rename_(plot_id = "cluster_x")
  return(data)
}


Assign_clusterlabel_clusterid_for_patchseq <- function(patchseq_dir, dend, ref_cluster_label_color){
  
  anno <- read_feather(paste0 (patchseq_dir , "/anno.feather"))
  anno <- as.data.frame(anno)
  dend <- dend
  ref_cluster_label_color <- ref_cluster_label_color

  #modifying cluster label and id and color
  anno$cluster_label <- anno$topLeaf_label
  anno$dendcluster_label <- anno$cluster_label
  
  anno <- anno[,setdiff(colnames(anno), "cluster_color")]
  anno  <- inner_join(anno, ref_cluster_label_color[,c("cluster_color","cluster_label")], by = "cluster_label")
  anno$dendcluster_color <- anno$cluster_color
  
  anno <- anno[,setdiff(colnames(anno), "cluster_id")]
  anno  <- inner_join(anno, ref_cluster_label_color[,c("cluster_label","cluster_id")], by = "cluster_label")
  anno  <- inner_join(anno, ref_cluster_label_color[,c("dendcluster_id","dendcluster_label")], by = "dendcluster_label")
  anno
}


Modify_layer_label_GABAcells <- function(GABAanno) {
  L1 <- c("VISp1")
  all_L1 <- unique(GABAanno$structure_label)[grep(c("VIS.*1"), unique(GABAanno$structure_label))]
  L1_ho <- setdiff(all_L1, L1)
  L23 <- c("VISp2/3")
  all_L23 <-  unique(GABAanno$structure_label)[grep(c("VIS.*2/3"), unique(GABAanno$structure_label))]
  L23_ho <- setdiff(all_L23, L23)
  L4 <- c("VISp4")
  all_L4 <- unique(GABAanno$structure_label)[grep(c("VIS.*4"), unique(GABAanno$structure_label))]
  L4_ho <- setdiff(all_L4, L4)
  L5 <- c("VISp5")
  all_L5 <- unique(GABAanno$structure_label)[grep(c("VIS.*5"), unique(GABAanno$structure_label))]
  L5_ho <- setdiff(all_L5, L5)
  L6a <- c("VISp6a")
  all_L6a <- unique(GABAanno$structure_label)[grep(c("VIS.*6a"), unique(GABAanno$structure_label))]
  L6a_ho <- setdiff(all_L6a, L6a) 
  L6b <- c("VISp6b")
  all_L6b <- c(unique(GABAanno$structure_label)[grep(c("VIS.*6b"), unique(GABAanno$structure_label))])
  L6b_ho <- setdiff(all_L6b, L6b)
  
  GABAanno <- GABAanno[GABAanno$structure_label %in% c(L1, L23, L4, L5, L6a, L6b, L1_ho, L23_ho, L4_ho, L5_ho, L6a_ho, L6b_ho),] %>%
    rownames_to_column("id") %>%
    mutate(Revisited_layer_label = case_when(structure_label %in% L1 ~ "VISp1",
                                           structure_label %in% L23 ~ "VISp2/3",
                                           structure_label %in% L4 ~ "VISp4", 
                                           structure_label %in% L5 ~ "VISp5",
                                           structure_label %in% L6a ~ "VISp6a",
                                           structure_label %in% L6b ~ "VISp6b",
                                           structure_label %in% L1_ho ~ "VIS1_ho",
                                           structure_label %in% L23_ho ~ "VIS23_ho",
                                           structure_label %in% L4_ho ~ "VIS4_ho",
                                           structure_label %in% L5_ho ~ "VIS5_ho",
                                           structure_label %in% L6a_ho ~ "VIS6a_ho",
                                           structure_label %in% L6b_ho ~ "VIS6b_ho",
                                           TRUE ~ structure_label)) %>% 
    mutate(Revisited_layer_color = case_when(structure_label %in% L1 ~ "#3A1799",
                                           structure_label %in% L1_ho ~ "#941799",
                                           structure_label %in% L23 ~ "#5300FF", 
                                           structure_label %in% L23_ho ~ "#FF00FA",
                                           structure_label %in% L4 ~ "#875CCC",
                                           structure_label %in% L4_ho ~ "#CC5CC3",
                                           structure_label %in% L5 ~ "#5D2E99",
                                           structure_label %in% L5_ho ~ "#992E8B",
                                           structure_label %in% L6a ~ "#9326FF",
                                           structure_label %in% L6a_ho ~ "#FF26D5",
                                           structure_label %in% L6b ~ "#7200CC",
                                           structure_label %in% L6b_ho ~ "#FF73B3",
                                           TRUE ~ "grey")) %>%
    mutate(Revisited_layer_id = case_when(structure_label %in% L1 ~ 75,
                                        structure_label %in% L1_ho ~ 76,
                                        structure_label %in% L23 ~ 77, 
                                        structure_label %in% L23_ho ~ 78,
                                        structure_label %in% L4 ~ 79,
                                        structure_label %in% L4_ho ~ 80,
                                        structure_label %in% L5 ~ 81,
                                        structure_label %in% L5_ho ~ 82,
                                        structure_label %in% L6a ~ 83,
                                        structure_label %in% L6a_ho ~ 84,
                                        structure_label %in% L6b ~ 85,
                                        structure_label %in% L6b_ho ~ 86,
                                        TRUE ~ 90)) %>% 
    column_to_rownames("id")
  GABAanno 
}

build_dotplot_comparison_FACS_patch_plot <- function(facs.anno,
                                                      patch.anno,
                                                      dend,
                                                      dendcluster_ids,
                                                      seed_val = 42,
                                                      right_pad = 10){
  
  facs.anno <- as.data.frame(facs.anno)
  patch.anno <- as.data.frame(patch.anno)
  patch.ho.keep_layers <- c("VIS1_ho", "VIS23_ho", "VIS4_ho", "VIS5_ho", "VIS6_ho")
  patch.keep_layers <- c("VISp1", "VISp2/3", "VISp4", "VISp5","VISp6")
  facs.keep_layers <- c("L1","L2/3","L4","L5","L6")
  
  xpad <- 0.1
  ypad <- 0.05

  patch.anno$layer_label <- patch.anno$Revisited_layer_label
  patch.anno$layer_id <- patch.anno$Revisited_layer_id
  patch.anno$layer_color <- patch.anno$Revisited_layer_color
  
  facs.filtered_anno <- facs.anno %>%
    filter(dendcluster_id %in% dendcluster_ids) %>%
    filter(layer_label %in% facs.keep_layers)
  
  patch.filtered_anno <- patch.anno %>%
    filter(dendcluster_id %in% dendcluster_ids) %>%
    filter(layer_label %in% patch.keep_layers)
  
  #Layer range is the same for both patch and facs
  facs.layer_ranges <- data.frame(layer_label = rev(facs.keep_layers),
                                  ymin = seq(1, 3, by=0.5) -1 + ypad,
                                  ymax = seq(1, 3, by=0.5) -0.5 - ypad) %>% mutate(ymid = (ymin + ymax)/2) 
  
  patch.layer_ranges <- data.frame(layer_label = rev(patch.keep_layers),
                                   ymin = seq(1, 3, by=0.5) -1 + ypad,
                                   ymax = seq(1, 3, by=0.5) -0.5 - ypad) %>% mutate(ymid = (ymin + ymax)/2) 
  
  facs.cluster_ranges <- facs.filtered_anno %>%
    select(cluster_id, cluster_color, cluster_label, dendcluster_id) %>%
    unique() %>%
    arrange(dendcluster_id) %>%
    mutate(xmin = 1:n() - 1 + xpad,
           xmax = 1:n()     - xpad,
           xmid = 1:n() - 0.5) 
  

  set.seed(seed_val)
  
  dim(facs.filtered_anno)
  dim(patch.filtered_anno)
  
  facs.plot_data <- as.data.frame(facs.filtered_anno %>%
                                    select(sample_id,
                                           dendcluster_id, cluster_color, cluster_label,
                                           layer_id, layer_color, layer_label) %>%
                                    left_join(facs.layer_ranges) %>%
                                    left_join(facs.cluster_ranges) %>%
                                    group_by(dendcluster_id, layer_id, xmid, ymid) %>%
                                    mutate(ly_n = n()) %>%
                                    ungroup() %>%
                                    group_by(dendcluster_id) %>%
                                    arrange(layer_id) %>%
                                    mutate(cluster_n = n(),
                                           ly_frac = ly_n/cluster_n))

  patch.plot_data <- as.data.frame(patch.filtered_anno %>%
                                     select(sample_id,
                                            dendcluster_id, cluster_color, cluster_label,
                                            layer_id, layer_color, layer_label) %>%
                                     left_join(patch.layer_ranges) %>%
                                     left_join(facs.cluster_ranges) %>%
                                     group_by(dendcluster_id, layer_id, xmid, ymid) %>%
                                     mutate(ly_n = n()) %>%
                                     ungroup() %>%
                                     group_by(dendcluster_id) %>%
                                     arrange(layer_id) %>%
                                     mutate(cluster_n = n(),
                                            ly_frac = ly_n/cluster_n))  
  
  
  # Layer background rectangles
  #layer.rect.colors <- c("#ECE09C","#ECE09C", "#ECE09C","#ECE09C", "#ECE09C", "#ECE09C")
  #layer.rect.colors <- c("#ECE09C","#ECE09C", "#ECE09C","#ECE09C", "#ECE09C")

  facs.layer_rects <- data.frame()
  i = 1
  for (l in facs.keep_layers) {
    print(l)
    tmp <- facs.cluster_ranges %>% select(cluster_label,xmin, xmax) %>% 
      mutate(xmin = facs.cluster_ranges$xmin , 
             xmax = facs.cluster_ranges$xmax ,
             layer_label = l) %>% 
      left_join(facs.layer_ranges)
    tmp$fill <- rep(c("#E0E0E0", "#FFFFFF"), dim(facs.cluster_ranges)[1])[1:dim(facs.cluster_ranges)[1]]
    facs.layer_rects <- rbind.data.frame(facs.layer_rects, tmp)
    i <- i+1
  }
  
  #patch.layer_rects <- facs.layer_rects %>% mutate(xmin = xmin +1,
  #                                                 xmax = xmax +1)
  
  # Cluster color highlights at bottom of the plot
  cluster_rects <- facs.cluster_ranges %>%
    mutate(ymin = -ypad, ymax = ypad)
  
  pad_rect <- data.frame(ymin = min(facs.layer_rects$ymin),
                         ymax = max(facs.layer_rects$ymax),
                         xmin = max(facs.layer_rects$xmax),
                         xmax = max(facs.layer_rects$xmax) + max(facs.layer_rects$xmax)*(right_pad/100)/(1 - right_pad/100))
  
  p <- ggplot() +
    # right side padding for alignment
    geom_rect(data = pad_rect,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = "#FFFFFF",
                  color = "#FFFFFF")) +
    # layer background rectangles
    geom_rect(data = facs.layer_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fill)) +
    geom_point(data = facs.plot_data,
               aes(x = xmid - 0.2,
                   y = ymid,
                   color = "#0000CC",
                   size = ly_frac)) +
    geom_point(data = patch.plot_data,
               aes(x = xmid + 0.2,
                   y = ymid,
                   color = "#FF0000",
                   size = ly_frac)) +
    # cluster name labels
    geom_rect(data = facs.cluster_ranges,
              aes(xmin = xmid - 0.5 + xpad/2,
                  xmax = xmid + 0.5 - xpad/2 ,
                  ymax = 0 - ypad,
                  ymin = -2),
              fill = "#CCE5FF")+
    geom_text(data = facs.cluster_ranges,
              aes(x = xmid,
                  y = -2 + ypad,
                  label = cluster_label),
              angle = 90,
              vjust = 0.3,
              hjust = 0,
              size = 4) +
    # Plot settings
    scale_color_identity() +
    scale_size_area(max_size = 5,
                    breaks = c(1,10,50,100,200,500)) +
    #scale_size(range = c(0.5,1), guide = FALSE) +
    scale_fill_identity() +
    scale_y_continuous(limits = c(-2.1, 4)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_void(base_size = 7) +
    theme(text = element_text(size = 6),
          legend.box.spacing = unit(0,"pt"))
  
  p
}

build_layer_comparison_FACS_patch_plot <- function(facs.anno,
                                                   patch.anno,
                                                   dend,
                                                   dendcluster_ids,
                                                   seed_val = 42,
                                                   right_pad = 10) {
  facs.anno <- as.data.frame(facs.anno)
  patch.anno <- as.data.frame(patch.anno)
  #patch.ho.keep_layers <- c("VIS1_ho", "VIS23_ho", "VIS4_ho", "VIS5_ho", "VIS6a_ho", "VIS6b_ho")
  #patch.keep_layers <- c("VISp1", "VISp2/3", "VISp4", "VISp5","VISp6a" ,"VISp6b")
  #facs.keep_layers <- c("L1","L2/3","L4","L5","L6", "L6b")
  patch.ho.keep_layers <- c("VIS1_ho", "VIS23_ho", "VIS4_ho", "VIS5_ho", "VIS6_ho")
  patch.keep_layers <- c("VISp1", "VISp2/3", "VISp4", "VISp5","VISp6")
  facs.keep_layers <- c("L1","L2/3","L4","L5","L6")
  
  xpad <- 0.1
  ypad <- 0.05
  
  patch.anno$layer_label <- patch.anno$Revisited_layer_label
  patch.anno$layer_id <- patch.anno$Revisited_layer_id
  patch.anno$layer_color <- patch.anno$Revisited_layer_color
  
  facs.filtered_anno <- facs.anno %>%
    filter(dendcluster_id %in% dendcluster_ids) %>%
    filter(layer_label %in% facs.keep_layers)
  
  patch.filtered_anno <- patch.anno %>%
    filter(dendcluster_id %in% dendcluster_ids) %>%
    filter(layer_label %in% patch.keep_layers)
  
  #Layer range is the same for both patch and facs
  facs.layer_ranges <- data.frame(layer_label = rev(facs.keep_layers),
                             ymin = (1:5 - 1) + ypad,
                             ymax = (1:5) - ypad)
  
  patch.layer_ranges <- data.frame(layer_label = rev(patch.keep_layers),
                                  ymin = (1:5 - 1) + ypad,
                                  ymax = (1:5) - ypad)
  
  facs.cluster_ranges <- facs.filtered_anno %>%
    select(cluster_id, cluster_color, cluster_label, dendcluster_id) %>%
    unique() %>%
    arrange(dendcluster_id) %>%
    mutate(xmin = seq(1, 2 * n(), by=2) - 1 + xpad,
           xmax = seq(1, 2 * n(), by=2)     - xpad,
           xmid = seq(1, 2 * n(), by=2) - 0.5)
  
  patch.cluster_ranges <- patch.filtered_anno %>%
    select(cluster_id, cluster_color, cluster_label, dendcluster_id) %>%
    unique() %>%
    arrange(dendcluster_id) %>%
    mutate(xmin = seq(2, 2 * n(), by=2) - 1 + xpad,
           xmax = seq(2, 2 * n(), by=2)     - xpad,
           xmid = seq(2, 2 * n(), by=2) - 0.5)
  
  set.seed(seed_val)
  
  dim(facs.filtered_anno)
  dim(patch.filtered_anno)
  
  facs.plot_data <- as.data.frame(facs.filtered_anno %>%
    select(sample_id,
           dendcluster_id, cluster_color, cluster_label,
           layer_id, layer_color, layer_label) %>%
    left_join(facs.layer_ranges) %>%
    left_join(facs.cluster_ranges) %>%
    group_by(dendcluster_id, layer_id) %>%
    mutate(x = runif(n(),xmin + xpad, xmax - xpad),
           y = runif(n(),ymin + ypad, ymax - ypad),
           fill_color = cluster_color))
  
  patch.plot_data <- as.data.frame(patch.filtered_anno %>%
                                    select(sample_id,
                                           dendcluster_id, cluster_color, cluster_label,
                                           layer_id, layer_color, layer_label) %>%
                                    left_join(patch.layer_ranges) %>%
                                    left_join(patch.cluster_ranges) %>%
                                    group_by(dendcluster_id, layer_id) %>%
                                    mutate(x = runif(n(),xmin + xpad, xmax - xpad),
                                           y = runif(n(),ymin + ypad, ymax - ypad),
                                           fill_color = cluster_color))
  
  
  
  # Layer background rectangles
  #layer.rect.colors <- c("#C1E5E7","#C1E5E7", "#FDE4DF","#ECE09C", "#F7F2DA", "#A7D7DF")
  layer.rect.colors <- c("#ECE09C","#ECE09C", "#ECE09C","#ECE09C", "#ECE09C", "#ECE09C")
  layer.rect.colors <- c("#ECE09C","#ECE09C", "#ECE09C","#ECE09C", "#ECE09C")
  #layer.rect.colors <- c("#000000","#000000", "#000000","#ECE09C", "#ECE09C", "#ECE09C")
  
  facs.layer_rects <- data.frame()
  i = 1
  for (l in facs.keep_layers) {
    print(l)
    tmp <- facs.cluster_ranges %>% select(cluster_label,xmin, xmax) %>% 
      mutate(xmin = facs.cluster_ranges$xmin , 
             xmax = facs.cluster_ranges$xmax ,
             layer_label = l) %>% 
      left_join(facs.layer_ranges)
    tmp$fill <- rep(c("#E0E0E0", "#FFFFFF"), dim(facs.cluster_ranges)[1])[1:dim(facs.cluster_ranges)[1]]
    facs.layer_rects <- rbind.data.frame(facs.layer_rects, tmp)
    i <- i+1
  }
  
  patch.layer_rects <- facs.layer_rects %>% mutate(xmin = xmin +1,
                                                   xmax = xmax +1)
    
  # Cluster color highlights at bottom of the plot
  cluster_rects <- facs.cluster_ranges %>%
    mutate(ymin = -ypad, ymax = ypad)
  
  # Filter the dendrogram
  #prune_dend_labels <- labels(dend)[!labels(dend) %in% filtered_anno$cluster_label]
  #if (patch) {prune_dend_labels <- labels(dend)[!labels(dend) %in% filtered_anno$topLeaf_label]}
  #filtered_dend <- dend %>%
  #  prune.dendrogram(prune_dend_labels)
  #dend_seg <- as.ggdend(filtered_dend)$segments %>%
  #  mutate(y = (y/max(y))*3 + max(layer_rects$ymax) + ypad,
  #         yend = (yend/max(yend))*3 + max(layer_rects$ymax) + ypad,
  #         x = x - 0.5,
  #         xend = xend - 0.5)
  
  pad_rect <- data.frame(ymin = min(facs.layer_rects$ymin),
                         ymax = max(facs.layer_rects$ymax),
                         xmin = max(facs.layer_rects$xmax),
                         xmax = max(facs.layer_rects$xmax) + max(facs.layer_rects$xmax)*(right_pad/100)/(1 - right_pad/100))
  
  p <- ggplot() +
    # right side padding for alignment
    geom_rect(data = pad_rect,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = "#FFFFFF",
                  color = "#FFFFFF")) +
    # dendrogram segments
    #geom_segment(data = dend_seg,
    #             aes(x = x, xend = xend,
    #                 y = y, yend = yend,
    #                 size = lwd,
    #                 color = col),
    #             lineend = "square") +
    # layer background rectangles
    geom_rect(data = facs.layer_rects,
              aes(xmin = xmin, xmax = xmax,
                 ymin = ymin, ymax = ymax,
                 fill = fill)) +
    geom_rect(data = patch.layer_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fill))+
    # cluster label rectangles
    #geom_rect(data = cluster_rects,
    #          aes(xmin = xmin, xmax = xmax,
    #              ymin = ymin, ymax = ymax,
    #              fill = cluster_color)) +
    #geom_rect(data = cluster_rects,
    #          aes(xmin = xmin, xmax = xmax,
    #              ymin = -2 - ypad, ymax = -2,
    #              fill = cluster_color)) +
    # jittered cell points
    geom_point(data = facs.plot_data,
               aes(x = x,
                   y = y,
                   color = "#0000CC"),
               size = 0.1) +
    geom_point(data = patch.plot_data,
               aes(x = x,
                   y = y,
                   color = "#FF0000"),
               size = 0.1) +
    # cluster name labels
    geom_rect(data = facs.cluster_ranges,
              aes(xmin = xmid - 0.5 + xpad/2,
                  xmax = xmid + 0.5 - xpad/2 + 1,
                  ymax = 0 - ypad,
                  ymin = -2),
              fill = "#CCE5FF")+
    geom_text(data = facs.cluster_ranges,
              aes(x = xmid + 0.5,
                  y = -2 + ypad,
                  label = cluster_label),
              angle = 90,
              vjust = 0.3,
              hjust = 0,
              size = 2.5) +
    # Plot settings
    scale_color_identity() +
    scale_size(range = c(0.5,1), guide = FALSE) +
    scale_fill_identity() +
    scale_y_continuous(limits = c(-2.1, 8)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_void(base_size = 7) +
    theme(text = element_text(size = 6),
          legend.box.spacing = unit(0,"pt"))
  
  p
}

get_violin_data <- function(anno.feather = anno.feather, data.feather= data.feather, 
                            FACS.or.patch= "FACS", dend=dend, genes = genes, logscale = logscale, group_by= group_by, clusters = clusters){
  
  anno <- anno.feather %>% dplyr::mutate_if(is.factor, as.character)
  if (FACS.or.patch == "FACS") {
    anno <- anno[anno$cluster_label %in% labels(dend) & anno$region_label == "VISp",]
  }
  data <- get_feather_data(anno = anno, data = data.feather, genes = genes, group_by = group_by,
                           group_ids =  clusters) 
  
  genes <- sub("-", ".", genes)
  genes <- genes[genes %in% names(data)]
  if (FACS.or.patch == "FACS") {
    data <- data %>% select(-xpos) %>% mutate(xpos = plot_id + (plot_id -1 ))
  } else {
    data <- data %>% select(-xpos) %>% mutate(xpos = plot_id  * 2)
  }
  
  genes[grepl("^[0-9]", genes)] <- paste0("X", genes[grepl("^[0-9]", 
                                                           genes)])
  names(data)[grepl("^[0-9]", genes)] <- paste0("X", names(data)[grepl("^[0-9]", 
                                                                       genes)])
  
  max_vals <- data %>% select(one_of(genes)) %>% summarise_each(funs(max)) %>% 
    unlist()
  data[genes] <- data[genes] + runif(nrow(data), 0, 1e-05)
  for (i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- max_vals[i]
    if (logscale) {
      data[gene] <- log10(data[gene] + 1)/log10(gene_max + 
                                                  1) * 0.9 + i
    }
    else {
      data[gene] <- data[gene]/gene_max * 0.9 + i
    }
  }
  return(list(data = data,max_vals= max_vals, genes = genes))
}



group_violin_FACS_patch_plot  <- function (genes = c("Hspa8", "Snap25", "Gad2", "Vip"), group_by = "final", 
                                                              clusters = 1:10, data_source = "internal", sort = F, logscale = F, 
                                                              showcounts = T, rotatecounts = F, fontsize = 7, labelheight = 25, 
                                                              max_width = 10, REF.anno.feather = REF.anno.feather, REF.data.feather = REF.data.feather, 
                                                              map.anno.feather = map.anno.feather, map.data.feather = map.data.feather,dend = dend ) 
{
  library(dplyr)
  library(ggplot2)

  genes <- rev(genes)
  
  tmp <- get_violin_data(anno.feather = REF.anno.feather, data.feather = REF.data.feather, 
                          FACS.or.patch = "FACS", dend = dend, genes = genes, logscale = logscale, group_by = group_by, clusters = clusters)
  REF.data <- tmp$data
  REF.max_vals <- tmp$max_vals
  tmp <- get_violin_data(anno.feather = map.anno.feather, data.feather = map.data.feather, 
                              FACS.or.patch = "Patch", dend = dend, genes = genes, logscale = logscale, group_by=group_by, clusters= clusters)
  map.data <- tmp$data
  map.max_vals <- tmp$max_vals
  
  genes <- tmp$genes
  ngenes <- length(genes)
  nclust <- length(union(unique(REF.data$plot_id), unique(map.data$plot_id)))

  header_labels <- build_header_labels(data = REF.data, grouping = "plot", 
                                       ymin = ngenes + 1,
                                       label_height = labelheight, 
                                       label_type = "simple")
  
  REF.max_labels <- data.frame(x = (nclust + 0.5) * 1.01, y = 1:ngenes + 
                             0.5, label = sci_label(REF.max_vals))
  
  map.max_labels <- data.frame(x = (nclust + 0.5) * 1.01, y = 1:ngenes + 
                                 0.5, label = sci_label(map.max_vals))
  
  max_header <- data.frame(x = (nclust + 0.5) * 1.01, y = ngenes + 
                             1, label = "Max value")
  max_width <- nclust * (max_width/100)/(1 - max_width/100)
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  cluster_data <- REF.data %>% group_by(plot_label, plot_color, 
                                    plot_id) %>% summarise(cn = n()) %>% as.data.frame(stringsAsFactors = F) %>% 
    arrange(plot_id) %>% mutate(labely = ngenes + label_y_size * 
                                  0.05, cny = max(header_labels$ymax) - 0.1 * label_y_size, 
                                xpos = plot_id)
  
  vertical_background_rect <- header_labels %>% mutate(ymin = 1,
                                                       ymax = length(genes) + 1,
                                                       xmin = seq(0.5,nclust*4,2)[1:nclust],
                                                       xmax = seq(2.5,nclust*4,2)[1:nclust], 
                                                       color = rep(c("#E0E0E0", "#FFFFFF"), dim(header_labels)[1])[1:dim(header_labels)[1]])

  
p <- ggplot() +
scale_fill_identity() +
   scale_y_continuous("",
                      breaks = 1:length(genes) + 0.45,
                      labels = genes,
                      expand = c(0,0)) +
   scale_x_continuous("", expand = c(0, 0, 0, 5)) +
   theme_classic(fontsize) +
   theme(axis.text = element_text(size = rel(1)),
         axis.text.x = element_blank(),
         axis.ticks = element_blank(),
         axis.text.y = element_text(face = "italic", color = "#000000"),
         legend.position = "none",
         axis.line = element_line(size = 0.1)) +
  geom_rect(data = vertical_background_rect,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill = color))


  for (i in 1:length(genes)) {
    p <- p + geom_violin(data = REF.data, 
                         aes_string(x = "xpos", 
                                    y = genes[i],
                                    group = "xpos"),
                         fill = "#0000CC",
                         scale = "width",
                         width = 0.9, 
                         size = 0.05,
                         adjust = 2) +
      stat_summary(data = REF.data, aes_string(x = "xpos", 
                                           y = genes[i]), 
                   fun.y = "median", 
                   fun.ymin = "median", 
                   fun.ymax = "median", 
                   geom = "point", 
                   size = 0.1) +
      geom_violin(data = map.data, 
                  aes_string(x = "xpos", 
                             y = genes[i],
                             group = "xpos"),
                  fill = "#FF0000",
                  scale = "width",
                  width = 0.9, 
                  size = 0.05,
                  adjust = 2) +
      stat_summary(data = map.data, aes_string(x = "xpos", 
                                               y = genes[i]), 
                   fun.y = "median", 
                   fun.ymin = "median", 
                   fun.ymax = "median", 
                   geom = "point", 
                   size = 0.1) 
      
  }

    p <- p + geom_text(data = REF.max_labels, 
              aes(x = 2*x  , 
                  y = y, 
                  label = label ),
              color = "#0000CC",
              hjust = 0, 
              vjust = 0.65, 
              size = pt2mm(fontsize), 
              parse = TRUE)  +
      geom_text(data = map.max_labels, 
                aes(x = 2*x +5 , 
                    y = y, 
                    label = label ),
                color = "#FF0000",
                hjust = 0, 
                vjust = 0.65, 
                size = pt2mm(fontsize), 
                parse = TRUE)
   
   if (showcounts) {
   if (rotatecounts) {
     p <- p + geom_text(data = cluster_data, aes(y = cny,
                                                 x = xpos, label = cn), angle = 90, vjust = 0.35,
                        hjust = 1, size = pt2mm(fontsize))
   }
   else {
     p <- p + geom_text(data = cluster_data, aes(y = cny,
                                                 x = xpos, label = cn), size = pt2mm(fontsize))
   }
 }
p
}


