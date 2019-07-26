library(dendextend)
library(matrixStats)

resolve_cl <- function(cl.g, cl.med, markers, dat, map.dat, select.cells, p=0.7, low.th=0.2)
  {
    ##
    genes=names(markers)[markers > 0]
    tmp.cl= unlist(cl.g)
    #print(genes)
    #print(cl.g)
    ###For each branch point, find the highest expression cluster.
    #save(markers, file="markers.rda")
    #save(genes, file = "genes.rda")
    #save(tmp.cl, file = "tmp.cl.rda")
    tmp.med = sapply(cl.g, function(g)rowMaxs(cl.med[genes, g,drop=F]))
    #save(tmp.med, file = "tmp.med.rda")
    
    #print(tmp.med)
    row.names(tmp.med)= genes
    #print(tmp.med)
    ###Make sure the genes are discriminative between all the branches. 
    #print(rowMaxs(tmp.med)- rowMins(tmp.med))
    genes = genes[rowMaxs(tmp.med)- rowMins(tmp.med)>1]
    #print(genes)
    ###Sample the markers based on the weigts.
    ##TO DO: randomforest sometimes give importance value of 0. adjust for that.
    #print(length(markers))
    genes= sample(genes, round(length(genes) * p), prob=markers[genes])
    #print(c("genes", length(genes)))
    
    ###Compute the correlation with the median cluster profile.
    ###add drop=F
    cl.cor = cor(as.matrix(map.dat[genes,select.cells,drop=F]), cl.med[genes,tmp.cl,drop=F])
    #print(c("cl.cor", length(cl.cor)))
    cl.cor[is.na(cl.cor)]=0
    ###Compute the best match in each branch. 
    tmp.score = do.call("cbind",sapply(cl.g, function(x)rowMaxs(cl.cor[,x,drop=F]),simplify=F))
    row.names(tmp.score)= row.names(cl.cor)
    ####Determine the best match. 
    best.score=setNames(rowMaxs(tmp.score), row.names(tmp.score))
    ###determine the difference from the best match. 
    diff.score = best.score - tmp.score
    
    ####Give up on cells can't be discriminated,choose one branch randomly.
    unresolved.cl = row.names(tmp.score)[rowSums(diff.score < low.th)==ncol(diff.score)]
    #print(c("unresolved.cl", length(unresolved.cl)))
    mapped.cl= setNames(sample(colnames(tmp.score), length(unresolved.cl),replace=T), unresolved.cl)
    #print(c("mapped.cl", length(mapped.cl)))
    ###Cells mapped to one or more branches. 
    mapped.cells= setdiff(row.names(cl.cor),unresolved.cl)
    ###For binary branch, done already
    if(length(cl.g) == 2){
      mapped.cl= c(mapped.cl,setNames(colnames(diff.score)[apply(diff.score[mapped.cells,,drop=F], 1, which.min)], mapped.cells))
      return(mapped.cl)
    }
    ##The remaining options for mapped cells 
    tmp.cl= sapply(mapped.cells, function(x) colnames(diff.score)[which(diff.score[x,] < low.th)],simplify=F)
    ###cells with multiple options
    resolve.cells = names(tmp.cl)[sapply(tmp.cl, length)> 1]
    ###cells with only one option. Not further job. 
    mapped.cells = setdiff(mapped.cells, resolve.cells)
    if(length(mapped.cells)>0){
      mapped.cl = c(mapped.cl, setNames(unlist(tmp.cl[mapped.cells]), mapped.cells))
    }
    ###Resolve further options. 
    if(length(resolve.cells)>0){
      tmp.cat = sapply(tmp.cl[resolve.cells], function(x)paste(x, collapse=" "))
      for(cat in unique(tmp.cat)){
        tmp.cl= unlist(strsplit(cat, " "))
        select.cells = names(tmp.cat)[tmp.cat==cat]
        mapped.cl=c(mapped.cl,resolve_cl(cl.g[tmp.cl], cl.med, markers, dat, map.dat, select.cells, p=p,low.th=low.th))
      }
    }
    return(mapped.cl)
  }


map_dend <- function(dend,cl, dat, map.dat, select.cells, p=0.8, low.th=0.2, default.markers=NULL, cl.mean= NULL)
  {
    final.cl= c(setNames(rep(attr(dend,"label"),length(select.cells)),select.cells))
    #print(attr(dend,"label"))
    #print(attr(dend, "markers"))
    if(length(dend)<=1){
      return(final.cl)
    }
    #print(dend)
    markers = attr(dend, "markers")
    #print(attr(dend, "markers"))
    markers= markers[names(markers)%in% row.names(map.dat)]
    cl.g = sapply(dend,labels, simplify=F)
    #print(cl.g)
    names(cl.g)=1:length(cl.g)
    select.cl = cl[cl %in% unlist(cl.g)]
    ###Sampling the cells from the reference cluster
    cells = unlist(tapply(names(select.cl), select.cl, function(x)sample(x, round(length(x)*p))))
    genes=names(markers)
    genes=union(genes, default.markers)
    # ###Compute reference cluster median based on subsampled cells
    if(is.null(cl.mean)){
      print("This should not happen!")
      cl.med = do.call("cbind",tapply(cells, droplevels(cl[cells]),function(x)rowMeans(dat[genes, x,drop=F])))
      row.names(cl.med)= genes
    }else{
      print("Reading in cl means as input! to save time")
      cl.med <- cl.mean
    }
    ###determine which branch to take.
    #save(cl.g, file ="cl.g.rda")
    #save(markers, file ="markers.rda")
    mapped.cl = resolve_cl(cl.g, cl.med, markers, dat, map.dat, select.cells, p=p, low.th=low.th)
    if(length(mapped.cl)>0){
      for(i in unique(mapped.cl)){
        select.cells= names(mapped.cl)[mapped.cl==i]
        if(length(select.cells)>0){
          final.cl=c(final.cl,map_dend(dend[[as.integer(i)]],cl, dat, map.dat, select.cells, p=p,low.th=low.th, cl.mean = cl.mean))
        }
      }
    }
    return(cl=final.cl)
  }



get.markers.num <- function(markers.cl.list, map.dat)
  {
    all.markers=unique(unlist(markers.cl.list))
    markers.num = sapply(markers.cl.list, function(x){
      colSums(map.dat[x, ]>0.5)
    })
  }
 
summarize_cl <- function(dend, memb, map.dat, exp.th = 1, conf.th=0.7, min.genes.ratio = 0.3,min.genes=3)
  {
    node.height=setNames(get_nodes_attr(dend, "height"),get_nodes_attr(dend, "label"))
    dend.list = dend_list(dend)
    tmp = sapply(dend.list, length)>1
    tmp.dend.list = dend.list[tmp]
    markers.cl.list = lapply(tmp.dend.list, function(x){
      tmp = attr(x,"markers.byCl")
	  if(!is.null(tmp))  # So it doesn't crash on the leaves
        names(tmp) = sapply(x, function(y)attr(y, "label"))
      tmp
    })
    names(markers.cl.list)=NULL
    markers.cl.list=do.call("c",markers.cl.list)

    all.markers=unique(unlist(markers.cl.list))
    memb.th= lapply(row.names(memb),function(cell){
      ###Check all the node with confidence > conf.th
      x = memb[cell,]
      mapped.node = colnames(memb)[which(x>conf.th)]
      
      ###Check for detected markers at the given cell
      det.genes = all.markers[map.dat[all.markers,cell]>= exp.th]
      
      ##compute detected markers at every branch point.   
      gene.olap = sapply(mapped.node, function(i)intersect(markers.cl.list[[i]], det.genes),simplify=F)
      gene.olap.num =sapply(gene.olap, length)
      #TO DO: weight markers instead of using absolute counts/ratio. 

      ####set the root, so that root always succeed. 
      gene.olap.num[attr(dend,"label")] = min.genes
      gene.olap[attr(dend,"label")] = ""
      gene.olap.ratio =gene.olap.num/sapply(markers.cl.list[names(gene.olap.num)],length)
      gene.olap.ratio[is.na(gene.olap.ratio)] = 1
      
      ###mapped nodes not met the minimal gene number/ratio constraints 
      fail.node= mapped.node[gene.olap.ratio < min.genes.ratio | gene.olap.num < min.genes]
      if(length(fail.node)>0){
        ###choose the mapped nodes above any failed nodes  
        mapped.node = mapped.node[node.height[mapped.node] > max(node.height[fail.node])]
      }
      ###Choose the deepest nodes that pass all the criteria. 
      mapped.node=mapped.node[order(node.height[mapped.node])]
      best.node = mapped.node[1]
      ###Get the markers on every mapped nodes. 
      gene.anno = sapply(mapped.node, function(x){paste0(x,":",paste0(gene.olap[[x]], collapse=" "))})
      c(cl=best.node, score=x[best.node], marker.num=gene.olap.num[best.node],markers=paste(gene.anno, collapse=","))
    })
    memb.th = do.call("rbind",memb.th)
    row.names(memb.th) = row.names(memb)
    memb.df = data.frame(cl=memb.th[,1], score=as.numeric(memb.th[,2]), marker.num = as.integer(memb.th[,3]), stringsAsFactors=F)
    memb.df$resolution.index = 1- (node.height[memb.df$cl]/attr(dend,"height"))
    memb.df$cl = factor(memb.df$cl, names(node.height))
    #tmp = t(t(memb) * (1-node.height[colnames(memb)]/max(node.height)))
    #memb.df$h.score = rowMaxs(tmp) 
    memb.df$h.score = memb.df$resolution.index * memb.df$score
    memb.df$markers = memb.th[,4]
    ord = order(-memb.df$resolution.index, memb.df$cl,-memb.df$h.score)
    memb.df = memb.df[ord,]
    
    # adding resolution.index.percentile
    qq <- sort(quantile(memb.df$resolution.index, probs = seq(0, 1, by= 0.01)))
    memb.df$resolution.index.percentile <- findInterval(memb.df$resolution.index, qq) - 1
    
    
    return(memb.df) 
  }

map_dend_membership <- function(dend, cl, dat, map.dat, map.cells, bs.num=100, cl.mean = NULL, ...)
  {
    mem=sapply(1:bs.num, function(i){
      print(i)
      tmp=map_dend(dend, cl, dat, map.dat, map.cells, cl.mean = cl.mean, ...)
    },simplify=F)
    memb = unlist(mem)
    memb = data.frame(cell=names(memb),cl=memb)
    memb = table(memb$cell,memb$cl)
    memb = memb/bs.num
    tmp=get_nodes_attr(dend, "label")
    tmp = tmp[tmp %in% colnames(memb)]
    memb = memb[,tmp]
    return(memb)
  }

plot_map_num <- function(dend, map.num,fn,height=4,width=7){
  label= dend %>% get_nodes_attr("label")
  size = sqrt(map.num)*0.5
  tmp.dend=dend %>% set("nodes_pch", 16) %>% set("nodes_cex", size[label]) %>% set("nodes_col","black")
  pdf(fn, height=height,width=width)
  plot(tmp.dend)
  dev.off()
}


plot_mapping <- function(dend, mapping_results)
  {
    require(dendextend)
    require(dplyr)
    require(ggplot2)
    dend_data <- as.ggdend(dend)
    dend_nodes <- dend_data$nodes %>% mutate(cl = get_nodes_attr(dend, "label"))
    node_counts <- mapping_results %>% group_by(cl) %>% summarise(n = n())
    dend_nodes <- dend_nodes %>% left_join(node_counts)
      
    mapping_plot <- ggplot() +
      geom_segment(data = dend_data$segments,
                   aes(x = x, xend = xend,
                       y = y, yend = yend),
                   color = "#808080") +
    geom_point(data = dend_nodes,
               aes(x = x,
                   y = y,
                   size = n),
               fill = "dodgerblue",
               color = "black",
               alpha = 0.7,
               pch = 21) +
                 scale_size_area() +
                   theme_void()
    mapping_plot
  }



####Work in progress. Not tested 
map_dend_membership_knn<- function(dend, cl, dat, map.dat, map.cells, bs.num=100, ...)
  {
    mem=sapply(1:bs.num, function(i){
      print(i)
      tmp=map_dend_knn(dend, cl, dat, map.dat, map.cells,...)
    },simplify=F)
    memb = unlist(mem)
    memb = data.frame(cell=names(memb),cl=memb)
    memb = table(memb$cell,memb$cl)
    memb = memb/bs.num
    tmp=get_nodes_attr(dend, "label")
    tmp = tmp[tmp %in% colnames(memb)]
    memb = memb[,tmp]
    return(memb)
  }



map_dend_knn <- function(dend,cl, dat, map.dat, select.cells, cell.p=0.7,gene.p=0.7, k=10)
  {
    require(FNN)
    final.cl= c(setNames(rep(attr(dend,"label"),length(select.cells)),select.cells))
    if(length(dend)<=1){
      return(final.cl)
    }
    markers = attr(dend, "markers")
    markers= markers[names(markers)%in% row.names(map.dat)]
    cl.g = sapply(dend,labels, simplify=F)
    names(cl.g)=1:length(cl.g)
    select.cl = cl[cl %in% unlist(cl.g)]
    cells = unlist(tapply(names(select.cl), select.cl, function(x)sample(x, round(length(x)*cell.p))))
    
    genes= sample(names(markers), round(length(markers)*gene.p))
    knn.result= knn(t(dat[genes,cells]), t(map.dat[genes,select.cells]), cl[cells],  k= k, prob=TRUE)
    prob = attributes(.Last.value)
    if(sum(rownames(cl.med) %in% genes) != length(genes)){
      print("cl.med should have rownames as genes when it is input!!!")
      row.names(cl.med)= genes
    }
    mapped.cl = resolve_cl(cl.g, cl.med, markers, dat, map.dat, select.cells, p=p,low.th=low.th)
    if(length(mapped.cl)>0){
      for(i in unique(mapped.cl)){
        select.cells= names(mapped.cl)[mapped.cl==i]
        if(length(select.cells)>0){
          final.cl=c(final.cl,map_dend(dend[[as.integer(i)]],cl, dat, map.dat, select.cells, p=p,low.th=low.th))
        }
      }
    }
    return(cl=final.cl)
  }



baysian_map_dend <- function(dend,cl, dat, map.dat, map.cells,pr, cl.med = NULL)
  {
    markers = attr(dend, "markers")
    markers= markers[names(markers)%in% row.names(map.dat)]
    cl.g = sapply(dend,labels, simplify=F)
    names(cl.g)=1:length(cl.g)
    select.cl = cl[cl %in% unlist(cl.g)]
    cells = names(select.cl)
    genes=names(markers)
    if(is.null(cl.med)){
      cl.med =  do.call("cbind",tapply(cells, droplevels(cl[cells]),function(x)rowMeans(dat[genes, x,drop=F])))
    }
    low.th = setNames(pmax(rep(1,length(genes)), rowMaxs(cl.med)-10), genes)
    cl.present = do.call("cbind",tapply(cells, droplevels(cl[cells]),function(x)(rowSums(dat[genes, x,drop=F] > low.th)+0.5)/(length(x)+1)))
    log.cl.present = log2(cl.present)
    log.cl.absent = log2(1-cl.present)
    cl.mix.present = rowMeans(sapply(cl.g, function(x)rowMeans(cl.present[,x,drop=F])))
    log.cl.mix.present = log2(cl.mix.present)
    log.cl.mix.absent = log2(1-cl.mix.present)
    
    logodds = t(apply(map.dat[genes,map.cells],2, function(x){
      pos = which(x > low.th)
      neg = which(x <= low.th)      
    }))
    logodds = logodds - rowMins(logodds)
    odds = 2^logodds
    odds = odds/rowSums(odds)
    branch.odds = sapply(cl.g, function(x)rowSums(odds[,x,drop=F]))
    branch.odds = cbind(branch.odds, mix= odds[,"mix"])
  }


