#' @include generics.R
#' @include Classes.R
NULL

#' @param object
#'
.separateBic <- function(object = NULL){
  tmp.expression <- object@Processed_count
  bic.number<-length(unique(object@BiCluster@CoCond_cell$Condition))
  Bic <- c()
  Bic.name <-c()
  for (i in 1:bic.number){
    gene.name <- object@BiCluster@CoReg_gene$Gene[object@BiCluster@CoReg_gene$Condition==i]
    cell.name <- object@BiCluster@CoCond_cell$cell_name[object@BiCluster@CoCond_cell$Condition==i]
    tmp.Bic <- tmp.expression[gene.name,cell.name]
    Bic.name <- c(Bic.name, paste0("Bicluster",i))
    Bic <- c(Bic, list(tmp.Bic))
  }
  names(Bic) <- Bic.name
  return(Bic)
}


#' Title PlotNetwork
#'
#' @description This function is for building the module- module network. Nodes mean individual module. Edges mean overlapped gene (or cell).
#' The weight of line show the number of gene (or cell).
#' @param object Input object
#' @param edge.by This parameter decide the edge by "cell" or by "gene". The default value is by gene.
#' @param lay.out The type of layout to create, include linear (default), circle, kk.
#' @param ...
#'
#' @name PlotNetwork
#' @return
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom ggraph ggraph geom_node_point geom_edge_arc scale_edge_width geom_node_text theme_graph
#' @import ggplot2
#' @examples
.plotnetwork <- function(object, edge.by = "gene",lay.out = "linear",N.bicluster =c(1:20), ...){
  Bic.list <- .separateBic(object)
  Bic.list.select<- Bic.list[N.bicluster]
  ntwork.adjacency.mtx <- matrix(1:(length(N.bicluster)*length(N.bicluster)),nrow = length(N.bicluster))
  rownames(ntwork.adjacency.mtx) <- colnames(ntwork.adjacency.mtx)<- names(Bic.list.select)
  for (i in 1:nrow(ntwork.adjacency.mtx)){
    for (j in 1:ncol(ntwork.adjacency.mtx)){
      tmp.block.1 <- Bic.list.select[[rownames(ntwork.adjacency.mtx)[i]]]
      tmp.block.2 <- Bic.list.select[[colnames(ntwork.adjacency.mtx)[j]]]
      if (edge.by == "gene"){
        n.intersect <- length(intersect(rownames(tmp.block.1),rownames(tmp.block.2)))
      } else if(edge.by == "cell"){
        n.intersect <- length(intersect(colnames(tmp.block.1),colnames(tmp.block.2)))
      } else {stop(paste("Please select 'gene' or 'cell' to edge.by parameter"))}
      ntwork.adjacency.mtx[i,j] <- n.intersect
    }
  }
  edge.list <- graph_from_adjacency_matrix(adjmatrix =ntwork.adjacency.mtx ,mode = "undirected",weighted = T,diag = F)
  label = rownames(ntwork.adjacency.mtx)
  p<- ggraph(edge.list, layout = lay.out) +
    geom_node_point() +
    geom_edge_arc(aes(width = weight), alpha = 0.8) +
    scale_edge_width(range = c(0.2, 2)) +
    geom_node_text(aes(label = label), repel = TRUE) +
    labs(edge_width = edge.by) +
    theme_graph()
  print(p)
}

#' @param BRIC
#'
#' @rdname PlotNetwork
#' @export
#'
setMethod("PlotNetwork","BRIC", .plotnetwork)


#' Title
#' @importFrom  stats cor
.generateNetObject  <- function(object, N.bicluster = c(1,5),method = "spearman"){
  groups <- N.bicluster
  tmp.data <- object@Processed_count
  if (length(N.bicluster) < 1) {stop("please choose number of biclusters")}
  rownamelist <- list()
  colnamelist <- list()
  for (i in seq_along(N.bicluster)){
    idx <- i
    rownamelist[[paste0("Bicluster_",N.bicluster[idx])]] <- object@BiCluster@CoReg_gene$Gene[object@BiCluster@CoReg_gene$Condition==N.bicluster[idx]]
    colnamelist[[paste0("Bicluster_",N.bicluster[idx])]] <- object@BiCluster@CoCond_cell$cell_name[object@BiCluster@CoCond_cell$Condition==N.bicluster[idx]]
  }
  allrownames <- Reduce(union, rownamelist)
  allcolnames <- Reduce(union, colnamelist)
  un <- tmp.data[allrownames, allcolnames]
  rowidlist <- list()
  index <- paste0("Bicluster_",N.bicluster)
  # Bic.list <- .separateBic(object)
  # bics <- Bic.list
  if (length(groups) > 2)
    stop("length(group) > 2")
  if (length(groups) == 1) {
    rowidlist[[index[[1]]]] <-
      match(rownamelist[[index[[1]]]],
            rownames(un))
  } else if (length(groups) == 2) {
    rowidlist[[paste0(index[[1]], " & ", index[[2]])]] <-
      match(intersect(rownamelist[[index[[1]]]], rownamelist[[index[[2]]]]),
            rownames(un))
    rowidlist[[index[[1]]]] <-
      match(setdiff(rownamelist[[index[[1]]]],
                    rownamelist[[index[[2]]]]), rownames(un))
    rowidlist[[index[[2]]]] <-
      match(setdiff(rownamelist[[index[[2]]]],
                    rownamelist[[index[[1]]]]), rownames(un))
    rowidlist[["Others"]] <-
      match(setdiff(allrownames, union(rownamelist[[index[[1]]]],
                                       rownamelist[[index[[2]]]])), rownames(un))
  }
  cort <- stats::cor(t(un), method = method)
  return(list(cort, rowidlist))
}

#' Title PlotModuleNetwork
#' @description This function will visualize
#' @param object Input object.
#' @param N.bicluster number of biclsuter to plot.
#' @param Node.color color of nodes. This parameter also accepts color codes, e.g. "#AE1503" or "darkred."
#' @importFrom igraph graph.adjacency degree
#' @import qgraph
#' @return
#' @name PlotModuleNetwork
#'
#' @examples
.plotmodulenetwork <- function(object = NULL, N.bicluster = c(1,5), Node.color = "#E8E504", node.label.cex = 1 ){
  my.list <- .generateNetObject(object =object,N.bicluster=N.bicluster )
  cort <- my.list[[1]]
  my.adjacency <- ifelse(abs(cort)<0.1,0,cort)
  g <- graph.adjacency(my.adjacency,weighted = T,diag = F,mode ="undirected" )
  degree.normalize <- 4*(degree(g)/max(degree(g)))
  a<- my.list
  if (length(N.bicluster)==1){
    qgraph(a[[1]], groups = a[[2]],
           theme = "classic",cut =0,
           layout = "circle", minimum = 0.5, posCol=c("grey"), negCol="darkred",
           # barLength = 0.5,
           legend.cex = 0.7, color = Node.color,vsize = degree.normalize,vsize2= degree.normalize,
           vTrans = 200,label.cex =node.label.cex,labels = rownames(a[[1]]),label.scale=F)
  }
  if (length(N.bicluster)==2){
    qgraph(a[[1]], groups = a[[2]],
           theme = "classic",cut =0,
           layout = "spring", minimum = 0.5, posCol="grey", negCol="darkred",
           # barLength = 0.5,
           legend.cex = 0.7, color = c("#AE1503", "#012290", "#908E03", "#808080"),vsize = degree.normalize,vsize2= degree.normalize, vTrans = 200,label.cex =node.label.cex,label.scale=F)
  }
}

#' @param BRIC
#'
#' @rdname PlotModuleNetwork
#' @export
#'
setMethod("PlotModuleNetwork","BRIC", .plotmodulenetwork)



