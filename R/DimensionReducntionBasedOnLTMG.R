#' @include generics.R
#' @include Classes.R
NULL

#' @param object
#'
#' @param reduction
#' @param dims
#' @name RundimensionReduction
#' @importFrom Seurat CreateSeuratObject ScaleData RunPCA RunTSNE RunUMAP FindVariableFeatures
.runDimensionReduction <- function(object, reduction = "tsne", dims = 1:15 ,perplexity = 15, seed = 1){
  Tmp.seurat <- CreateSeuratObject(object@LTMG@LTMG_discrete)
  Tmp.seurat<- ScaleData(Tmp.seurat)
  Tmp.seurat <- suppressMessages(RunPCA(Tmp.seurat, features = rownames(Tmp.seurat@assays$RNA)))
  object@LTMG@DimReduce@PCA <- Tmp.seurat@reductions$pca@cell.embeddings
  if(grepl("tsne", reduction, ignore.case = T) || grepl("umap", reduction, ignore.case = T)){
    if(grepl("tsne", reduction, ignore.case = T)){
      Tmp.seurat <- RunTSNE(Tmp.seurat,dims=dims , seed.use = seed,
                            perplexity = perplexity)
      object@LTMG@DimReduce@TSNE <- Tmp.seurat@reductions$tsne@cell.embeddings
    }
    if(grepl("umap", reduction, ignore.case = T)){
      Tmp.seurat <- suppressMessages(RunUMAP(Tmp.seurat,dims = dims))
      object@LTMG@DimReduce@UMAP <- Tmp.seurat@reductions$umap@cell.embeddings
    }
  } else {stop("choose a dimension reduction method between umap or tsne")}
object@LTMG@Tmp.seurat<-Tmp.seurat
return(object)
}


#' @rdname RunDimensionRecution
#' @export
setMethod("RunDimensionReduction", "BRIC", .runDimensionReduction)



#' @param object
#' @param k.param
#' @param resolution
#' @param algorithm
#' @param dims
#' @name RunClassification
#' @importFrom Seurat FindNeighbors FindClusters
#' @import ggplot2
.runClassification <- function(object,dims = 1:15, k.param = 20, resolution = 0.6, algorithm = 1 ){
  if ( is.null(object@LTMG@Tmp.seurat)) {stop("There is no temporary seurat obejct getting detected. \n Try to run RundimensionRuduction first.")}
  Tmp.seurat <- object@LTMG@Tmp.seurat
  Tmp.seurat <- FindNeighbors(Tmp.seurat,dims=dims, k.param = k.param)
  # Seurat::ElbowPlot(Tmp.seurat, ndims = length(dims)+5)
  Tmp.seurat <- FindClusters(Tmp.seurat, resolution = resolution, algorithm = algorithm)
  tmp.meta <- object@MetaInfo
  tmp.colname <- colnames(tmp.meta)
  res.index <- grep(paste0("res.",resolution),colnames(Tmp.seurat@meta.data))
  tmp.colname <- c(tmp.colname,paste0("Seurat",resolution) )
  tmp.meta <- cbind(tmp.meta, Tmp.seurat@meta.data[,res.index])
  colnames(tmp.meta) <- tmp.colname
  object@MetaInfo <- tmp.meta
  object@LTMG@Tmp.seurat <- Tmp.seurat
  return(object)
}

#' @rdname RunDimensionRecution
#' @export
setMethod("RunClassification", "BRIC", .runClassification)




#' Title Visualize dimension reduction results
#'
#' @param object Input Object
#' @param reduction Choose one of approaches for dimension reduction, including "pca", "tsne", "umap".
#' @param pt_size point size, default is 0.
#'
#' @return
#' @export
#' @name PlotDimension
#' @examples
.plotDimension <- function(object, reduction = "umap", pt_size = 1){

  if(grepl("tsne", reduction, ignore.case = T) || grepl("umap", reduction, ignore.case = T)||grepl("pca", reduction, ignore.case = T)){

    if(grepl("tsne", reduction, ignore.case = T)){
      tmp.plot.table <- object@LTMG@DimReduce@TSNE[,c(1,2)]
    }
      else if(grepl("umap", reduction, ignore.case = T)){
      tmp.plot.table <- object@LTMG@DimReduce@UMAP[,c(1,2)]
     }
      else if(grepl("tsne", reduction, ignore.case = T)){
      tmp.plot.table <- object@LTMG@DimReduce@PCA[,c(1,2)]
    }
  } else {stop("choose a dimension reduction method from pca, umap, or tsne")}
  message("select condition to present")
  message(paste0(c(1:ncol(object@MetaInfo))," : ",c(colnames(object@MetaInfo)),"\n"))
  ident.index <-  readline(prompt="select index of cell condition: ")
  ident.index <- as.numeric(ident.index)
  tmp.ident <- object@MetaInfo[,ident.index]
  names(tmp.ident) <- rownames(object@MetaInfo)
  # check name later add
  tmp.plot.table <- cbind.data.frame(tmp.plot.table,Cell_type=as.character(tmp.ident))
  p.cluster <- ggplot(tmp.plot.table, aes(x= tmp.plot.table[,1],y = tmp.plot.table[,2]))

  p.cluster <- p.cluster+geom_point(stroke=pt_size,size=pt_size,aes(col=tmp.plot.table[,"Cell_type"]))

  p.cluster <- p.cluster + guides(colour = guide_legend(override.aes = list(size=5))) + labs( color ="cell type")+ xlab("Dimension 1") + ylab("Dimentsion 2")
  p.cluster <- p.cluster+theme_classic()
  print(p.cluster)
}

#' @rdname PlotDimension
#' @export
setMethod("PlotDimension", "BRIC", .plotDimension)





