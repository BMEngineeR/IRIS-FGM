#' @include generics.R
#' @include Classes.R
NULL

.RunGO <- function(genes.use = NULL, species = "mouse"){
  if(grepl("mouse", species, ignore.case = T)){
    pathway<- invisible(suppressMessages(enrichGO(gene=genes.use,OrgDb = org.Mm.eg.db,
                      ont="ALL",keyType = "SYMBOL",pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)))
  } else if(grepl("human", species, ignore.case = T)){
    pathway<- invisible(suppressMessages(enrichGO(gene=genes.use,OrgDb = org.Hs.eg.db,
                      ont="ALL",keyType = "SYMBOL",pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)))
  }
  return(pathway)
}

.IDConvert <- function(genes.use = NULL,species = NULL){
  if(grepl("mouse", species, ignore.case = T)){
    gene.convert <- bitr(genes.use, fromType = "SYMBOL",
                         toType = c("ENSEMBL", "ENTREZID"),
                         OrgDb = org.Mm.eg.db)
  }
  else if(grepl("human", species, ignore.case = T)){
    gene.convert <- bitr(genes.use, fromType = "SYMBOL",
                         toType = c("ENSEMBL", "ENTREZID"),
                         OrgDb = org.Hs.eg.db)
  }
  return(gene.convert)
}

.RunKEGG <- function(genes.use = NULL,species = "mouse"){
  if(grepl("mouse", species, ignore.case = T)){
    pathway<-invisible(suppressMessages(enrichKEGG(gene=genes.use,
                        organism = "mmu",
                        keyType = "kegg",pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)))
  }
  else if(grepl("human", species, ignore.case = T)){
    pathway<-invisible(suppressMessages(enrichKEGG(gene=genes.use,
                        organism = "hsa",
                        keyType = "kegg",pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05)))
  }
  return(pathway)
}
#' @param object
#'
#' @param module.number
#' @param selected.gene.cutoff
#' @param species "Human" "Mouse"
#' @param database "GO" "KEGG"
#' @param genes.source CTS Bicluster
#'
#' @importFrom clusterProfiler enrichKEGG enrichGO bitr
#' @import org.Mm.eg.db org.Hs.eg.db
#' @importFrom AnnotationDbi select
#' @name RunPathway
.runPathway <- function(object = NULL,module.number = NULL, selected.gene.cutoff = 0.05,
                        species = "Human", database = "GO", genes.source = c("CTS","Bicluster")){
  if (genes.source == "CTS"){
    tmp.table<- as.data.frame(object@LTMG@MarkerGene)
    genes.use.LTMG <- rownames(tmp.table)[tmp.table$pvalue.adj.FDR < selected.gene.cutoff]
    if(database == "GO"){
      pathway <- .RunGO(genes.use = genes.use.LTMG , species = species)
    }else if(database == "KEGG"){
      ID.covert <- .IDConvert(genes.use.LTMG,species =species)
      pathway <- .RunKEGG(genes.use = ID.covert$ENTREZID , species = species)
    }
    object@LTMG@Pathway <- pathway@result
  } else if (genes.source == "Bicluster" ){
    block.number <- module.number
    genes.use.module <- object@BiCluster@CoReg_gene$Gene[object@BiCluster@CoReg_gene$Condition == block.number]
    # run on Bicluster marker gene
      if(is.null(object@BiCluster@MarkerGene)){message("There is no gene in MarkerGene slot.
                                                     \n Ignore pathway analysis based on marker gene derived from MC defined cell type. ")
        genes.use.MC <- NULL
      } else{
        tmp.table<- object@BiCluster@MarkerGene
        genes.use.MC <- rownames(tmp.table)[tmp.table$pvalue.adj.FDR < selected.gene.cutoff]
      }
    gene.list <-c(genes.use.module = list(genes.use.module), genes.use.MC = list(genes.use.MC))
    if(database == "GO"){
      pathway <- lapply(gene.list, function(x) .RunGO(genes.use = x , species = species))
    }else if(database == "KEGG"){
      ID.covert<- lapply(gene.list, function(x) .IDConvert(genes.use = x , species = species))
      pathway <- lapply(ID.covert, function(x) .RunKEGG(genes.use = x$ENTREZID , species = species))
    }
    object@BiCluster@PathwayFromModule <- pathway$genes.use.module@result
    if(is.null(genes.use.MC)){
      message("not find MC based marker genes, ignore pathway enrichment.......")
    } else{ object@BiCluster@PathwayFromMC <- pathway$genes.use.MC@result}
  }

return(object)

}

#' @rdname RunPathway
#' @export
setMethod("RunPathway","BRIC",.runPathway)

