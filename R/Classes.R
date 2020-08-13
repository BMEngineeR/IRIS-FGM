#' Create DimReduce object
#'
#' @slot PCA ANY.
#' @slot UMAP ANY.
#' @slot TSNE ANY.
setClass("DimReduce",slots = c(
  PCA = "ANY",
  UMAP = "ANY",
  TSNE = "ANY"
)
)
# create Bicluster object
#'
#' @slot block ANY.
#' @slot CoReg_gene ANY.
#' @slot CoCond_cell ANY.
setClass("Bicluster", slots = c(
  CoReg_gene = "ANY",
  CoCond_cell = "ANY",
  MarkerGene = "ANY",
  PathwayFromMC = "ANY",
  PathwayFromModule = "ANY"
)
)
# create LTMG object
#'
#' @slot LTMG_discrete matrix.
#' @slot LTMG_BinarySingleSignal matrix.
#' @slot DimReduce DimReduce
#' @slot Cluster ANY
#' @slot MarkerGene ANY
#' @slot Pathway ANY
#' @slot LTMG_BinaryMultisignal matrix.
setClass("LTMGr", slots = c(
  LTMG_discrete = "matrix",
  LTMG_BinarySingleSignal = "matrix",
  LTMG_BinaryMultisignal = "matrix",
  DimReduce = "DimReduce",
  MarkerGene = "ANY",
  Pathway = "ANY",
  Tmp.seurat = "ANY"
)
)
# set BRIC class
#'
#' @slot raw_count matrix.
#' @slot MetaInfo ANY.
#' @slot processed_count ANY.
#' @slot Discretization matrix.
#' @slot LTMG LTMGr.
#' @slot Bicluster Bicluster.
#'
#' @return
#' @export
#' @name BRIC
#' @rdname BRIC
#' @exportClass BRIC
setClass("BRIC",
         slots=c(Raw_count = "ANY",
                 Processed_count = "ANY",
                 MetaInfo = "ANY",
                 Discretization = "matrix",
                 LTMG = "LTMGr",
                 BiCluster = "Bicluster"
         )

)
