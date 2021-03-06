% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/EnrichPathway.R
\name{RunPathway}
\alias{RunPathway}
\alias{.runPathway}
\alias{RunPathway,BRIC-method}
\title{Functional enrichment analysis}
\usage{
RunPathway(object, ...)

.runPathway(
  object = NULL,
  module.number = 1,
  selected.gene.cutoff = 0.05,
  species = "Human",
  database = "GO",
  genes.source = c("CTS", "Bicluster")
)

\S4method{RunPathway}{BRIC}(
  object = NULL,
  module.number = 1,
  selected.gene.cutoff = 0.05,
  species = "Human",
  database = "GO",
  genes.source = c("CTS", "Bicluster")
)
}
\arguments{
\item{object}{Input IRIS-FGM object}

\item{module.number}{Select the numebr of bicluster to perform this function.}

\item{selected.gene.cutoff}{Set up a statistical significance cutoff for all identified DEGs.}

\item{species}{You can choose either "Human" or "Mouse"}

\item{database}{You can choose either "GO" or "KEGG" database}

\item{genes.source}{You can choose a gene list source, either "CTS" or "Bicluster." "CTS" means from cell-type-specific DEGs, 
and "Bicluster means using gene module from the selected bicluster."}
}
\value{
It will reture a function enrichment analysis.
}
\description{
This function will perform enrichment analysis based on a gene module or identified differentially expressed genes (DEG).
This function is also depended on clusterProfiler, AnnotationDbi, org.Mm.eg.db, and org.Hs.eg.db package.
}
\examples{
# If you want to perform this function based on identified DEGs, you should use: 
\dontrun{object <- RunPathway(object = NULL,module.number = NULL, selected.gene.cutoff = 0.05,
species = "Human", database = "GO", genes.source = "CTS" }
# If you want to perform this function based on the gene module from an identified bicluster, you should use: 
\dontrun{object <- RunPathway(object = NULL,module.number = NULL, selected.gene.cutoff = 0.05,
species = "Human", database = "GO", genes.source = "Bicluster" }
}
