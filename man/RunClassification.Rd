% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/DimensionReducntionBasedOnLTMG.R
\name{RunClassification}
\alias{RunClassification}
\alias{.runClassification}
\title{Classify cell type prediction}
\usage{
RunClassification(object, ...)

.runClassification(
  object,
  dims = 1:15,
  k.param = 20,
  resolution = 0.6,
  algorithm = 1
)
}
\arguments{
\item{object}{input IRIS-FGM object.}

\item{dims}{Dimensions of reduction to use as input.}

\item{k.param}{Defines k for the k-nearest neighbor algorithm.}

\item{resolution}{Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.}

\item{algorithm}{Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.}
}
\value{
It will generate cell type inforamtion.
}
\description{
This function is based on Seurat package.
}
\examples{
\dontrun{object <- RunClassification(object, dims = 1:15, k.param = 20, resolution = 0.6, algorithm = 1)}
}
