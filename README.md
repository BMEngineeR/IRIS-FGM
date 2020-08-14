# IRIS-FGM

IRIS-CEM provides an R environment for integrating [LTMG](https://academic.oup.com/nar/article/47/18/e111/5542876) and [QUBIC2](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz692/5567116)

## Abstract
Single-cell RNA-Seq data is useful in discovering cell heterogeneity and signature genes in specific cell populations in cancer and other complex diseases. Specifically, the investigation of functional gene modules (FGM) can help to understand gene interactive networks and complex biological processes. QUBIC2 is recognized as one of the most efficient and effective tool for FGM identification from scRNA-Seq data. However, its availability is limited to a C implementation and its applicative power is affected by only a few downstream analyses functionalities. We developed an R package named IRIS-FGM (integrative scRNA-Seq interpretation system for functional gene module analysis) to support the investigation of FGMs and cell clustering using scRNA-Seq data. Empowered by QUBIC2, IRIS-FGM can identify co-expressed and co-regulated FGMs, predict types/clusters, identify differentially expressed genes, and perform functional enrichment analysis. It is noteworthy that IRIS-FGM also applies Seurat objects that can be easily used in the Seurat vignettes.

## Workflow
<img src="https://user-images.githubusercontent.com/26455910/90200476-f64ff900-dda5-11ea-902c-c91726c22eac.jpg" alt="IRISFGM_working_flow" width="400" height="700">

## Environment

1. System environment.

We recommend users to install IRIS-FGM on large memory (32GB) based UNIX/Linux operation system if the user aims at analyzing bicluster-based co-expression analysis. If the user seeks to analyze data by quick mode, we recommend installing IRIS-FGM on a small memory (8GB) based Windows or Linux operation system; IRIS-FGM does not support MAC. 

2. R environment:
R (equal or greater than 3.5)

# Usage

## Pre-installation
1. Install required packages from CRAN. 

```install.packages(c('BiocManager','devtools', 'AdaptGauss', "pheatmap", 'mixtools','MCL', 'anocva', 'qgraph','Rtools','ggpubr',"ggraph","Seurat"))```
                   
2. Install required package from Bioconductor.  

```BiocManager::install(c('org.Mm.eg.db','multtest','org.Hs.eg.db','clusterProfiler','DEsingle', 'DrImpute','scater', 'scran'))```

3. Install IRIS-FGM from github.

```devtools::install_github("BMEngineeR/IRIS-FGM", force = T)```

Now the MetaQUBIC is successfully installed and ready for use. 
***
## Data preparation

1. Download example data from [Yan's RPKM 90 cell embroyonic single cell RNA-Seq data](https://bmbl.bmi.osumc.edu/downloadFiles/Yan_expression.txt), and [paper](https://www.nature.com/articles/nsmb.2660).

2. Set the working directory where you put data and import the IRIS-FGM library.
```{r setwd, eval =FALSE, echo = TRUE}
setwd("./my_working_dir/")
library(IRISFGM)
```
## Read in data.
```InputMatrix <- read.table("./my_working_dir/Yan_expression.txt",header = T, row.names = 1)```

IRIS-FGM also provides the function to read in 10X scRNA-Seq data format.

Read HDF5 file```ReadFrom10X_h5()``` or read the folder which contain three files (barcode, sparse matrix, genes)```ReadFrom10X_folder()```

## Analysis data
### Data preprocessing and LTMG modeling
1. **Create IRIS-FGM object**.
```{r create_object, eval= FALSE, echo = TRUE,message=FALSE}
object <- CreateIRISCEMObject(InputMatrix)
```
2. **Add meta information for your cell**: 

This step can add the customized cell label by the user. The format of a file passing to `meta.info` is a data frame of which row name should be cell ID, and column name should be cell type.
```{r add_metadata, eval= FALSE, echo = TRUE}
object <- AddMeta(object, meta.info = NULL)
```

3. **Filter out low quality data**: 

Use `PlotMeta` and `SubsetData` together to filter out low-quality cells. 
```{r plot_metadata,eval= TRUE, echo = TRUE}
PlotMeta(object)
```
<img src="https://user-images.githubusercontent.com/26455910/90202006-6eb8b900-ddaa-11ea-9bcd-38ec1401e88b.png" alt="metadata" width="400" height="400">

```
object <- SubsetData(object , nFeature.upper=15000,nFeature.lower=8000,
                         Counts.upper=700000,Counts.lower=400000)
```
After cells get filtered, we can check the distribution again.
```{r plot_metadata,eval= TRUE, echo = TRUE}
PlotMeta(object)
```
<img src="https://user-images.githubusercontent.com/26455910/90202006-6eb8b900-ddaa-11ea-9bcd-38ec1401e88b.png" alt="metadata" width="400" height="400">

4. Remove low quality data based on the previous plot.
```{r subset_data, echo=FALSE,eval= FALSE, echo =  TRUE}
object <- SubsetData(object , nFeature.upper=15000,nFeature.lower=8000,
                         Counts.upper=700000,Counts.lower=400000)
```

5








