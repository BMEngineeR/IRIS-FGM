# IRIS-FGM

IRIS-CEM provides an R environment for integrating [LTMG](https://academic.oup.com/nar/article/47/18/e111/5542876) and [QUBIC2](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz692/5567116)

## Abstract
Single-cell RNA-Seq data is useful in discovering cell heterogeneity and signature genes in specific cell populations in cancer and other complex diseases. Specifically, the investigation of functional gene modules (FGM) can help to understand gene interactive networks and complex biological processes. QUBIC2 is recognized as one of the most efficient and effective tools for FGM identification from scRNA-Seq data. However, its availability is limited to a C implementation, and its applicative power is affected by only a few downstream analyses functionalities. We developed an R package named IRIS-FGM (integrative scRNA-Seq interpretation system for functional gene module analysis) to support the investigation of FGMs and cell clustering using scRNA-Seq data. Empowered by QUBIC2, IRIS-FGM can identify co-expressed and co-regulated FGMs, predict types/clusters, identify differentially expressed genes, and perform functional enrichment analysis. It is noteworthy that IRIS-FGM also applies Seurat objects that can be easily used in the Seurat vignettes.

## Workflow
<img src="https://user-images.githubusercontent.com/26455910/90200476-f64ff900-dda5-11ea-902c-c91726c22eac.jpg" alt="IRISFGM_working_flow">

## Environment

1. System environment:

We recommend users to install IRIS-FGM on large memory (32GB) based UNIX/Linux operation system if the user aims at analyzing bicluster-based co-expression analysis. If the user seeks to analyze data by quick mode, we recommend installing IRIS-FGM on a small memory (8GB) based Windows or Linux operation system; IRIS-FGM does not support MAC. 

2. R environment:
R (equal or greater than 3.5)

# Usage

## Pre-installation
1. Install required packages from CRAN: 

```install.packages(c('BiocManager','devtools', 'AdaptGauss', "pheatmap", 'mixtools','MCL', 'anocva', 'qgraph','Rtools','ggpubr',"ggraph","Seurat"))```
                   
2. Install required package from Bioconductor:  

```BiocManager::install(c('org.Mm.eg.db','multtest','org.Hs.eg.db','clusterProfiler','DEsingle', 'DrImpute','scater', 'scran'))```

3. Install IRIS-FGM from github:

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
1. **Create IRIS-FGM object**:
```{r create_object, eval= FALSE, echo = TRUE,message=FALSE}
object <- CreateIRISCEMObject(InputMatrix)
```
2. **Adding meta information for your cell**:

This step can add the customized cell label by the user. The format of a file passing to `meta.info` is a data frame of which row name should be cell ID, and column name should be cell type.
```{r add_metadata, eval= FALSE, echo = TRUE}
object <- AddMeta(object, meta.info = NULL)
```

3. **Filtering out low quality data**:

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
```
PlotMeta(object)
```
<img src="https://user-images.githubusercontent.com/26455910/90202006-6eb8b900-ddaa-11ea-9bcd-38ec1401e88b.png" alt="metadata" width="400" height="400">

4. **Normalization**:

Users can choose to perform normalization based on their needs. The normalization method has two options, one is the simplest CPM normalization (default `normalization = 'LibrarySizeNormalization'`). The other is from package scran and can be opened by using parameter `normalization = 'scran'`. Compared to the CPM normalization method, scran will be more accurate but takes more time.
```
object <- ProcessData(object, normalization = "LibrarySizeNormalization")
```
5. **LTMG modeling**:

Here, we will use Left-truncated Mixture Gaussian distribution to model the regulatory signal of each gene. Parameter, 'Gene_use', decides number of top high variant gene for LTMG modeling, and here we use all genes.
```{r run_LTMG, echo = TRUE,eval = FALSE}
object <- RunLTMG(object, Gene_use = "all", seed = 123)
```
### Biclustering

IRIS-FGM can provide biclustering function, which is based on our in-house novel algorithm, [QUBIC2] (https://github.com/maqin2001/qubic2) to predict functional gene module. 

**Discretization & biclustering**  
If you have more cells to analyze functional gene module you can use LTMG or Quantile based discretization ([QUBIC](https://academic.oup.com/nar/article/37/15/e101/2409951)). In this step, IRIS-FGM will generate some files in the local working directory.

1. LTMG based discretization: 

It might take a lot of time if you have a large size of cells.
```
object <- CalBinaryMultiSignal(object)
object <- RunBicluster(object, DiscretizationModel = "LTMG",OpenDual = FALSE,
                          NumBlockOutput = 100, BlockOverlap = 0.7, BlockCellMin = 15)
```
2. Quantile based discretization: 

It might require users to adjust parameter q for deciding cutoff of binarizing. 
```
object <- RunDiscretization(object, q = 0.06)
object <- RunBicluster(object, DiscretizationModel = "Quantile",OpenDual = TRUE, Extension = 0.90,
                          NumBlockOutput = 1000, BlockOverlap = 0.7, BlockCellMin = 15)
```
**Analyzing functional gene module**
This section is based on the quantile based discretization and biclustering results.

1. Visualize the gene module-to-module relationship globally:

The figure via the number of overlap genes (controlled `edge.by = "gene"`) shows that bicluster 20, 23, 24, and 25 have similar gene components. Bicluster 21, 22, 25, and until bicluster 35 have similar gene components. In this figure, we can see gene modules in bicluster 20, 23, 24, and 25 may have a similar component and similar functionality, whereas the gene components from this gene module may differ from the other gene modules from the other biclusters (i.e., bicluster 21, 22, 25, and until bicluster 35)
```
PlotNetwork(object,N.bicluster =c(20:30), edge.by = "gene")
```

<img src="https://user-images.githubusercontent.com/26455910/90253364-0819b680-de0f-11ea-91a6-a61e1df576e8.png" alt="metadata" width="600" height="300">

2. Visualize bicluster 20 & bicluster 35 on heatmap.

From step 1 in "Analysis functional gene module," we postulate gene components from bicluster 20, 23, 24, and 25 may differ from gene components from bicluster 21, 22, 25, and until bicluster 35. Therefore, in this section, we will focus on how the difference is. Therefore, we use bicluster 20 and bicluster 35 to generate heatmap and show such a difference.
```
PlotHeatmap(object,N.bicluster = c(20,35),show.annotation = F)
```
<img src="https://user-images.githubusercontent.com/26455910/90255235-01d90980-de12-11ea-8469-8992f578ee4b.png" alt="metadata" width="400" height="300">


3. Visualize local co-expression gene module network: 

Since we already know the bicluster 20 and bicluster 35 showing the difference at the global level. We then will focus on a local gene module and investigate the co-expression gene network with in the module. Yellow nodes represent the gene module network from bicluster #20. The nodes' size indicates the degree of presence (the number of connected edges for one node). The thickness of edges suggests the value of the correlation coefficient. 

3.1 From this figure (bicluster 20, click the figure for zooming in), we can tell the EIFAD gene show a negative correlation (red color edge) to GOSR1 & BBS5, and show a positive correlation (grey edge) to ZNF394 & POTEM.
```
PlotModuleNetwork(object, N.bicluster = 20, cutoff=0.6, Node.color = "#E8E504")
```
<img src="https://user-images.githubusercontent.com/26455910/90256029-1a95ef00-de13-11ea-87a4-a4302396df8e.png" alt="metadata" width="600" height="400">

3.2 From this figure (bicluster 35, click the figure for zooming in), we can tell the ROBO1 gene shows a negative correlation (red color edge) to NLRP4 & BPGM. GNPDA2 gene shows a positive correlation to BPGM, KIT, and CCDC25.
```
PlotModuleNetwork(object, N.bicluster = 35, cutoff=0.6, Node.color = "#E8E504")
```
<img src="https://user-images.githubusercontent.com/26455910/90259954-d6a5e880-de18-11ea-8ade-8d822fde3053.png" alt="metadata" width="600" height="400">

4. Functional enrichment analysis.

ISIR-FGM provide a functional enrichment analysis for a selected gene module. 

1. For gene module from bicluster 20: 

```
object <- RunPathway(object ,module.number =20, selected.gene.cutoff = 0.05, species = "Human", database = "GO", genes.source = "Bicluster")
```


|            | ONTOLOGY | ID         | Description                   | GeneRatio | BgRatio  | pvalue               | p.adjust           | qvalue             | geneID          | Count |
|------------|----------|------------|-------------------------------|-----------|----------|----------------------|--------------------|--------------------|-----------------|-------|
| GO:0008278 | CC       | GO:0008278 | cohesin complex               | 2/34      | 16/19717 | 0.000341142521037922 | 0.0313851119354888 | 0.0280095964641662 | STAG3L3/STAG3L2 | 2     |
| GO:0005689 | CC       | GO:0005689 | U12-type spliceosomal complex | 2/34      | 27/19717 | 0.000986050510033922 | 0.0453583234615604 | 0.0404799683066557 | ZMAT5/SNRNP35   | 2     |
| GO:0000062 | MF       | GO:0000062 | fatty-acyl-CoA binding        | 2/32      | 22/17697 | 0.000715377512795969 | 0.0470414266188569 | 0.038513448693801  | ACBD5/SOAT2     | 2     |
| GO:1901567 | MF       | GO:1901567 | fatty acid derivative binding | 2/32      | 28/17697 | 0.00116271921939764  | 0.0470414266188569 | 0.038513448693801  | ACBD5/SOAT2     | 2     |
| GO:0005484 | MF       | GO:0005484 | SNAP receptor activity        | 2/32      | 31/17697 | 0.001425497776329    | 0.0470414266188569 | 0.038513448693801  | STX11/GOSR1     | 2     |


2. For gene module from bicluster 35: 

```
object <- RunPathway(object ,module.number =35, selected.gene.cutoff = 0.05, species = "Human", database = "GO", genes.source = "Bicluster")
```

|            | ONTOLOGY | ID         | Description                                                  | GeneRatio | BgRatio   | pvalue               | p.adjust           | qvalue             | geneID                            | Count |
|------------|----------|------------|--------------------------------------------------------------|-----------|-----------|----------------------|--------------------|--------------------|-----------------------------------|-------|
| GO:0001667 | BP       | GO:0001667 | ameboidal-type cell migration                                | 6/29      | 461/18670 | 6.43045336163777e-05 | 0.0376874533879247 | 0.0301960094109077 | ROBO1/PTPRG/ERBB4/KIT/PTPRR/HDAC9 | 6     |
| GO:0021889 | BP       | GO:0021889 | olfactory bulb interneuron differentiation                   | 2/29      | 12/18670  | 0.000152281291259449 | 0.0376874533879247 | 0.0301960094109077 | ROBO1/ERBB4                       | 2     |
| GO:0010631 | BP       | GO:0010631 | epithelial cell migration                                    | 5/29      | 351/18670 | 0.000186973532900314 | 0.0376874533879247 | 0.0301960094109077 | ROBO1/PTPRG/KIT/PTPRR/HDAC9       | 5     |
| GO:0090132 | BP       | GO:0090132 | epithelium migration                                         | 5/29      | 354/18670 | 0.000194519550646714 | 0.0376874533879247 | 0.0301960094109077 | ROBO1/PTPRG/KIT/PTPRR/HDAC9       | 5     |
| GO:0090130 | BP       | GO:0090130 | tissue migration                                             | 5/29      | 360/18670 | 0.000210309449709401 | 0.0376874533879247 | 0.0301960094109077 | ROBO1/PTPRG/KIT/PTPRR/HDAC9       | 5     |
| GO:0042692 | BP       | GO:0042692 | muscle cell differentiation                                  | 5/29      | 385/18670 | 0.000286907770833753 | 0.0428448937778404 | 0.0343282631067754 | KIT/SYNE1/TMOD2/UCHL1/HDAC9       | 5     |
| GO:0005001 | MF       | GO:0005001 | transmembrane receptor protein tyrosine phosphatase activity | 2/32      | 17/17697  | 0.000423558866693072 | 0.0222368405013863 | 0.0164965032290986 | PTPRG/PTPRR                       | 2     |
| GO:0019198 | MF       | GO:0019198 | transmembrane receptor protein phosphatase activity          | 2/32      | 17/17697  | 0.000423558866693072 | 0.0222368405013863 | 0.0164965032290986 | PTPRG/PTPRR                       | 2     |




