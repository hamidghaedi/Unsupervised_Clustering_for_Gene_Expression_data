# Gene-Expression-Unsupervised-Clusteing
Using consensus clustering approach to find pattern in expression data sets

Unsupervised class discovery is a data mining method to identify unknown possible groups (clusters) of items solely based on intrinsic features and no external variables. Basically clustering includes four steps:

#### [1] Data preparation,
#### [2] Dissimilarity matrix calculation,
#### [3] applying clustering algorithms, 
#### [4] Assessing cluster assignment
_________________________________________________________________________________________________________________________________________________________________________________________

### [1] Data preparation

In this step we need to filter out incomplete cases and low expressed genes, then transforming/normalizing gene expression values. Usually analysis would start from a raw count matrix comming from an RNA-seq experiment. Because there are a large number of features (gene) in such matrix , a feature selection step should be done to limit the analysis to those genes that possibly explain variation between samples in the cohort.  
To do so ;
- For filteration: I keep those genes  that  have expression in 10% of samples with a count of 10 or higher. 
- For transformation/normalization : I  use variance stabilizing transformation (VST). I use ```vst``` function from ```DESeq2 packages``` which at the same time will normalize the raw count also. Using other type of transformation like Z score, log transformation are also quiet commmon.
- For feature selection: I  select 2k, 4k and 6k top genes based on median absolute deviation (MAD) . A number of other methods like "feature selection based on the most variance", "feature dimension reduction and extraction based on Principal Component Analysis (PCA)", and "feature selection based on Cox regression model" could be applied. See bioconductor package [CancerSubtypes manual](http://www.bioconductor.org/packages/release/bioc/html/CancerSubtypes.html). 

Other approaches also could be used for example: 
 - Using log2(RSEM) gene expression value and removing genes with NA values more than 10% across samples. Then selected top 25% most-varying genes by standard deviation of gene expression across samples ([A. Gordon Robertson et al., Cell,2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5687509/)). 

 - Using FPKM matrix and keeping genes with (log2(FPKM+1)>2 in at least 10% of samples and selecting a subsets (2K, 4K, 6K ) of MAD ranked genes to identify stable classes ([Jakob Hedegaard et al, Cancer Cell, 2016](https://www.sciencedirect.com/science/article/pii/S1535610816302094#mmc1)). 
 

```R
#_________________________________reading data and processing_________________________________#
# reading count data
rna <- read.table("Uromol1_CountData.v1.csv", header = T, sep = ",", row.names = 1)
head(rna[1:5, 1:5], 5)
#                   U0001 U0002 U0006 U0007 U0010
#ENSG00000000003.13  1458   228  1800  3945   293
#ENSG00000000005.5      0     0     9     0     0
#ENSG00000000419.11   594    23   792  1378   139
#ENSG00000000457.12   548    22  1029   976   148
#ENSG00000000460.15    53     2   190   136    47

dim(rna)
# [1] 60483   476  this is a typical output from hts-seq count matrix with more than 60,000 genes

# dissecting dataset based on gene-type
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(attributes= c("ensembl_gene_id","hgnc_symbol", "gene_biotype"), 
               mart= mart)

# see gene types returned by biomaRt
data.frame(table(genes$gene_biotype))[order(- data.frame(table(genes$gene_biotype))[,2]), ]

# we will continue with protein coding here:
rna <- rna[substr(rownames(rna),1,15) %in% genes$ensembl_gene_id[genes$gene_biotype == "protein_coding"],]

dim(rna)
#[1] 19581   476 These are protein coding genes

# reading associated clinical data
clinical.exp <- read.table("uromol_clinic.csv", sep = ",", header = T, row.names = 1)
head(clinical.exp[1:5,1:5], 5)
#      UniqueID   CLASS  BASE47    CIS X12.gene.signature
#U0603    U0603 luminal luminal no-CIS          high_risk
#U0497    U0497 luminal   basal no-CIS           low_risk
#U0839    U0839 luminal luminal    CIS           low_risk
#U1043    U1043 luminal luminal no-CIS          high_risk
#U0566    U0566 luminal   basal no-CIS           low_risk

# making sure about sample order in rna and clincal.exp dataset
all(rownames(clinical.exp) %in% colnames(rna))
#[1] TRUE
all(rownames(clinical.exp) == colnames(rna))
#[1] FALSE
# reordering rna dataset
rna <- rna[, rownames(clinical.exp)]
#______________________ Data tranformation & Normalization ______________________________#
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = rna,
                              colData = clinical.exp,
                              design = ~ 1) # 1 passed to the function because of no model
# pre-filteration, however while using DESeq2 package it is not necessary, because the function automatically will filter out low count genes
# Keeping genes with expression in 10% of samples with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= round(ncol(rna)*0.1)
dds <- dds[keep,]

# vst tranformation
vsd <- vst(dds) # For a fully unsupervised transformation one can set blind = TRUE (which is the default).

# Inspecting deseq object
#head(assay(vsd), 3)
#colData(vsd)
#_________________________________# Feature Selection _________________________________#
# top 5K based on MAD 

mads=apply(assay(vsd),1,mad)
# check data distribution
hist(mads, breaks=nrow(assay(vsd))*0.1)
# selecting features
mad2k=assay(vsd)[rev(order(mads))[1:2000],]
#mad4k=assay(vsd)[rev(order(mads))[1:4000],]
#mad6k=assay(vsd)[rev(order(mads))[1:6000],]

```

_________________________________________________________________________________________________________________________________________________________________________________________

### [2&3] Dissimilarity matrix calculation & applying clustering algorithms:

Clustering is grouping similar samples into one cluster and keeping them far from disimilar samples based on distance measure. There are a diverse list of dissimilarity matrix calculation methods (distance measures). To see list of available distance measures please check ```?stats::dist``` and ```?vegan::vegdist()```. The latter need to have the ```vegan``` package to be installed:

   * For log-transformed gene expression, Euclidean based measures can be applied.
   * For RNA-seq normalised counts, correlation based measures (Pearson, Spearman) or a Poisson-based distance can be used.  

**Clustering algorithms**: To see a list of algorithms please check ```diceR``` package [vignettes](https://cran.r-project.org/web/packages/diceR/vignettes/overview.html). Most widely used algorithms are *partitional clustering* and *hierarchical clustering*.

*Partitional clustering* are clustering methods used to classify samples into multiple clusters based on their similarity. The algorithms required to specify the number of clusters to be generated. The commonly used partitional clustering, including:

 * K-means clustering (KM): each cluster is represented by the center or means of the data points belonging to the cluster. This method is sensitive to outliers.
 * K-medoids clustering or PAM (Partitioning Around Medoids), in which, each cluster is represented by one of the objects in the cluster. PAM is less sensitive to outliers.

*Hierarchical clustering (HC)* in contrast to partitional clustering, this method does not require to pre-specify the number of clusters to be generated. HC can be grouped into two classes.
* Agglomerative:"bottom-up" approach, each observation is initially considered as a cluster of its own (leaf), and pairs of clusters are merged as one moves up the hierarchy.
* Divisive: "top-down" approach,  This begins with the root so all observations start in one cluster, and splits are performed recursively as one moves down the hierarchy.

Some words on how to compute distances between clusters: Indeed there are diffrent way to do cluster agglomeration (i.e, linkage ). When one using ```stat``` package, possible methods include “ward.D”, “ward.D2”, “single”, “complete”, “average”, “mcquitty”, “median” or “centroid”. Generally complete linkage and Ward’s method are preferred. More reading materials on this can be find [here](https://www.datanovia.com/en/lessons/agglomerative-hierarchical-clustering/). As a result hierarchical clustering provides a tree-based representation of the objects, which is also known as dendrogram.

In practice HC, KM and PAM are commonly used for gene exoression data.

To quantitatively dtermine the number and membership of possible clusters within the dataset, I will use Consensus Clustering (CC) approach. Applying this method has proved to be effective in new cancer subclasses discoveries. For more information on the methodology please refere to the seminal paper by [Monti et al. (2003)](https://link.springer.com/article/10.1023/A:1023949509487)  and ```ConsensusClusterPlus``` package [manual](https://bioconductor.org/packages/release/bioc/html/ConsensusClusterPlus.html). 

To determine the number of cluster based on CC, there are several graphics which would help to this extent. (1) A Color-coded heatmap corresponding to the consensus matrix that represent consensus values from 0–1; white corresponds to 0 (never clustered togather) and dark blue to 1 (allways clustered togather). (2) Consensus Cumulative Distribution Function (CDF) plot, This graphic lets one to determine at what number of clusters,k, the CDF reaches an approximate maximum. So consensus and cluster confidence is at a maximum at this k. (3) Based on CDF plot, a chart for relative change in area under CDF curve can be plotted. This  plot  allows  a  user  to determine the relative increase in consensus and determine k at which there is no appreciable increase.

```R
#_________________________________# Clustering & Cluster assignmnet validation _________________________________#
# finding optimal clusters by CC
library(CancerSubtypes)
cc.res=ExecuteCC(clusterNum=4, ## Refer below for more details on this
                 d=mad2k,
                 maxK=10,# maximum cluster number for Consensus Clustering Algorithm to evaluate
                 clusterAlg="pam",
                 distance="spearman",
                 title="UROMOL",
                 plot= pdf)
```
By ```clusterNum``` argument I have to provide the number of cluster that I am intrested to get. This is needed by the package ```CancerSubtypes```. The main package to perform CC in R is ```ConsensusClusterPlus``` and it does not need to specify cluster number. From where I got the number 4? From inspecting the  CDF plots in the result folder (here "my/wd/UROMOL").Also I performed same analysis for datasets ```mad4k``` and ```mad6k```. So the number for clusters seems to be 4. This is not in agreement with the original [paper]( https://www.sciencedirect.com/science/article/pii/S1535610816302094), where the authors reported three subtypes. However in new report from the same group using a new pipeline they conculded that NMIBC to have four diffrent subtypes! [ref](https://www.nature.com/articles/s41467-021-22465-w). 

The above command with return two plots which is helful to make decision about cluster number:  consensus CDF and  relative change in area under CDF curve.
                                                                                              
![alt-text-1](https://github.com/hamidghaedi/Gene-Expression-Unsupervised-Clusteing/blob/main/consensus011.png "title-1") ![alt-text-2](https://github.com/hamidghaedi/Gene-Expression-Unsupervised-Clusteing/blob/main/consensus012.png "title-2")

### [4] Assessing cluster assignment
Assessing cluster assignment or cluster validation indicate to the  procedure of assessing the goodness of clustering  results. [Alboukadel Kassambara](https://www.datanovia.com/en/lessons/cluster-validation-statistics-must-know-methods/) has published a detailed pot on this topic. In thi tutorial I will use Silhouette method for cluster assessment.  this method can be used to investigate the separation distance between the obtained clusters. The silhouette plot reflects a measure of how close each data point in one cluster is to a points in the neighboring clusters. This measure, Silhouette width, has a range of -1 to +1. Value near +1 show that the sample is far away from the closeset data point from neighboring cluster. A negative value may indicate wrong cluster assignment and a value close to 0 means an arbitrary cluster assignment to that data point.
_________________________________________________________________________________________________________________________________________________________________________________________
### Refrences
1- Biostar posts:
https://www.biostars.org/p/321773/
https://www.biostars.org/p/225315/
https://www.biostars.org/p/74223/
https://www.biostars.org/p/281161/
https://www.biostars.org/p/273107/

2- https://www.datanovia.com/en/courses/partitional-clustering-in-r-the-essentials/

3- https://2-bitbio.com/2017/10/clustering-rnaseq-data-using-k-means.html


