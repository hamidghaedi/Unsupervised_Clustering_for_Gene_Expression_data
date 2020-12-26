# Gene-Expression-Unsupervised-Clusteing
Using consensus clustering approach to find pattern in expression data sets

Unsupervised class discovery is a data mining method to identify unknown possible groups (clusters) of items solely based on intrinsic features and no external variables. Basically clustering includes 4 steps:

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

Other approaches also could be use for example: 
 - Using log2(RSEM) gene expression value and removing genes with NA values more than 10% across samples. Then selected top 25% most-varying genes by standard deviation of gene expression across samples ([A. Gordon Robertson et al., Cell,2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5687509/)). 

 - Using FPKM matrix and keeping genes with (log2(FPKM+1)>2 in at least 10% of samples and selecting a subsets (2K, 4K, 6K ) of MAD ranked genes to identify stable classes ([Jakob Hedegaard et al, Cancer Cell, 2016](https://www.sciencedirect.com/science/article/pii/S1535610816302094#mmc1)). 
 

```R
# reading count data
rna <- read.table("Uromol1_CountData.v1.csv", header = T, sep = ",")
head(rna[1:5, 1:5], 5)
#                   U0001 U0002 U0006 U0007 U0010
#ENSG00000000003.13  1458   228  1800  3945   293
#ENSG00000000005.5      0     0     9     0     0
#ENSG00000000419.11   594    23   792  1378   139
#ENSG00000000457.12   548    22  1029   976   148
#ENSG00000000460.15    53     2   190   136    47
