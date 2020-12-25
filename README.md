# Gene-Expression-Unsupervised-Clusteing
Using consensus clustering approach to find pattern in expression data sets

Unsupervised class discovery is a data mining method to identify unknown possible groups (clusters) of items solely based on intrinsic features and no external variables. Basically clustering includes 4 steps:

#### [1] Data preparation,
#### [2] Dissimilarity matrix calculation,
#### [3] applying clustering algorithms, 
#### [4] Assessing cluster assignment
_________________________________________________________________________________________________________________________________________________________________________________________

### [1] Data preparation

In this step we need to convert/transform/normalize gene expression data. Usually analysis would start from a raw count matrix comming from an RNA-seq experiment. Because there are a large number of features (gene) in such matrix , a simple feature selection step should be done to limit the analysis to those genes that possibly explain variation between samples in the cohort.  
