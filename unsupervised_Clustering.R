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
vsd <- assay(vst(dds)) # For a fully unsupervised transformation one can set blind = TRUE (which is the default).
#_________________________________# Feature Selection _________________________________#
# top 5K based on MAD 

mads=apply(vsd,1,mad)
# check data distribution
hist(mads, breaks=nrow(vsd)*0.1)
# selecting features
mad2k=vsd[rev(order(mads))[1:2000],]
#mad4k=vsd[rev(order(mads))[1:4000],]
#mad6k=vsd[rev(order(mads))[1:6000],]

#_________________________________# Clustering & Cluster assignmnet validation _________________________________#
# finding optimal clusters by CC
library(ConsensusClusterPlus)
results = ConsensusClusterPlus(mad2k,
                               maxK=6,
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title= "geneExp",
                               clusterAlg="pam",
                               distance="spearman",
                               seed=1262118388.71279,
                               plot="pdf")

#_________________________________#Assessing cluster assignment _________________________________#
#Silhouette width analysis
cc4 = results[[4]]

# calcultaing Silhouette width using the cluster package 
library(cluster)
cc4Sil = silhouette(x = cc4[[3]], # x is a numeric vector that indicates cluster assignment for each data point
                    dist = as.matrix(1- cc4[[4]])) # dist should be a distance matrix and NOT a similarity matrix, so I subtract the matrix from one to get that dist matrix

#For visualization:
library(factoextra)
fviz_silhouette(cc4Sil, palette = "jco",
                ggtheme = theme_classic())

#_________________________________#Assessing cluster assignment _________________________________#
#Survival analysis
### preparing dataset for survival analysis
cc4Class = data.frame(cc4$consensusClass)
cc4Class$ID = rownames(cc4Class)
cc4Class = data.frame(cc4Class[match(rownames(clinical.exp),cc4Class$ID),])
all(cc4Class$ID == rownames(clinical.exp))

# new encoding for status, time and cluster
clinical.exp$status = ifelse(clinical.exp$Progression.to.T2. == "NO", 0,1)
clinical.exp$time = as.numeric(clinical.exp$Progression.free.survival..months. * 30)
clinical.exp$cluster = as.factor(cc4Class$cc4.consensusClass)

library(survival)

res.cox <- coxph(Surv(time, status) ~ cluster, data = clinical.exp)
res.cox

summary(res.cox)








