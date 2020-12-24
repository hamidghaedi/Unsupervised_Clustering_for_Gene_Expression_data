# Gene-Expression-Unsupervised-Clusteing
Using consensus clustering approach to find pattern in expression data sets

Unsupervised class discovery is a data mining method to identify unknown possible groups (clusters) of items solely based on intrinsic features and no external variables. Basically clustering includes 4 steps:

#### [1] Data preparation,
#### [2] Dissimilarity matrix calculation,
#### [3] applying clustering algorithms, 
#### [4] Assessing cluster assignment
_________________________________________________________________________________________________________________________________________________________________________________________
In the first step you should calculate a dissimilarity matrix for clustering . Again there is difficulty regarding to math calculation on categorical/binary data. In a dissimilarity matrix distances between all possible pair would be calculated and then subsequent algorithm will use these stats to make clusters. For doing this there are diffrent methods: Gower distance from cluster R base package, methods available in vegan R package a: binomial, chi-square, raup and jaccard . It depends on your data and your decision to chose which method.


For this technique, an investigator seeks to answer two questions: how many groups are present in a dataset, and what is the confidence in the number of groups and the group memberships.

# three  main  steps  to  useConsensusClusterPlus:   
# 1-preparing  inputdata, 
# 2- running the program, and 
# 3-generating cluster-consensus and item-consensus.

