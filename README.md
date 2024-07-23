# SRPairs
Thanks for supporting the SRPairs package written by Minyu Nong. WeChat official account: shengxinjy123
```
install.packages("devtools")
devtools::install_github("MYNong123/SRPairs")
```
Use this code to install the SRPairs package on your R.
# Log
The current version is V2.1.2
# Description
SRPairs package is more convenient for bioinformatics researchers to find stable reversal gene pairs in gene expression profile.
If GeneA is larger than GeneB in one group of samples, and GeneA is smaller than GeneB in another group of samples, then the gene symmetry of these two genes is a relative expression orderings(REOs).
Stable expression orders (REOs), that is, a pattern of gene order. Assuming that gene i (Gi) and gene j (Gj) meet Ga > Gb (or Ga < Gb) in 95% of normal samples, we call them stable gene pairs.
However, if it is the opposite in tumor samples, such as Ga > Gb in normal samples and Ga < Gb in tumor samples, then such genes are called Stable expression orders (REOs).
REOs can better solve the problem of batch effect in differential expression analysis, have stronger robustness, and can qualitatively study the relationship between gene expressions. Therefore, it can provide preliminary screening of characteristic genes for researchers.
# Example
Prepare the input data of Normal.txt and Tumor.txt (one gene for each behavior, and one sample for each column).
# Sample code
```
Normal = read.table("Normal.txt",header=T,sep="\t",check.names=F,row.names=1) #read Normal matrix.
Tumor = read.table("Tumor.txt",header=T,sep="\t",check.names=F,row.names=1) #read Tumor matrix.
#library packages.
library(SRPairs)
library(stringr)
PairTime=testRunTime(Normal,Tumor,0.9,1) #Test the running speed of the computer.#Here, 0.9 represent the proportion of HighPercent expressed in the sample, indicating that 90% of the samples of this gene are larger than the normal group. The default value is 0.9, which can also be set to 0.95.
GetPair=FindStabRvsPair(Normal,Tumor,0.9,1,showbar=T) #get REOs
GenePair=CalStabRvsPair(GetPair,Normal,Tumor) #get REOs matrix.
```
# Sample Data
```
data(Sample_Data)
```
Get Sample Data.
# Downstream analysis
It can be used as a machine learning algorithm for feature screening such as SVM, random forest and LASSO to provide input data, and can also be used to build an artificial neural network.
# Multicore operation
Bugs in multi-core computing are still being fixed.
