# Author: Eric Zillich, last change: 2024-04-04

source("/CIBERSORT_dir/CIBERSORT.R")

library(data.table)
library(variancePartition)
library(readxl)
library(BiocParallel)

#read pheno
pheno <- read.csv("/dir/pheno.txt", sep=";")
pheno$pH[is.na(pheno$pH)] <- 6.335 # mean impute missing pH

colnames(pheno)[8] <-"AUD"
colnames(pheno)[9] <-"MDD"

pheno$CUD <- as.factor(pheno$CUD)
pheno$AUD <- as.factor(pheno$AUD)
pheno$MDD<- as.factor(pheno$MDD)

rownames(pheno) <- pheno$rn

pheno <- pheno[,c("CUD","PMI","pH","RIN","Age","AUD","MDD")]

#specify model
mRNA_CUD <- get(load("dir/mRNA_CUD_counts.Rdata"))
colnames(mRNA_CUD$counts) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(mRNA_CUD$counts))
mRNA_CUD <- as.data.frame(mRNA_CUD$counts)
keep <- rowSums(mRNA_CUD) >= 1
mRNA_CUD <- mRNA_CUD[keep,]

form <- ~ Age + pH + PMI + RIN + (1|AUD) + (1|MDD) 

varPart <- fitExtractVarPartModel(mRNA_CUD, form, pheno,BPPARAM = MulticoreParam(20))
vp <- sortCols(varPart)

png("dir/Variance_partition.png",width = 10,height=6,units = "in",res=300)
plotVarPart(vp)
dev.off()

