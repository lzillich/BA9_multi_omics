# Differential expression analysis in DESeq2 based on vignette https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# Author: Eric Zillich, last change: 2024-04-04

library(ggplot2)
library(DESeq2)
library(biomaRt)
library(dplyr)
library(data.table)
library(readxl)
library(readr)
library(varhandle)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(qqman)


#### 1. Import phenotype data ####
seq_info <- read.delim("/dir/SampleSheet.txt", header=FALSE)
pheno <- read.csv2("/pheno_dir/pheno.csv")
pheno <- pheno[pheno$Brain_ID %in% seq_info$Brain_ID,]
RIN <- read_excel("/dir/RIN.xlsx")
RIN <- RIN[, c("Brain_ID","RIN")]
RIN$RIN <- as.numeric(RIN$RIN)
RIN <- RIN[RIN$Brain_ID %in% seq_info$Brain_ID,]
pheno2 <- merge(seq_info,pheno,by="Brain_ID")
pheno2 <- merge(pheno2,RIN,by="Brain_ID")

pheno2$Simplified_Axis_1 <- recode(pheno2$Simplified_Axis_1,"MDD" = 1,"Depressive disorder NOS"=1,"Nil"=0)
pheno2$Axis_1_Dependence <- recode(pheno2$Axis_1_Dependence,"Alcohol"=1,"Nil"=0)

# Remove irrelevant pheno data 
pheno2 <- pheno2[,c("Brain_ID","rn","CUD","Age","PMI","pH","RIN","Axis_1_Dependence","Simplified_Axis_1")]

#mean impute missing pH data 
mean(pheno2$pH,na.rm=T)
pheno2$pH[is.na(pheno2$pH)] <- 6.386957

DT <- data.table
DF <- data.frame

#### 2. Import counts ####

# CUD
mRNA_CUD <- get(load("/dir/mRNA_CUD_counts.Rdata"))
colnames(mRNA_CUD$counts) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(mRNA_CUD$counts))
mRNA_CUD$targets <- gsub("Aligned.sortedByCoord.out.bam","",mRNA_CUD$targets)
colnames(mRNA_CUD$stat) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(mRNA_CUD$stat))

mRNA_CUD_counts <- DF(mRNA_CUD$counts)

#### START ANALYSIS ####

#pheno

pheno_mRNA <- pheno2
pheno_mRNA <- pheno_mRNA[order(pheno_mRNA$rn,decreasing = T),]
mRNA_CUD_counts <- mRNA_CUD_counts[,order(colnames(mRNA_CUD_counts),decreasing = T)]

pheno_mRNA$rn  == colnames(mRNA_CUD_counts)
# TRUE, continue 

#coldata
coldata <- as.data.frame(pheno_mRNA)
coldata$CUD <- as.factor(coldata$CUD)
coldata$Age <- scale(coldata$Age)
coldata$PMI <- scale(coldata$PMI)
coldata$RIN <- scale(coldata$RIN)
coldata$pH <- scale(coldata$pH)
coldata$Simplified_Axis_1 <- as.factor(coldata$Simplified_Axis_1)
coldata$Axis_1_Dependence <- as.factor(coldata$Axis_1_Dependence)
rownames(coldata) <- coldata$rn
  
#sort
coldata <- coldata[colnames(mRNA_CUD_counts),]
identical(colnames(mRNA_CUD_counts),rownames(coldata)) # must be TRUE

## Specify model
dds <- DESeqDataSetFromMatrix(countData = mRNA_CUD_counts,
                              colData = coldata,
                              design = ~ Age + RIN + pH + PMI + CUD)
# Prefiltering

keep <- rowSums(counts(dds) >= 2) >= 4

#Preprocessing of data for WGCNA
dds2 <- DESeqDataSetFromMatrix(countData = mRNA_CUD_counts,
                               colData = coldata,
                               design = ~1)
keep2 <- rowSums(counts(dds2) >= 2) >= 4
dds2 <- dds2[keep2,]

vsd_nc <- vst(dds2, blind = F)
vsd_mat_nc <- assay(vsd_nc)
save(vsd_mat_nc, file = "/WGCNA_dir/input_BA9_expression.Rdata")

# DE analysis
dds <- dds[keep,]
dds <- DESeq(dds)

## Define contrasts, extract results table and shrink log2 fold changes
contrast <-  c("CUD", "1", "0")

res <- results(dds, contrast=contrast, alpha = 0.05)

#Order genes by adjusted p-value
res_genes <- as.data.frame(res)
res_genes_ordered <- res_genes[order(res_genes$padj,decreasing = F), ]
write.table(res_genes_ordered,"/results_dir/DE_results.txt",quote=F,row.names = T,sep=",")

# qqplot
png("results_dir/DE_qqplot.png", width = 800, height = 800)
qq(res_genes_ordered$pvalue)
dev.off()

#calculate Lambda
Lambda <- qchisq(median(res_genes_ordered$pvalue,na.rm=T), df = 1, lower.tail = F)/qchisq(0.5, 1)
#Save model number, sample size & lambda to text file:
Model <- "BA9"
tab1 <- cbind(Model, Lambda)
write.table(tab1, file="/results_dir/LAMBDA.txt", row.names=F, quote=F, sep="\t")

# Save a gene list for drug repurposing analysis with 150 top down and up genes based on the test statistics as used for gsea
DE_results <- read.csv("/results_dir/DE_results.txt")
DE_top_up <- DE_results$Gene[order(DE_results$stat,decreasing = T)][1:150]
DE_top_down <- DE_results$Gene[order(DE_results$stat,decreasing = F)][1:150]

write.table(DE_top_up, file="/results_dir/DE_up_top150_drug_repo.txt", row.names=F,col.names = F, quote=F, sep="\t")
write.table(DE_top_down, file="/results_dir/DE_down_top150_drug_repo.txt", row.names=F, col.names = F, quote=F, sep="\t")

##### Model with AUD and MDD additional covariates
## Specify model
dds <- DESeqDataSetFromMatrix(countData = mRNA_CUD_counts,
                              colData = coldata,
                              design = ~ Age + RIN + pH + PMI + Simplified_Axis_1 + Axis_1_Dependence + CUD)

# Prefiltering

keep <- rowSums(counts(dds) >= 2) >= 4
dds <- dds[keep,]

# DE analysis
dds <- DESeq(dds)

## Define contrasts, extract results table and shrink log2 fold changes
contrast <-  c("CUD", "1", "0")

res <- results(dds, contrast=contrast, alpha = 0.05)

#Order genes by adjusted p-value
res_genes <- as.data.frame(res)
res_genes_ordered <- res_genes[order(res_genes$padj,decreasing = F), ]
write.table(res_genes_ordered,"/with_AUD_MDD/DE_results.txt",quote=F,row.names = T,sep=",")


# Sensitivtiy analysis for including MDD and AUD in DE testing

DE_results <- read.csv("/results_dir/DE_results.txt")
DE_results_MDDAUD <- read.csv("/with_AUD_MDD/DE_results.txt")

# Subset for nominal significant association 
DE_results <- DE_results[DE_results$pvalue <0.05,]
DE_results_MDDAUD <- DE_results_MDDAUD[DE_results_MDDAUD$pvalue <0.05,]
merged <- merge(DE_results,DE_results_MDDAUD,by="Gene")

png("/with_AUD_MDD/sensitivity_MDD_AUD.png",height=6,width=6,res=600,units="in")
plot(merged$log2FoldChange.x, merged$log2FoldChange.y, xlab = "log2FC - w/o AUD MDD", ylab = "log2FC - with AUD MDD",xlim=c(-5,5),ylim=c(-5,5))
dev.off()
cor.test(merged$log2FoldChange.x, merged$log2FoldChange.y)# r=0.94, p-value < 2.2e-16
