# Applied CIBERSORT based on documentation files from https://cibersortx.stanford.edu/
# Author: Eric Zillich, last change: 2024-04-04

source("/CIBERSORT_dir/CIBERSORT.R")

# DESeq2 data import 
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
library(BayesianFirstAid)


#### 1. Import phenotype data ####
pheno <- read.csv("/pheno_dir/pheno.txt", sep=";")

#### 2. Import counts ####

# CUD
mRNA_CUD <- get(load("/dir/mRNA_CUD_counts.Rdata"))
colnames(mRNA_CUD$counts) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(mRNA_CUD$counts))
mRNA_CUD$targets <- gsub("Aligned.sortedByCoord.out.bam","",mRNA_CUD$targets)
colnames(mRNA_CUD$stat) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(mRNA_CUD$stat))

mRNA_CUD_counts <- DF(mRNA_CUD$counts)

#pheno
rownames(pheno) <- pheno$rn

#coldata
coldata <- pheno[,-c(1,2)]
coldata$CUD <- as.factor(coldata$CUD)
coldata$Age <- scale(coldata$Age)
coldata$PMI <- scale(coldata$PMI)
coldata$RIN <- scale(coldata$RIN)
coldata$pH <- scale(coldata$pH)
coldata$Axis_1_Dependence <- as.factor(coldata$Axis_1_Dependence)
coldata$Simplified_Axis_1 <- as.factor(coldata$Simplified_Axis_1)

#sort
coldata <- coldata[colnames(mRNA_CUD_counts),]
identical(colnames(mRNA_CUD_counts),rownames(coldata)) # must be TRUE

## Specify model
dds <- DESeqDataSetFromMatrix(countData = mRNA_CUD_counts,
                              colData = coldata,
                              design = ~ Age + RIN + pH + PMI + Axis_1_Dependence + Simplified_Axis_1 + CUD)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds,normalized=T)

#### Yu & He reference datasets ####

# Import 10fold dataset
markers10 <-read_excel("reference/Yu_He_markers_10fold.xls")
markers10 <- as.data.frame(markers10[,-c(1)])
rownames(markers10) <- markers10$geneSymbol
markers10 <- markers10[,-c(1,2,5,6)]

# Apply CIBERSORT
results10 <- CIBERSORT(markers10,norm_counts, QN=F, absolute=F)
write.table(results10,"/CIBERSORT_dir/Markers_10fold_Yu_He_BA9.txt",sep=";",quote=F)

results10 <- read.csv("/CIBERSORT_dir/Markers_10fold_Yu_He_BA9.txt", sep=";")
res10_2 <- results10[,c(1:7)]
colnames(res10_2) <- gsub("fpkm_","",colnames(res10_2))
colnames(res10_2)[1]<-"rn"

# merge with pheno
pheno <- read.csv("/pheno_dir/pheno.txt", sep=";")
merged <- merge(pheno,res10_2,by="rn")
merged <- merged[order(merged$CUD,decreasing=T),]
merged$name <- c(paste0("CUD-",c(1:13)),paste0("Ctrl-",c(1:12)))
rownames(merged) <- merged$name
merged <- merged[,c(8:13)]
write.table(merged,"/CIBERSORT_dir/celltypes.txt",sep=";",quote=F)

perc_df <- data.frame(samples=rep(rownames(merged),each=6),celltypes=rep(colnames(merged),times=25),value=as.vector(t(merged)))
perc_df$samples <- factor(perc_df$samples,levels=rownames(merged))

# Plot the fraction of celltypes 
per <- ggplot(perc_df, aes(fill=celltypes, x=value, y=samples)) + 
  geom_bar(position="fill", stat="identity")+scale_y_discrete(limits=rev(rownames(merged)))+scale_fill_manual(values=c("#F7AB64","#74BA59","#70305A","#E8326D", "#3A9BCC","#85CEE4"))+xlab("celltype proportion")+ylab(NULL)+
  theme_minimal()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(plot = per,filename = "/CIBERSORT_dir/celltype_percentage.pdf", width = 6, height=6)


# Bayesian testing for the cell type distributions
markers <- read.csv("/CIBERSORT_dir/Markers_10fold_Yu_He_BA9.txt", sep=";")
pheno <- read.csv("pheno_dir/pheno.txt", sep=";")

colnames(pheno)[2] <- "sample"
markers2 <- merge(markers, pheno,by="sample")

# Astrocytes 
astro_CUD <- markers2$astrocytes[markers2$CUD == 1]
astro_Ctrl <- markers2$astrocytes[markers2$CUD == 0]
fit <- bayes.t.test(astro_CUD, astro_Ctrl, paired = F,alternative = "two.sided",mu=mean(astro_CUD)-mean(astro_Ctrl))
plot(fit) # 95% HDI -0.245  0.040

# endothelial 

endo_CUD <- markers2$endothelial[markers2$CUD == 1]
endo_Ctrl <- markers2$endothelial[markers2$CUD == 0]
fit <- bayes.t.test(endo_CUD, endo_Ctrl, paired = F,alternative = "two.sided",mu=mean(endo_CUD)-mean(endo_Ctrl))
plot(fit) # 95% HDI -0.0006  0.0024

# microglia
micro_CUD <- markers2$microglia[markers2$CUD == 1]
micro_Ctrl <- markers2$microglia[markers2$CUD == 0]
fit <- bayes.t.test(micro_CUD, micro_Ctrl, paired = F,alternative = "two.sided",mu=mean(micro_CUD)-mean(micro_Ctrl))
plot(fit) # 95% HDI -0.0065  0.010

# neurons 
neurons_CUD <- markers2$neurons[markers2$CUD == 1]
neurons_Ctrl <- markers2$neurons[markers2$CUD == 0]
fit <- bayes.t.test(neurons_CUD, neurons_Ctrl, paired = F,alternative = "two.sided",mu=mean(neurons_CUD)-mean(neurons_Ctrl))
plot(fit) # 95% HDI -0.039  0.28

# oligodendrocytes 
oligo_CUD <- markers2$oligodendrocytes[markers2$CUD == 1]
oligo_Ctrl <- markers2$oligodendrocytes[markers2$CUD == 0]
fit <- bayes.t.test(oligo_CUD, oligo_Ctrl, paired = F,alternative = "two.sided",mu=mean(oligo_CUD)-mean(oligo_Ctrl))
plot(fit) # 95% HDI -0.013  0.055

# OPC
OPC_CUD <- markers2$OPC[markers2$CUD == 1]
OPC_Ctrl <- markers2$OPC[markers2$CUD == 0]
fit <- bayes.t.test(OPC_CUD, OPC_Ctrl, paired = F,alternative = "two.sided",mu=mean(OPC_CUD)-mean(OPC_Ctrl))
plot(fit) # 95% HDI -0.0031  0.0052
