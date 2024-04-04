# Author: Eric Zillich, last change: 2024-04-04
## Data Preparation MOFA - adapted from Tutorial code https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/CLL.html and https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html

library(data.table)
library(dplyr)
library(readr)
library(varhandle)
library(matrixStats)
library(data.table)
library(ggplot2)
library(tidyverse)
library(psych)
library(ggpubr)
library(limma)
library(biomaRt)

DF <- data.frame

setwd(dir0)

# Import Illumina manifest
manifest <- read.csv("/Annotation_dir/EPIC_anno.txt")

promoter_cpgs <- manifest[manifest$UCSC_RefGene_Group %in% c("TSS200","TSS1500"),]

#read and create phenofile with combined meth and expr data
pheno_Meth <- read.csv("/EWAS_dir/pheno_EWAS.txt", sep=";")

pheno_Expr <- read.csv("pheno_dir/pheno.txt", sep=";")
pheno_Expr <-pheno_Expr[,c(1,2,7)]

all_merged <- merge(pheno_Expr, pheno_Meth, by = "Brain_ID", all.y = T)
all_merged <- all_merged[complete.cases(all_merged),]

write.table(all_merged , file="/MOFA_dir/input/pheno_meth_expr.txt", row.names=F, quote=F, sep="\t")

##### Methylation #####

mval <- get(load("/EWAS_dir/output/m_values_filtered.Rdata"))
meth <- get(load("/EWAS_dir/output/beta_QN.RData"))

# Extract beta values for relevant CpGs - remove filtered sites from the beta value matrix and remove subjects that do not have expression data
meth <- meth[rownames(meth) %in% rownames(mval),colnames(meth) %in% all_merged$rn.y]

# Extract promoter CpGs only 
meth2 <- meth[rownames(meth) %in% promoter_cpgs$Name,]

# Two-step approach: 1. select all promoter CpGs for the methylation data input, 2. select only the top 20000 promoter CpGs showing the highest variance to balance the sizes between expression and methylation datasets
  var_met <- rowVars(meth2)
  meth_var <- DF(rownames(meth2),var_met)
  
  #filter by variance
  names(meth_var)[names(meth_var) == "rownames.meth2."] <- "Name"
  meth_var_pro <- merge(meth_var,promoter_cpgs,by ="Name")
  meth_var_pro <- meth_var_pro[order(meth_var_pro$var_met, decreasing = T),]
  meth_var_pro <- meth_var_pro[1:20000,]
  
  meth3 <- meth2[rownames(meth2) %in% meth_var_pro$Name,]
  
  #get unique gene names for pasting them to the cg name
  gene_list <- str_split(promoter_cpgs$UCSC_RefGene_Name,";")

  for(i in 1:length(gene_list)){
    l <- gene_list[[i]]
    gene_list[[i]] <- unique(l)
  }
  promoter_cpgs$genes2 <- ""
  
  for(i in 1:length(promoter_cpgs$genes2)){
    promoter_cpgs$genes2[i] <- paste(unlist(gene_list[[i]]),collapse=";")
  }
  
  colnames(meth3) <- all_merged$Brain_ID[match(colnames(meth3),all_merged$rn.y)]
  rownames(meth3) <- paste0(rownames(meth3)," - ",promoter_cpgs$genes2[match(rownames(meth3),promoter_cpgs$Name)])
  
  
  
  ####### EXPRESSION #######
  norm_counts = get(load("/WGCNA_dir/input_BA9_expression.Rdata"))
  
  #filter out low variable genes - keep top 20000 variable genes
  ntop<-20000
  sds <- genefilter::rowSds(norm_counts)
  exprMat<-norm_counts[order(sds,decreasing = T)[1:ntop],]
  
  #scale and center RNAseq data
  exprMat <- jyluMisc::mscale(exprMat, center = TRUE, scale = TRUE, useMad = FALSE)
  colnames(exprMat) <- all_merged$Brain_ID[match(colnames(exprMat),all_merged$rn.x)]
  
  list_mofa <- list( meth3, exprMat)
  
  save(list_mofa, file ="/MOFA_dir/input/BA9_meth_expr_Rdata")
  

