# Author: Eric Zillich, last change: 2024-04-04
# Genes and Pathways associated with CUD in BA9 - evidence from multiple omics datasets

library(stringr)
library(UpSetR)
library(ComplexUpset)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)

# Import the different results files for methylation, expression, WGCNA and splicing

# EWAS
EWAS <- read.csv("/results_ranking/CUD_EWAS_results_lm.txt", sep="")

# WGCNA module methylation 
wgcna_meth <- read.csv("/results_ranking/Meth_geneInfo0_promoter.csv", row.names=1)

# RNA-Seq
DE <- read.csv("/results_ranking/DE_results.txt")

# WGCNA module expression 
wgcna_expr <- read.csv("/results_ranking/geneInfo0_Expr.txt", sep=";")

# Splicing
load("/leafviz_dir/leafviz.RData")
clusters <- merge(introns,clusters, by="clusterID")
clusters <- clusters[clusters$FDR <0.05 & abs(clusters$deltapsi) > 0.025,]
clusters$gene.x <- gsub("<i>","",clusters$gene.x)
clusters$gene.x <- gsub("</i>","",clusters$gene.x)
AS_genes <- unique(clusters$gene.x[clusters$gene.x != "." & clusters$FDR < 0.05])

# MOFA
MOFA <- read.table("/results_ranking/Factor9_mRNA_loadings.txt", quote="\"", comment.char="")

MOFA2 <- read.table("/results_ranking/Factor9_meth_loadings.txt", quote="\"", comment.char="")

# Extract genes from the individual analyses 

EWAS <- EWAS$Gene[EWAS$P_VAL < 0.001]
EWAS <- unique(unlist(str_split(EWAS,pattern=";")))[-1]

meth_module_blue <- wgcna_meth$geneSymbol[wgcna_meth$moduleColor == "blue"]
meth_module_blue <- unique(unlist(str_split(meth_module_blue,pattern=";")))
meth_module_steelblue <- wgcna_meth$geneSymbol[wgcna_meth$moduleColor == "steelblue"]
meth_module_steelblue<- unique(unlist(str_split(meth_module_steelblue,pattern=";")))
meth_module_brown <- wgcna_meth$geneSymbol[wgcna_meth$moduleColor == "brown"]
meth_module_brown<- unique(unlist(str_split(meth_module_brown,pattern=";")))
meth_module_brown4 <- wgcna_meth$geneSymbol[wgcna_meth$moduleColor == "brown4"]
meth_module_brown4<- unique(unlist(str_split(meth_module_brown4,pattern=";")))

DE_genes <- DE$Gene[DE$pvalue < 0.05]

expr_module_yellow <- unique(wgcna_expr$geneSymbol[wgcna_expr$moduleColor == "yellow"])

MOFA_Expr <- MOFA$V1
MOFA_Meth <- MOFA2$V1

# Create a list of the genes 

BA9_list <- list(EWAS=EWAS,meth_module_blue=meth_module_blue,meth_module_steelblue=meth_module_steelblue,meth_module_brown=meth_module_brown,meth_module_brown4=meth_module_brown4,DE_genes=DE_genes,expr_module_yellow=expr_module_yellow, AS_genes=AS_genes, MOFA_Expr = MOFA_Expr,MOFA_Meth=MOFA_Meth)

#create data frame from the list for export into the supplement
df <- as.data.frame(matrix(nrow=9201,ncol=10,""))
colnames(df) <- names(BA9_list)
df[,1][1:length(BA9_list[[1]])] <- BA9_list[[1]]
df[,2][1:length(BA9_list[[2]])] <- BA9_list[[2]]
df[,3][1:length(BA9_list[[3]])] <- BA9_list[[3]]
df[,4][1:length(BA9_list[[4]])] <- BA9_list[[4]]
df[,5][1:length(BA9_list[[5]])] <- BA9_list[[5]]
df[,6][1:length(BA9_list[[6]])] <- BA9_list[[6]]
df[,7][1:length(BA9_list[[7]])] <- BA9_list[[7]]
df[,8][1:length(BA9_list[[8]])] <- BA9_list[[8]]
df[,9][1:length(BA9_list[[9]])] <- BA9_list[[9]]
df[,10][1:length(BA9_list[[10]])] <- BA9_list[[10]]
write.table(df,"/results_ranking/all_datasets_gene_list_BA9.txt",quote=F,row.names = F,col.names = T,sep=";")

# GO term analysis 
  # BP
GO <- compareCluster(geneClusters = BA9_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pvalueCutoff=0.05,ont="BP")

GO<-pairwise_termsim(GO)

p1 <- emapplot(GO,legend_n=3,cex_line=0.1,cex_label_category=0.75,layout="nicely",cex_category=0.8,showCategory=20,pie="equal")+scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00" ,"#CAB2D6", "#6A3D9A"))
ggsave(paste0("/results_ranking/GO_BP_BA9_multiome.pdf"),p1,height=10,width=10)
ggsave(paste0("/results_ranking/GO_BP_BA9_multiome_legend.pdf"),p1,height=25,width=25)

tab <- table(GO@compareClusterResult$Description)
tab <- tab[order(tab, decreasing=T)]
write.table(tab ,"/results_ranking/high_score_pathways_GO_BP.txt",quote=F,row.names = F,col.names = F)

 