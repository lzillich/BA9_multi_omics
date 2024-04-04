# Author: Eric Zillich, last change: 2024-04-04
# Overlap between DE results in discovery and replication cohorts

# DBCBB discovery
DE_results <- read.csv("/DE_results/DE_results.txt")
DE_results_sig <- DE_results[DE_results$pvalue <0.05,]
DE_results_sig20 <- DE_results_sig$Gene[1:20]

# BA46 replication
DE_results_Ribeiro <- read.csv("/BA46_rep_data/DE_results.txt")
DE_results_Ribeiro_sig <- DE_results_Ribeiro[DE_results_Ribeiro$pvalue <0.05,]

# BA9 replication
BA9_rep_data <- read.csv("/BA9_rep_data/CUD.11mar2024.rep.14S.csv")
BA9_rep_data_sig <- BA9_rep_data[BA9_rep_data$pvalue<0.05,]

# Venn Diagram 

library(VennDiagram)
venn.diagram(list(BA9_discovery = DE_results_sig$Gene,BA9_replication=BA9_rep_data_sig$hgnc_symbol,BA46 = DE_results_Ribeiro_sig$Gene),filename="/replication/DE_genes_overlap.png", imagetype = "png",fill=c("#2b8cbe","#a6bddb","#ece7f2"),width=6,height=6,units="in",resolution=600)

# without labels
venn.diagram(list(BA9_discovery = DE_results_sig$Gene,BA9_replication=BA9_rep_data_sig$hgnc_symbol,BA46 = DE_results_Ribeiro_sig$Gene),filename="/replication/DE_genes_overlap_no_lab.png", imagetype = "png",fill=c("#2b8cbe","#a6bddb","#ece7f2"),width=4,height=4,units="in",resolution=600,category.names = c("","",""))

intersect(intersect( DE_results_sig$Gene,DE_results_Ribeiro_sig$Gene),BA9_rep_data_sig$hgnc_symbol)# "HSPA6" "FKBP4"

## RRHO for BA46 replication sample
# RRHO code based on documentation in http://htmlpreview.github.io/?https://github.com/RRHO2/RRHO2/blob/master/vignettes/RRHO2.html
library(RRHO2)
DE_results <- DE_results[order(abs(DE_results$log2FoldChange),decreasing = T),]
DE_results$RRHO_score <- -log10(DE_results$pvalue)*sign(DE_results$log2FoldChange)
DE_results_Ribeiro <- DE_results_Ribeiro[order(abs(DE_results_Ribeiro$log2FoldChange),decreasing = T),]
DE_results_Ribeiro$RRHO_score <- -log10(DE_results_Ribeiro$pvalue)*sign(DE_results_Ribeiro$log2FoldChange)

genes <- intersect(DE_results$Gene,DE_results_Ribeiro$Gene)
DE_results<- DE_results[DE_results$Gene %in% genes, c("Gene","RRHO_score")]
DE_results_Ribeiro <- DE_results_Ribeiro[DE_results_Ribeiro$Gene %in% genes, c("Gene","RRHO_score")]
RRHO_obj <-  RRHO2_initialize(DE_results, DE_results_Ribeiro , labels = c("BA9 discovery", "BA46 replication"), log10.ind=TRUE,method = "hyper")

pdf("/replication/RRHO_DBCBB_BA46_rep.pdf",height=6,width=6)
RRHO2_heatmap(RRHO_obj)
dev.off()

# RRHO for BA9 replication sample
DE_results <- read.csv("/BA46_rep_data/DE_results.txt")
DE_results <- DE_results[order(abs(DE_results$log2FoldChange),decreasing = T),]
DE_results$RRHO_score <- -log10(DE_results$pvalue)*sign(DE_results$log2FoldChange)
BA9_rep_data <- BA9_rep_data[order(abs(BA9_rep_data$log2FoldChange),decreasing = T),]
BA9_rep_data$RRHO_score <- -log10(BA9_rep_data$pvalue)*sign(BA9_rep_data$log2FoldChange)

genes <- intersect(DE_results$Gene,BA9_rep_data$hgnc_symbol)
DE_results<- DE_results[DE_results$Gene %in% genes, c("Gene","RRHO_score")]
colnames(BA9_rep_data)[10]<-"Gene"
BA9_rep_data <- BA9_rep_data[BA9_rep_data$Gene %in% genes, c("Gene","RRHO_score")]
RRHO_obj <-  RRHO2_initialize(DE_results, BA9_rep_data , labels = c("BA9 discovery", "BA9 replication"), log10.ind=TRUE,method = "hyper")

pdf("/replication/RRHO_DBCBB_BA9_rep.pdf",height=6,width=6)
RRHO2_heatmap(RRHO_obj)
dev.off()

