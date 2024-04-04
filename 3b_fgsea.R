# GSEA using clusterProfiler based on documentation from https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
# Author: Eric Zillich, last change: 2024-04-04

library(fgsea)
library(msigdbr)
library(ggplot2)

DE_results <- read.csv("/results_dir/DE_results.txt")
DE_results <- DE_results[complete.cases(DE_results),]
DE_results <- DE_results[!duplicated(DE_results$Gene),]


# Generate gene lists

DE_list <- DE_results$stat
names(DE_list) <- DE_results$Gene
DE_list <- sort(DE_list,decreasing = TRUE)


# GSEA with clusterProfiler
  
library(clusterProfiler)
library(enrichplot)
  library(DOSE)
  library(org.Hs.eg.db)


gse <- gseGO(
  DE_list,
  ont = "BP",
  OrgDb=org.Hs.eg.db,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  by = "fgsea")

gse<- pairwise_termsim(gse)
gse_pos <- gse@result[gse@result$NES > 0,]
gse_pos2 <- gse_pos[order(gse_pos$NES,decreasing=T),]
gse_pos2 <- gse_pos2[gse_pos2$p.adjust<0.05,]
gse_neg <- gse@result[gse@result$NES < 0,]
gse_neg2 <- gse_neg[order(gse_neg$NES,decreasing=F),]
gse_neg2 <- gse_neg2[gse_neg2$p.adjust<0.05,]
write.table(gse_pos2,"/DE_downstream/fgsea/BP_up_gsea.txt",sep=";",quote=F,row.names = F,col.names = T)
write.table(gse_neg2,"/DE_downstream/fgsea/BP_down_gsea.txt",sep=";",quote=F,row.names = F,col.names = T)

gse_up <- gse
gse_up@result <- gse_pos

gse_down <- gse
gse_down@result <- gse_neg


p1 <- emapplot(gse_up,showCategory=15,cex.params = list(line = 0.1,category_label = 0.8,category_node = 0.8))+scale_fill_gradient(low="#E43F3F",high="gray60",name="p.adjust")

p2 <- emapplot(gse_down,showCategory=15,cex.params = list(line = 0.1,category_label = 0.8,category_node = 0.8))+scale_fill_gradient(low="#268989",high="gray60",name="p.adjust")

ggsave("DE_downstream/fgsea/GO_GSEA_BA9_upregulated_NES.pdf", plot=p1, width = 7, height = 6, units = "in")
ggsave("DE_downstream/fgsea/GO_GSEA_BA9_downregulated_NES.pdf", plot=p2, width = 7, height = 6, units = "in")
