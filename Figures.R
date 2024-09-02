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

#### Figure 1 ####

# 1a 
# Volcano plot
DE_results <- read.csv("/results_dir/DE_results.txt")
DE_results$DE <- 
  with(DE_results, ifelse(padj < 0.05, "DE 5% FDR",
                          ifelse(pvalue < 0.05 & log2FoldChange > 0.07, "DE up p<0.05",
                                 ifelse(pvalue < 0.05 & log2FoldChange < -0.07, "DE down p<0.05", "n.s."))))
DE_results$DE <- 
  factor(DE_results$DE, 
         ordered = TRUE, 
         levels = c("DE 5% FDR","DE up p<0.05","DE down p<0.05" ,"n.s."))

customPlot <- list(
  theme_minimal(base_size = 12), 
  scale_fill_manual(values=c("#7C0000","#E43F3F","#268989","gray80")), 
  scale_colour_manual(values=c("#7C0000","#E43F3F","#268989","gray80"))
)


pdf("/results_dir/DE_volcano_manual.pdf", width = 10, height = 8)
ggplot(data = DE_results, aes(log2FoldChange, -log10(pvalue), colour = DE)) +
  geom_point(size=1) + geom_hline(aes(yintercept = -log10(0.05)),size=0.2) + geom_hline(aes(yintercept = -log10(2.340386e-06)),size=0.2,linetype="dashed") + geom_vline(aes(xintercept=0),size=0.2)+
  geom_text(aes(label = Gene), 
            data = subset(DE_results, DE != "n.s." | abs(log2FoldChange) >= 0.5 & pvalue < 0.05), 
            vjust = 0, nudge_y = 0.1, size = 5, check_overlap = TRUE) +
  xlab("log2FoldChange") +
  customPlot+ theme(legend.position = "none",text=element_text(size=16))
dev.off()


# 1b
# Gene Overlap for marker genes used in CIBERSORT
DE_results <- read.csv("/results_dir/DE_results.txt")
markers10 <-read_excel("/ref_dir/Yu_He_markers_10fold.xls")

astro <- markers10$geneSymbol[markers10$enrichedCellType == "astrocytes"]
neuron<-markers10$geneSymbol[markers10$enrichedCellType == "neurons"] 
endo<-markers10$geneSymbol[markers10$enrichedCellType == "endothelial"]
oligo<-markers10$geneSymbol[markers10$enrichedCellType == "oligodendrocytes"]
micro<-markers10$geneSymbol[markers10$enrichedCellType == "microglia"]
OPC <- markers10$geneSymbol[markers10$enrichedCellType == "OPC"]


DE_up <- DE_results$Gene[DE_results$pvalue<0.05 & DE_results$log2FoldChange > 0]
DE_down<-DE_results$Gene[DE_results$pvalue<0.05 & DE_results$log2FoldChange < 0]

list <- list(astrocytes = astro, neurons = neuron,endothelial=endo,oligodendrocytes=oligo,microglia=micro,OPC=OPC,DE_genes_up = DE_up, DE_genes_down=DE_down)

library(GeneOverlap)
gom.obj <- newGOM(list, genome.size = 21364)

pdf("/zi-flstorage/group_genepi/data/EP/SysMed/Cocaine/BA9_multiome/mRNA/3_DE_results/GeneOverlap_DE_celltypes_CUD_associated_bl.pdf", width = 14, height = 14)
drawHeatmap(gom.obj,cutoff = 0.05,note.col = "black")
dev.off()

pdf("/zi-flstorage/group_genepi/data/EP/SysMed/Cocaine/BA9_multiome/mRNA/3_DE_results/GeneOverlap_DE_celltypes_CUD_associated_wh.pdf", width = 14, height = 14)
drawHeatmap(gom.obj,cutoff = 0.05,note.col = "white")
dev.off()

# 1c

library(fgsea)
library(msigdbr)

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

gse_up <- gse
gse_up@result <- gse_pos

gse_down <- gse
gse_down@result <- gse_neg

p1 <- emapplot(gse_up,showCategory=15,cex.params = list(line = 0.1,category_label = 0.8,category_node = 0.8))+scale_fill_gradient(low="#E43F3F",high="gray60",name="p.adjust")
p2 <- emapplot(gse_down,showCategory=15,cex.params = list(line = 0.1,category_label = 0.8,category_node = 0.8))+scale_fill_gradient(low="#268989",high="gray60",name="p.adjust")
ggsave("4_DE_downstream/fgsea/GO_GSEA_BA9_downregulated_NES.pdf", plot=p2, width = 7, height = 6, units = "in")

# 1d
ggsave("4_DE_downstream/fgsea/GO_GSEA_BA9_upregulated_NES.pdf", plot=p1, width = 7, height = 6, units = "in")


#### Figure 2 ####

#2a
# Splicing volcano plot
load("/Splicing/leafviz/leafviz.RData")
clusters <- merge(clusters,introns,by="clusterID")
clusters$gene.x <- gsub("<i>","",clusters$gene.x)
clusters$gene.x <- gsub("</i>","",clusters$gene.x)

clusters$dPSI <- 
  with(clusters, ifelse(FDR < 0.05 & deltapsi > 0.025, "dPSI up 5% FDR",
                        ifelse(FDR < 0.05 & deltapsi < -0.025, "dPSI down 5% FDR","n.s.")))
clusters$dPSI <- 
  factor(clusters$dPSI, 
         ordered = TRUE, 
         levels = c("dPSI up 5% FDR","dPSI down 5% FDR","n.s."))

customPlot <- list(
  theme_minimal(base_size = 12), 
  scale_fill_manual(values=c("#E43F3F","#268989","gray80")), 
  scale_colour_manual(values=c("#E43F3F","#268989","gray80"))
)

p1 <- ggplot(data = clusters, aes(deltapsi, -log10(FDR), colour = dPSI)) +
  geom_point(size=1) + geom_hline(aes(yintercept = -log10(0.05)),size=0.2) +geom_vline(aes(xintercept=0),size=0.2) + xlab("dPSI")  + ylab("-log10(FDR q-value)") + geom_text(aes(label = gene.x), data = subset(clusters, dPSI != "n.s." | abs(deltapsi) >= 0.025 & FDR < 0.05), 
                                                                                                                                                                             vjust = 0, nudge_y = 0.1, size = 4, check_overlap = TRUE) +
  customPlot+ theme(legend.position = "none",text=element_text(size=16))
ggsave("/Splicing/Splicing_volcano_manual.pdf",p1, width = 6, height = 8)


#2b 
# Generated using the leafCutter Shiny App: $ ./run_leafviz.R leafviz.RData in command line
# Then entered BIN1 gene, download as Figure is implemented in the leafviz.R Shiny app

#2c
load("/Splicing/leafviz/leafviz.RData")
clusters <- merge(introns,clusters, by="clusterID")
clusters <- clusters[clusters$FDR <0.05 & abs(clusters$deltapsi) > 0.025,]
clusters$gene.x <- gsub("<i>","",clusters$gene.x)
clusters$gene.x <- gsub("</i>","",clusters$gene.x)
genes <- unique(clusters$gene.x[clusters$gene.x != "." & clusters$FDR < 0.05])

# GO enrichment analysis
BA9_AS <- enrichGO(genes,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",readable = T,ont="ALL",pvalueCutoff = 0.25)
BA9_AS <-pairwise_termsim(BA9_AS)

enr_BA9_AS <- emapplot(BA9_AS,showCategory=20,cex_line=0.1,cex_label_category=1.2,layout="nicely",cex_category=0.8,label_format = 30)+theme(legend.position = "none")+scale_fill_gradient(low="#b2df8a",high="gray60",name="p.adjust")
ggsave("/Splicing/All_AS_genes_emap.pdf",enr_BA9_AS,width=10,height=10)

#2d
# Venn Diagram overlap DE and AS genes
DE_results <- read.csv("/DE_results/DE_results.txt")
DE_results <- DE_results[DE_results$pvalue <0.05,]
DE_genes <- unique(DE_results$Gene)
intersect(DE_genes,genes) 
# "ITPKB"  "CPLX1"  "HLA-F"  "INPP5E" "GALNT8" "IGFBP6" "ZBTB4"  "BCAT2" 
library(VennDiagram)
venn.diagram(list(DE=DE_genes, AS = genes),filename="/Splicing/DE_genes_AS_genes_overlap.png", imagetype = "png",fill=c("gray60","#b2df8a"),width=3.5,height=3.5,units="in",resolution=600,category.names =c("",""))

#### Figure 3 ####

#3a
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
venn.diagram(list(BA9_discovery = DE_results_sig$Gene,BA9_replication=BA9_rep_data_sig$hgnc_symbol,BA46 = DE_results_Ribeiro_sig$Gene),filename="/replication/DE_genes_overlap.png", imagetype = "png",fill=c("#2b8cbe","#a6bddb","#ece7f2"),width=6,height=6,units="in",resolution=600)
intersect(intersect( DE_results_sig$Gene,DE_results_Ribeiro_sig$Gene),BA9_rep_data_sig$hgnc_symbol)# "HSPA6" "FKBP4"

#3b
# Volcano plot for spliceosome-associated genes 
library(KEGGREST)
# Spliceosome
genes_spli <- keggGet("hsa03040")[[1]]$GENE
genes_spli <-  genes_spli[seq(0,length(genes_spli),2)]
genes_spli <- gsub("\\;.*","",genes_spli)
genes_spli <- unique(genes_spli)

DE_results$Spliceosome <- 
  with(DE_results, ifelse(pvalue < 0.05 & log2FoldChange > 0.07 & Gene %in% genes_spli, "Spliceosome gene up p<0.05",
                          ifelse(pvalue < 0.05 & log2FoldChange < -0.07 & Gene %in% genes_spli, "Spliceosome gene down p<0.05", 
                                 ifelse(pvalue > 0.05 & Gene %in% genes_spli,"Spliceosome gene","other"))))
DE_results$Spliceosome <- 
  factor(DE_results$Spliceosome, 
         ordered = TRUE, 
         levels = c("Spliceosome gene up p<0.05","Spliceosome gene down p<0.05","Spliceosome gene" ,"other"))

customPlot <- list(
  theme_minimal(base_size = 12), 
  scale_fill_manual(values=c("#E43F3F","#268989","black","gray80")), 
  scale_colour_manual(values=c("#E43F3F","#268989","black","gray80"))
)

DE_2 <- DE_results[DE_results$Spliceosome %in% c("Spliceosome gene up p<0.05","Spliceosome gene down p<0.05","Spliceosome gene"),]
DE_1 <- DE_results[!(DE_results$Spliceosome %in% c("Spliceosome gene up p<0.05","Spliceosome gene down p<0.05","Spliceosome gene")),]
DE_3 <- rbind(DE_1,DE_2)

pdf("/DE_results/DE_volcano_spliceosome.pdf", width = 5, height = 5)
ggplot(data = DE_3, aes(log2FoldChange, -log10(pvalue), colour = Spliceosome)) +
  geom_point(size=1) + geom_hline(aes(yintercept = -log10(0.05)),size=0.2) + geom_hline(aes(yintercept = -log10(2.340386e-06)),size=0.2,linetype="dashed") + geom_vline(aes(xintercept=0),size=0.2)+
  geom_text_repel(aes(label = Gene), 
                  data = subset(DE_results,Gene %in% genes_spli & Spliceosome != "Spliceosome gene"), 
                  vjust = 0, nudge_y = 0.1, size = 3,max.overlaps=10) +
  xlab("log2FoldChange") +
  customPlot +theme(legend.position = "none",text=element_text(size=12))
dev.off()


#3c - manually generated table 

#3d
library(RRHO2)
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

#3e 
# RRHO for BA9 replication sample
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

#### Figure 4 ####

#4a
# Venn Diagram Methylation-Expression-Splicing 
DE_results <- read.csv("/DE_results/DE_results.txt")
DE_results <- DE_results$Gene[DE_results$pvalue <0.05]
EWAS_results <- read.csv("/EWAS/results/CUD_EWAS_results_20220920_lm.txt", sep="")
EWAS_results <- EWAS_results$Gene[EWAS_results$P_VAL <0.01]
EWAS_results <- unique(unlist(str_split(EWAS_results,pattern=";")))[-1]

# Splicing results 
load("/Splicing/leafviz/leafviz.RData")
clusters <- merge(introns,clusters, by="clusterID")
clusters <- clusters[clusters$FDR <0.05 & abs(clusters$deltapsi) > 0.025,]
clusters$gene.x <- gsub("<i>","",clusters$gene.x)
clusters$gene.x <- gsub("</i>","",clusters$gene.x)
AS_genes <- unique(clusters$gene.x[clusters$gene.x != "." & clusters$FDR < 0.05])
intersect(intersect(EWAS_results, DE_results),AS_genes) # "ZBTB4"  "INPP5E"

library(VennDiagram)
venn.diagram(list(DE = DE_results, AS = AS_genes, Meth = EWAS_results),filename="/Splicing/DE_genes_DiffMeth_genes_AS_genes_overlap.png", imagetype = "png",fill=c("#1f78bf","#b2df8a","#a6cee3"),width=4.5,height=4.5,units="in",resolution=600,category.names = c("","",""))

#4b
library(stringr)
# Create correlation plot for DE genes that also show differential methylation 
DE_results2 <- read.csv("/DE_results/DE_results.txt")
DE_results2 <- DE_results2[DE_results2$pvalue <0.05 & DE_results2$Gene %in% diffmethdiffexp,]

EWAS_results2 <- read.csv("/EWAS/results/CUD_EWAS_results_20220920_lm.txt", sep="")
EWAS_results2 <- EWAS_results2[EWAS_results2$P_VAL <0.01,]
EWAS_results2 <- EWAS_results2[complete.cases(EWAS_results2$Gene),]
genenames <- str_split(EWAS_results2$Gene,pattern=";")

library(tidyverse)
EWAS_results2 <- EWAS_results2  %>%
  add_column(gene1 = "",
             gene2 = "",gene3 = "",
             gene4 = "",gene5 = "",
             gene6 = "",gene7 = "",
             gene8 = "",gene9 = "",
             gene10 = "",gene11 = "",
             gene12 = "",gene13= "",
             gene14 = "",gene15 = "",
             gene16 = "",gene17 = "",
             gene18 = "",gene19 = "",
             gene20 = "",gene21 = "",
             gene22 = "" )


for(i in 1:length(genenames)){
  genenames[[i]] <- unique(genenames[[i]])
  EWAS_results2[i,12] <- genenames[[i]][1]
  if(length(genenames[[i]]) > 1){
    EWAS_results2[i,13] <- genenames[[i]][2]
  }
  if(length(genenames[[i]]) > 2){
    EWAS_results2[i,14] <- genenames[[i]][3]
  }
  if(length(genenames[[i]]) >3){
    EWAS_results2[i,15] <- genenames[[i]][4]
  }
  if(length(genenames[[i]]) > 4){
    EWAS_results2[i,16] <- genenames[[i]][5]
  }
  if(length(genenames[[i]]) > 5){
    EWAS_results2[i,17] <- genenames[[i]][6]
  }
  if(length(genenames[[i]]) > 6){
    EWAS_results2[i,18] <- genenames[[i]][7]
  }
  if(length(genenames[[i]]) > 7){
    EWAS_results2[i,19] <- genenames[[i]][8]
  }
  if(length(genenames[[i]]) > 8){
    EWAS_results2[i,20] <- genenames[[i]][9]
  }
  if(length(genenames[[i]]) > 9){
    EWAS_results2[i,21] <- genenames[[i]][10]
  }
  if(length(genenames[[i]]) > 10){
    EWAS_results2[i,22] <- genenames[[i]][11]
  }
  if(length(genenames[[i]]) > 11){
    EWAS_results2[i,23] <- genenames[[i]][12]
  }
  if(length(genenames[[i]]) > 12){
    EWAS_results2[i,24] <- genenames[[i]][13]
  }
  if(length(genenames[[i]]) > 13){
    EWAS_results2[i,25] <- genenames[[i]][14]
  }
  if(length(genenames[[i]]) > 14){
    EWAS_results2[i,26] <- genenames[[i]][15]
  }
  if(length(genenames[[i]]) > 15){
    EWAS_results2[i,27] <- genenames[[i]][16]
  }
  if(length(genenames[[i]]) > 16){
    EWAS_results2[i,28] <- genenames[[i]][17]
  }
  if(length(genenames[[i]]) > 17){
    EWAS_results2[i,29] <- genenames[[i]][18]
  }
  if(length(genenames[[i]]) > 18){
    EWAS_results2[i,30] <- genenames[[i]][19]
  }
  if(length(genenames[[i]]) > 19){
    EWAS_results2[i,31] <- genenames[[i]][20]
  }
  if(length(genenames[[i]]) > 20){
    EWAS_results2[i,32] <- genenames[[i]][21]
  }
  if(length(genenames[[i]]) > 21){
    EWAS_results2[i,33] <- genenames[[i]][22]
  }
  
}


# Create results df with CpG-associated genes 
mat <- matrix(nrow=0,ncol=33,NA)
colnames(mat) <- colnames(EWAS_results2)
for(i in diffmethdiffexp){
  res <- EWAS_results2[EWAS_results2$gene1 == i  | EWAS_results2$gene2 == i | EWAS_results2$gene3 == i | EWAS_results2$gene4 == i |  EWAS_results2$gene5 == i| EWAS_results2$gene6 == i |  EWAS_results2$gene7 == i|  EWAS_results2$gene8 == i|
                         EWAS_results2$gene9 == i|EWAS_results2$gene10 == i  |EWAS_results2$gene11 == i  |EWAS_results2$gene12 == i  |EWAS_results2$gene13 == i  |EWAS_results2$gene14 == i  |EWAS_results2$gene15 == i  | EWAS_results2$gene16 == i |EWAS_results2$gene17 == i  |EWAS_results2$gene18 == i  |EWAS_results2$gene19 == i  |EWAS_results2$gene20 == i| EWAS_results2$gene21 == i |EWAS_results2$gene22 == i, ]
  nrow(res)
  mat <- rbind(mat,res)
}

# non-unique genes:keep CpG site with strongest significance 

mat2 <- mat[!(mat$cg %in% c("cg23347250","cg11903404","cg22598810","cg11863717","cg00999152","cg01985396","cg16543958","cg19894264","cg12476579","cg09579928","cg13137321","cg02091100","cg19214184","cg14453704","cg24380053","	
cg10764626","cg16691888","cg14572252","cg04367164","cg24849439","	
cg08747692","cg26848184","cg01189092"	,"cg17962547","cg10764626","cg08747692","cg10426483","cg01124988",
                            "cg10426483","cg25011720","	
cg01189092","cg08607710","cg08607710","cg21893960","cg09033472","cg26804406","cg22985016","cg15384365","cg05496603","cg15420683","cg04221833","cg17612239","cg05479968","cg24392197","cg21796322","cg06125962","cg07491858","")),]

# rename genes with multiple names so that gene_final is the right label
mat2$gene_final <- mat2$gene1
mat2$gene_final[mat2$gene2 %in% diffmethdiffexp] <- mat2$gene2[mat2$gene2 %in% diffmethdiffexp]
mat2$gene_final[mat2$gene_final == "PCDHGA1"] <- "PCDHGB3"
mat2$gene_final %in% diffmethdiffexp #TRUE
mat3 <- mat2[,c("gene_final","BETA","P_VAL")]


DE_results2 <- DE_results2[,c("Gene","log2FoldChange","pvalue")]
DE_results2$Gene %in% mat3$gene_final # TRUE
colnames(DE_results2)[1]<-"gene_final"

#merge methylation and expression 
methexpr <- merge(DE_results2,mat3,by="gene_final")

methexpr$corr <- 
  with(methexpr, ifelse(log2FoldChange < 0 & BETA > 0, "methylation up - expression down",
                        ifelse(log2FoldChange > 0 & BETA > 0, "methylation up - expression up", 
                               ifelse(log2FoldChange < 0 & BETA < 0,"methylation down - expression down","methylation down- expression up"))))
methexpr$corr <- 
  factor(methexpr$corr, 
         ordered = TRUE, 
         levels = c("methylation down- expression up","methylation up - expression down","methylation down - expression down" ,"methylation up - expression up"))

customPlot <- list(
  theme_minimal(base_size = 12), 
  scale_fill_manual(values=c("#E43F3F","#268989","gray20","gray60")), 
  scale_colour_manual(values=c("#E43F3F","#268989","gray20","gray60"))
)

library(ggrepel)
pdf("/DE_results/DE_DiffMeth_pabel.pdf", width = 9, height = 6)
ggplot(methexpr , aes(x=BETA, y=log2FoldChange,colour=corr)) + geom_point(size=1) +xlim(-5.5,5.5)+ylim(-5.5,5.5)+ xlab("Differential Methylation (effect size)")+ylab("Differential Expression (log2FC)")+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+geom_text_repel(aes(label = gene_final), vjust = 0, nudge_y = 0.1, size = 3, max.overlaps=10)+customPlot+geom_hline(aes(yintercept = 0),size=0.1) + geom_vline(aes(xintercept = 0),size=0.1) +theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()


#4c - created using SparK (Kurtenbach & Harbour, bioRxiv, 2019) as indicated in the methods part of the manuscript

#### Supplementary Figures ####

#### Figure S1 ####

# S1a
# variance partition 
pheno <- read.csv("/dir/pheno.txt", sep=";")
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

# S1b
# celltype proportion plot based on CIBERSORT estimates
results10 <- read.csv("/CIBERSORT/Markers_10fold_Yu_He_BA9.txt", sep=";")
res10_2 <- results10[,c(1:7)]
colnames(res10_2) <- gsub("fpkm_","",colnames(res10_2))
colnames(res10_2)[1]<-"rn"

# merge with pheno data to get the sample ID infos 
pheno <- read.csv("/DE_results/pheno.txt", sep=";")
merged <- merge(pheno,res10_2,by="rn")
merged <- merged[order(merged$CUD,decreasing=T),]
merged$name <- c(paste0("CUD-",c(1:13)),paste0("Ctrl-",c(1:12)))
rownames(merged) <- merged$name
merged <- merged[,c(8:13)] # subset for celltypes only

perc_df <- data.frame(samples=rep(rownames(merged),each=6),celltypes=rep(colnames(merged),times=25),value=as.vector(t(merged)))
perc_df$samples <- factor(perc_df$samples,levels=rownames(merged))

# Plot the fraction of celltypes 
per <- ggplot(perc_df, aes(fill=celltypes, x=value, y=samples)) + 
  geom_bar(position="fill", stat="identity")+scale_y_discrete(limits=rev(rownames(merged)))+scale_fill_manual(values=c("#F7AB64","#74BA59","#70305A","#E8326D", "#3A9BCC","#85CEE4"))+xlab("celltype proportion")+ylab(NULL)+
  theme_minimal()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(plot = per,filename = "/CIBERSORT/celltype_percentage.pdf", width = 6, height=6)

# S1c
DE_results <- read.csv("/DE_results/DE_results.txt")
# qqplot
png("results_dir/DE_qqplot.png", width = 800, height = 800)
qq(DE_results$pvalue)
dev.off()

# S1d
# correlation plot
png("/with_AUD_MDD/sensitivity_MDD_AUD.png",height=6,width=6,res=600,units="in")
plot(merged$log2FoldChange.x, merged$log2FoldChange.y, xlab = "log2FC - w/o AUD MDD", ylab = "log2FC - with AUD MDD",xlim=c(-5,5),ylim=c(-5,5))
dev.off()
cor.test(merged$log2FoldChange.x, merged$log2FoldChange.y)# r=0.94, p-value < 2.2e-16


#### Figure S2 ####
library(WGCNA)
library(biomaRt)
DF <- data.frame
setwd(dir0)
options(stringsAsFactors = FALSE)
#read pheno
pheno <- read.csv("/pheno_dir/pheno.txt", sep=";")
pheno$pH[is.na(pheno$pH)] <- 6.386957  #mean impute
options(stringsAsFactors = FALSE);

# Load the expression and trait data saved in the first part
setwd("/WGCNA_dir/")
lnames = load(file ="BA9_Expr_dataInput_processed.RData") 
# Load network data saved in the second part.
lnames = load(file = "/dir2/BA9_networkConstruction_Expr.RData");
# Define numbers of genes and samples
setwd(dir0)
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

# Star visualization in the heatmap
stars_gen <- function(x){
  stars <- c( "***", "**", "*","")
  vec <- c(0, 0.001, 0.01, 0.05,1.01)
  i <- findInterval(x, vec)
  stars[i]
}

textmat2 <- moduleTraitPvalue

for(i in 1:length(colnames(textmat2))){
  textmat2[,i] <- as.character(stars_gen(textmat2[,i]))
}

#S2a
pdf("/WGCNA_dir/Heatmap_modules_Expr_filtered_stars.pdf", width = 5, height = 10)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor[,c(1:7)],
               xLabels = names(datTraits)[c(1:7)],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textmat2[,c(1:7)],
               setStdMargins = FALSE,
               cex.text = 1.2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#S2b
# S2b was generated during network construction - see script 5b for more details, code below in comments, takes long computation time 

# setwd("WGCNA_dir/")
# lnames = load(file ="BA9_Expr_dataInput_processed.RData")    
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# 
# sft$powerEstimate
# net = blockwiseModules(datExpr, power = sft$powerEstimate,
#                        TOMType = "signed", minModuleSize = 10, maxBlockSize = 36000,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "/dir2/BA9_mRNA_TOM",
#                        verbose = 3)
# table(net$colors)
# setwd("/dir2/")

# Convert labels to colors for plotting
#mergedColors = labels2colors(net$colors)

# Plot the dendrogram
png("BA9_dendogram_Expr.png", width = 1200, height = 1800)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#S2c
# Define variable weight containing the weight column of datTrait
CUD = as.data.frame(datTraits$CUD)
names(CUD) = "CUD"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, CUD, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(CUD), sep="");
names(GSPvalue) = paste("p.GS.", names(CUD), sep="")

mod <- DF(moduleTraitCor)

# GS-MM plot for module yellow
module <- "yellow"

for(i in module){
  column = match(i, modNames);
  moduleGenes = moduleColors==i;
  
  
  png(paste0("/WGCNA_dir/module_membership_gene_significance_mRNA_",i,".png",sep=""),height=6,width=6,res=600,units="in")
  par(mar = c(4,4,4,4));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", i, "module"),
                     ylab = "Gene significance for CUD status",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2 ,col = i)
  dev.off()
}
  

#S2d
probes = colnames(datExpr)
geneInfo0 = data.frame(genx = probes,
                       geneSymbol =probes,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
GO_gene_list <- geneInfo0$geneSymbol[geneInfo0$moduleColor == "yellow"]

# GO enrichment analysis
module_GO_BP <- enrichGO(GO_gene_list,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pvalueCutoff=1,ont="BP")
module_GO_BP <-pairwise_termsim(module_GO_BP)

write.table(module_GO_BP@result[module_GO_BP@result$p.adjust<0.05,],"/WGCNA_dir/GO_statistics_yellow.txt", sep = ";", quote = F, row.names = F)

p1 <- emapplot(module_GO_BP,legend_n=3,cex_line=0.1,cex_label_category=0.75,layout="nicely",cex_category=0.8,showCategory=30)
ggsave("/WGCNA_dir/GO_modules_BP.pdf",p1,height=9,width=9)


#### Figure S3 ####
#Created using Cytoscape STRING as indicated in the methods part of the manuscript

#### Figure S4 ####
library(data.table)
library(variancePartition)
library(readxl)
library(BiocParallel)

#S4a
#read pheno
pheno <- read.csv("/dir/pheno.txt", sep=";")
pheno$pH[is.na(pheno$pH)] <- 6.335

# Add cibersort cell type estimates
ciber <- read.csv("/CIBERSORT/Markers_10fold_Yu_He_BA9.txt", sep=";")
ciber <- ciber[,c(1:7)]
colnames(ciber)[1] <- "rn"

pheno_cc <- merge(pheno,ciber,by="rn")
rownames(pheno_cc) <- pheno_cc$rn

pheno_cc <- pheno_cc[,c("PMI","pH","RIN","Age")]

#specify model
mRNA_CUD <- get(load("/dir//mRNA_CUD_counts.Rdata"))
colnames(mRNA_CUD$counts) <- gsub("Aligned.sortedByCoord.out.bam","",colnames(mRNA_CUD$counts))
mRNA_CUD <- as.data.frame(mRNA_CUD$counts)
keep <- rowSums(mRNA_CUD) >= 1
mRNA_CUD <- mRNA_CUD[keep,]

form <- ~ Age + pH + PMI + RIN 

varPart <- fitExtractVarPartModel(mRNA_CUD, form, pheno_cc,BPPARAM = MulticoreParam(20))
vp <- sortCols(varPart)

png("/dir/Variance_partition.png",width = 8,height=8,units = "in",res=300)
plotVarPart(vp,col=c( "#C77CFF","#00BFC4","#7CAE00","orange","grey85"))
dev.off()

#S4b
#### 1. Import phenotype data ####
Ribeiro_pheno <- read_excel("/dir/Ribeiro_pheno.xlsx")

DT <- data.table
DF <- data.frame

#### 2. Import counts ####

# CUD
mRNA_CUD <- get(load("/dir/Ribeiro_mRNA_CUD_counts.Rdata"))
colnames(mRNA_CUD$counts) <- gsub("_name_sorted.bam","",colnames(mRNA_CUD$counts))
mRNA_CUD$targets <- gsub("_name_sorted.bam","",mRNA_CUD$targets)
colnames(mRNA_CUD$stat) <- gsub("_name_sorted.bam","",colnames(mRNA_CUD$stat))

mRNA_CUD_counts <- DF(mRNA_CUD$counts)
pheno_cc <- Ribeiro_pheno[,c("PMI","pH","RIN","Age")]

#specify model
keep <- rowSums(mRNA_CUD_counts) >= 1
mRNA_CUD_counts <- mRNA_CUD_counts[keep,]

form <- ~ Age + pH + PMI + RIN

varPart <- fitExtractVarPartModel(mRNA_CUD_counts, form, pheno_cc,BPPARAM = MulticoreParam(20))
vp <- sortCols(varPart)

png("/dir/Variance_partition_Ribeiro.png",width = 8,height=8,units = "in",res=300)
plotVarPart(vp,col=c("orange" ,"#C77CFF","#7CAE00","#00BFC4","grey85"))
dev.off()

#### Figure S5 ####
#S5a - created using BioRender.com

#S5b - based on downloaded results from CMap (clue.io)
library(forcats)
# Individual perturbagens 
query_result <- read.delim("/CMap_BA9/arfs/TAG/query_result.gct")
query_result <-query_result[-1,] # remove descriptive header
#-log10(0.05) -> 1.30103

query_result_sig <-query_result[abs(as.numeric(query_result$fdr_q_nlog10)) > 1.30103,]

pos_corr_perturbagens <- query_result_sig[order(as.numeric(query_result_sig$norm_cs),decreasing = T),][c(1:20),]
neg_corr_perturbagens <- query_result_sig[order(as.numeric(query_result_sig$norm_cs),decreasing = F),][c(1:20),]

# waterfall plot for perturbagens, subset top 20 up and down results for visualization
neg_corr_perturbagens <- neg_corr_perturbagens[,c("pert_iname","fdr_q_nlog10","norm_cs")]
colnames(neg_corr_perturbagens) <- c("compound","pval","cscore_norm")
neg_corr_perturbagens$pval <- as.numeric(neg_corr_perturbagens$pval)
neg_corr_perturbagens$cscore_norm <- as.numeric(neg_corr_perturbagens$cscore_norm)

pos_corr_perturbagens <- pos_corr_perturbagens[,c("pert_iname","fdr_q_nlog10","norm_cs")]
colnames(pos_corr_perturbagens) <- c("compound","pval","cscore_norm")
pos_corr_perturbagens$pval <- as.numeric(pos_corr_perturbagens$pval)
pos_corr_perturbagens$cscore_norm <- as.numeric(pos_corr_perturbagens$cscore_norm)

pert <- rbind(pos_corr_perturbagens,neg_corr_perturbagens)
pert$fdr_sig <- ifelse(pert$pval > 1.30103, "yes", "no")
pert <- pert[!duplicated(pert$compound),]
pert <- pert[order(pert$cscore_norm,decreasing = T),]

p1 <- ggplot(pert, aes(reorder(compound,cscore_norm), cscore_norm)) +
  geom_col(aes(fill=fdr_sig))+ scale_fill_manual(values=c("#1f78bf","gray60")) +
  coord_flip() +
  labs(x="", y="Normalized connectivity score",
       title="perturbagens",legend="sig") + theme_minimal()
ggsave("/CMap_BA9/perturbagen.pdf",p1,height=6,width=8)

#S5c
gsea_result <- read.delim("/CMap_BA9/gsea/TAG/arfs/NORM_CS/gsea_result.gct")
gsea_result <- gsea_result[-1,]

gsea_result_sig_PCL <- gsea_result[gsea_result$set_type == "PCL",]
gsea_result_sig_PCL <- gsea_result_sig_PCL[order(gsea_result_sig_PCL$fdr_q_nlog10,decreasing=T),]

# PCL
pos_corr_gsea_PCL <- gsea_result_sig_PCL[order(as.numeric(gsea_result_sig_PCL$norm_cs),decreasing = T),][c(1:30),]
pos_corr_gsea_PCL$src_set_id <- gsub("^...","",pos_corr_gsea_PCL$src_set_id)
pos_corr_gsea_PCL$cell_iname[pos_corr_gsea_PCL$cell_iname=="-666"]<-""
pos_corr_gsea_PCL$name <- paste0(pos_corr_gsea_PCL$src_set_id,"_",pos_corr_gsea_PCL$cell_iname)
pos_corr_gsea_PCL$name <- gsub("*_$","",pos_corr_gsea_PCL$name)

neg_corr_gsea_PCL <- gsea_result_sig_PCL[order(as.numeric(gsea_result_sig_PCL$norm_cs),decreasing = F),][c(1:30),]
neg_corr_gsea_PCL$src_set_id <- gsub("^...","",neg_corr_gsea_PCL$src_set_id)
neg_corr_gsea_PCL$cell_iname[neg_corr_gsea_PCL$cell_iname=="-666"]<-""
neg_corr_gsea_PCL$name <- paste0(neg_corr_gsea_PCL$src_set_id,"_",neg_corr_gsea_PCL$cell_iname)
neg_corr_gsea_PCL$name <- gsub("*_$","",neg_corr_gsea_PCL$name)

# GSEA PCL

neg_corr_gsea_PCL <- neg_corr_gsea_PCL[,c("name","fdr_q_nlog10","norm_cs")]
colnames(neg_corr_gsea_PCL) <- c("PCL","pval","cscore_norm")
neg_corr_gsea_PCL$pval <- as.numeric(neg_corr_gsea_PCL$pval)
neg_corr_gsea_PCL$cscore_norm <- as.numeric(neg_corr_gsea_PCL$cscore_norm)

pos_corr_gsea_PCL <- pos_corr_gsea_PCL[,c("name","fdr_q_nlog10","norm_cs")]
colnames(pos_corr_gsea_PCL) <- c("PCL","pval","cscore_norm")
pos_corr_gsea_PCL$pval <- as.numeric(pos_corr_gsea_PCL$pval)
pos_corr_gsea_PCL$cscore_norm <- as.numeric(pos_corr_gsea_PCL$cscore_norm)

PCL <- rbind(pos_corr_gsea_PCL,neg_corr_gsea_PCL)
PCL$fdr_sig <- ifelse(PCL$pval > 1.30103, "yes", "no")
PCL <- PCL[!duplicated(PCL$PCL),]
PCL <- PCL[order(PCL$cscore_norm,decreasing = T),]
PCL <- PCL[c(1:20,41:60),]

p2 <-ggplot(PCL, aes(fct_reorder(PCL,cscore_norm), cscore_norm)) +
  geom_col(aes(fill=fdr_sig)) +
  coord_flip() + scale_fill_manual(values=c("gray60","#1f78bf")) + theme_minimal() +
  labs(x="", y="Normalized connectivity score",
       title="GSEA perturbagen class",legend="sig")
ggsave("/CMap_BA9/gsea_PCL.pdf",p2,height=6,width=8)

#S5d
# subset for glucocorticoid inside labels
glucocort_perturbagens <- query_result_sig[grep("Glucocor",query_result_sig$moa),]
glucocort_perturbagens <- glucocort_perturbagens[,c("moa","cell_iname","pert_idose","pert_itime","pert_iname","fdr_q_nlog10","norm_cs")]
glucocort_perturbagens$name <- paste0(glucocort_perturbagens$pert_iname,"_",glucocort_perturbagens$cell_iname,"_",glucocort_perturbagens$pert_idose,"_",glucocort_perturbagens$pert_itime)
glucocort_perturbagens$fdr_q_nlog10 <- as.numeric(glucocort_perturbagens$fdr_q_nlog10)
glucocort_perturbagens$norm_cs <- as.numeric(glucocort_perturbagens$norm_cs)
glucocort_perturbagens$fdr_sig <- ifelse(glucocort_perturbagens$fdr_q_nlog10 > 1.30103, "yes", "no")
glucocort_perturbagens <- glucocort_perturbagens[order(glucocort_perturbagens$norm_cs,decreasing = F),]

p3 <- ggplot(glucocort_perturbagens, aes(reorder(name,norm_cs), norm_cs)) +
  geom_col(aes(fill=fdr_sig))+ylim(-2,2) +scale_fill_manual(values=c("#1f78bf","gray60")) +
  coord_flip() +
  labs(x="", y="Normalized connectivity score",
       title="Glucocorticoid receptor targeting drugs",legend="sig") + theme_minimal()
ggsave("/CMap_BA9/glucocort_modulators.pdf",p3,height=6,width=8)

#### Figure S6 ####
library(MOFA2)
library(MOFAdata)
BA9 <- readRDS("/MOFA_dir/output/MOFA_BA9_trained_model_20240102.rds")
# metadata
pheno_meth_expr <- read.delim("/MOFA_dir/input/pheno_meth_expr.txt")
colnames(pheno_meth_expr)[1] <- "sample"
pheno_meth_expr$sample <- as.character(pheno_meth_expr$sample)
samples_metadata(BA9) <- pheno_meth_expr

#S6a 
p1 <- plot_variance_explained(BA9)
ggsave("/MOFA_dir/model_factors.pdf",p1, width = 4, height = 5)

#S6b
p2 <- correlate_factors_with_covariates(BA9, abs=F,
                                        covariates = colnames(pheno_meth_expr[c("RIN","CUD","pH","Age","PMI","Axis_1_Dependence","Simplified_Axis_1")]), 
                                        plot="log_pval"
)
ggsave("/MOFA_dir/model_cor_cov.pdf",p2, width = 3, height = 5)

#S6c
p3 <- plot_factor(BA9, 
                  factors = 9, 
                  color_by = "CUD",
                  add_violin = TRUE,
                  dodge = TRUE
)+ylim(-0.16,0.16)+scale_fill_manual(values=c("#268989","#E43F3F"))
ggsave("/MOFA_dir/f9_CUD.pdf",p3, width = 5, height = 4)


#S6d
p4 <- plot_top_weights(BA9,
                       view = "Meth",
                       factor = 9,
                       nfeatures = 10,     # Top number of features to highlight
                       scale = T           # Scale weights from -1 to 1
)
ggsave("/MOFA_dir/meth_topweights.pdf",p4, width = 5, height = 4)

#S6e
p5 <- plot_top_weights(BA9,
                       view = "Expr",
                       factor = 9,
                       nfeatures = 10,     # Top number of features to highlight
                       scale = T           # Scale weights from -1 to 1
)
ggsave("/MOFA_dir/expr_topweights.pdf",p5, width = 5, height = 4)

#S6f
# Create input matrix for GSEA using GO BP terms for human
library(msigdbr)
# Extract the gene sets from the current GO BP repository to get the genes involved in the pathways
genesets = msigdbr(species = "human", category = "C5", subcategory = "BP") 
library(SAMBAR)
# convert the current GO BP gmt file from msigdb to a binary matrix with gene and geneset as columns and rows
gs <- convertgmt("/MOFA_dir/input/c5.go.bp.v2023.1.Hs.symbols.gmt",unique(genesets$gene_symbol))

# GSEA on positive weights, with default options
res.positive <- run_enrichment(BA9, 
                               feature.sets = gs, 
                               view = "Expr",
                               sign = "positive"
)

# GSEA on negative weights, with default options
res.negative <- run_enrichment(BA9, 
                               feature.sets = gs, 
                               view = "Expr",
                               sign = "negative"
)

# Extract enrichment results
pos_gsea1 <- as.data.frame(res.positive$set.statistics)
pos_gsea1$TERM <- rownames(pos_gsea1)
pos_gsea1 <- pos_gsea1[,c("TERM","Factor9")]
colnames(pos_gsea1)[2] <- "GSEA set statistic"
pos_gsea2 <- as.data.frame(res.positive$pval.adj)
pos_gsea2$TERM <- rownames(pos_gsea2)
pos_gsea2 <- pos_gsea2[,c("TERM","Factor9")]
colnames(pos_gsea2)[2] <- "padj"
pos_gsea <- merge(pos_gsea1,pos_gsea2,by="TERM")
pos_gsea <-pos_gsea[order(pos_gsea$padj,decreasing = F),]

neg_gsea1 <- as.data.frame(res.negative$set.statistics)
neg_gsea1$TERM <- rownames(neg_gsea1)
neg_gsea1 <- neg_gsea1[,c("TERM","Factor9")]
colnames(neg_gsea1)[2] <- "GSEA set statistic"
neg_gsea2 <- as.data.frame(res.negative$pval.adj)
neg_gsea2$TERM <- rownames(neg_gsea2)
neg_gsea2 <- neg_gsea2[,c("TERM","Factor9")]
colnames(neg_gsea2)[2] <- "padj"
neg_gsea <- merge(neg_gsea1,neg_gsea2,by="TERM")
neg_gsea <-neg_gsea[order(neg_gsea$padj,decreasing = F),]

# Visualize in waterfall plot
pos_ES <- as.data.frame(res.positive$set.statistics)
pos_padj <- as.data.frame(res.positive$pval.adj)
neg_ES <- as.data.frame(res.negative$set.statistics)
neg_padj <- as.data.frame(res.negative$pval.adj)

# create dataframe for plotting for factor 9
f9_GSEA <- data.frame(TERM=rep("",times=30),ES=rep(0,times=30), p_adj=rep(0,times=30),weights=c(rep("positive weights",times=15),rep("negative weights",times=15)))

f9_GSEA$p_adj[1:15] <- pos_padj$Factor9[order(pos_padj$Factor9,decreasing = F)]
f9_GSEA$TERM[1:15] <- rownames(pos_padj[order(pos_padj$Factor9,decreasing = F),])[1:15]
f9_GSEA$p_adj[16:30] <- neg_padj$Factor9[order(neg_padj$Factor9,decreasing = F)]
f9_GSEA$TERM[16:30] <- rownames(neg_padj[order(neg_padj$Factor9,decreasing = F),])[1:15]

pos_terms <- pos_ES[rownames(pos_ES) %in% rownames(pos_padj[order(pos_padj$Factor9,decreasing = F),])[1:15],]
pos_terms <- pos_terms[match(rownames(pos_padj[order(pos_padj$Factor9,decreasing = F),])[1:15],row.names(pos_terms)),]
neg_terms <- neg_ES[rownames(neg_ES) %in% rownames(neg_padj[order(neg_padj$Factor9,decreasing = F),])[1:15],]
neg_terms <- neg_terms[match(rownames(neg_padj[order(neg_padj$Factor9,decreasing = F),])[1:15],row.names(neg_terms)),]

f9_GSEA$ES[1:15] <- pos_terms$Factor9
f9_GSEA$ES[16:30] <- -(neg_terms$Factor9)

# rename pathways and round pvals for plotting
f9_GSEA$TERM <-  gsub("GOBP_","",f9_GSEA$TERM)
f9_GSEA$TERM <-  gsub("_"," ",f9_GSEA$TERM)
f9_GSEA$TERM <-  tolower(f9_GSEA$TERM)
f9_GSEA$p_adj <- signif(f9_GSEA$p_adj, digits=3)
f9_GSEA$TERM[f9_GSEA$TERM == "atp synthesis coupled electron transport"] <- "ATP synthesis coupled electron transport"
f9_GSEA$TERM[f9_GSEA$TERM == "mitochondrial electron transport nadh to ubiquinone"] <- "mitochondrial electron transport NADH to ubiquinone"
f9_GSEA$sig <- rep("***",times=length(f9_GSEA$TERM))

library(Seurat)
f9_gsea_plot <- ggbarplot(f9_GSEA, x = "TERM", y = "ES",
                          fill = "weights",               # change fill color by cyl
                          color = "black",            # Set bar border colors to white
                          sort.val = "desc",          # Sort the value in dscending order
                          sort.by.groups = FALSE, 
                          lab.size=4)+ylim(-12.5,12.5)+NoLegend()+xlab(NULL)+scale_fill_manual(values=c("#268989","#E43F3F"))+ylab("GSEA set statistic")+ geom_text(aes(y=2*abs(ES)/ES,label = sig), vjust = 0.8,hjust=0.3,size=8)+coord_flip()+theme_minimal() + theme(text=element_text(size=20))

ggsave("/MOFA_dir/GSEA_factor9_waterfall.pdf",f9_gsea_plot,width=14,height=10)

#S6g
#Methylation GO plots
f9_Meth <- get_weights(BA9,views="Meth")
f9_Meth<-f9_Meth$Meth
f9_Meth <- f9_Meth[,9]
f9_Meth <- sort(f9_Meth,decreasing = TRUE)

quantile(f9_Meth,probs=c(0.025,0.25,0.5,0.75,0.975))
#   2.5%          25%          50%          75%        97.5% 
# -0.145604531 -0.046109412  0.002832227  0.053496722  0.152205286 

f9_Meth_top <- f9_Meth[f9_Meth > 0.152205286 | f9_Meth < -0.145604531]
names_meth <- gsub("*...........-.","",names(f9_Meth_top))

# prepare for missMethyl 

f9_Meth_pos <- f9_Meth[f9_Meth > 0.152205286]
f9_Meth_neg <- f9_Meth[f9_Meth < -0.145604531]

names_f9_Meth_pos <- gsub(" -.*","",names(f9_Meth_pos))
names_f9_Meth_neg <- gsub(" -.*","",names(f9_Meth_neg))
names_meth_all <- gsub(" -.*","",names(f9_Meth))
names_meth <- unique(unlist(str_split(names_meth,";")))

# GO using missMethyl
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(cowplot)

go_res_pos <- gometh(names_f9_Meth_pos, names_meth_all, "GO", "EPIC", T,sig.genes = T)
go_res_pos <- go_res_pos[order(go_res_pos$P.DE),]
go_res_neg <- gometh(names_f9_Meth_neg, names_meth_all, "GO", "EPIC", T,sig.genes = T)
go_res_neg <- go_res_neg[order(go_res_neg$P.DE),]

# plot top GO terms
up_table <- go_res_pos[,c("TERM","P.DE")][c(1:15),]
up_table$P.DE <- -log10(up_table$P.DE)
position1 <- rev(up_table$TERM)

down_table <- go_res_neg[,c("TERM","P.DE")][c(1:15),]
down_table$P.DE <- -log10(down_table$P.DE)
position2 <- rev(down_table$TERM)

#S6g
p6 <- ggplot(data = down_table, aes(x = P.DE, y = TERM)) + geom_point(size=3,color="#268989") +theme_bw() + ylab("") + xlab("-log10(p)") + scale_y_discrete(limits = position2)+xlim(2,3.7) 
ggsave("/MOFA_dir/output/GO_factor9_meth_neg_weights.pdf",p6,width=8,height=4)

#S6h
p7 <- ggplot(data = up_table, aes(x = P.DE, y = TERM)) + geom_point(size=3,color="#E43F3F") +theme_bw() + ylab("") + xlab("-log10(p)") + scale_y_discrete(limits = position1)+xlim(2,3.7)
ggsave("/MOFA_dir/output/GO_factor9_meth_pos_weights.pdf",p7,width=8,height=4)

#### Figure S7 ####
# Integrative GO plot
# Import the different results files for methylation, expression, WGCNA, MOFA, and splicing
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

# GO term analysis 
GO <- compareCluster(geneClusters = BA9_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pvalueCutoff=0.05,ont="BP")
GO<-pairwise_termsim(GO)

p1 <- emapplot(GO,legend_n=3,cex_line=0.1,cex_label_category=0.75,layout="nicely",cex_category=0.8,showCategory=20,pie="equal")+scale_fill_manual(values=c("#A6CEE3", "#1F78B4", "#B2DF8A" ,"#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00" ,"#CAB2D6", "#6A3D9A"))
ggsave(paste0("/results_ranking/GO_BP_BA9_multiome.pdf"),p1,height=10,width=10)

#### Figure S8 ####
# manually created plot