# Author: Eric Zillich, last change: 2024-04-04
# Overlap between methylation and expression 

library(stringr)
DE_results <- read.csv("/results_dir/DE_results.txt")
DE_results <- DE_results$Gene[DE_results$pvalue <0.05]
EWAS_results <- read.csv("/results_dir/CUD_EWAS_results_lm.txt", sep="")
EWAS_results <- EWAS_results$Gene[EWAS_results$P_VAL <0.01]
EWAS_results <- unique(unlist(str_split(EWAS_results,pattern=";")))[-1]

diffmethdiffexp <- intersect(EWAS_results, DE_results)

# Create correlation plot for DE with DiffMeth for the involved genes 
DE_results2 <- read.csv("/results_dir/DE_results.txt")
DE_results2 <- DE_results2[DE_results2$pvalue <0.05 & DE_results2$Gene %in% diffmethdiffexp,]

EWAS_results <- read.csv("/results_dir/CUD_EWAS_results_lm.txt", sep="")
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


# Create results dataframe with the genes for CpGs 
mat <- matrix(nrow=0,ncol=33,NA)
colnames(mat) <- colnames(EWAS_results2)
for(i in diffmethdiffexp){
  res <- EWAS_results2[EWAS_results2$gene1 == i  | EWAS_results2$gene2 == i | EWAS_results2$gene3 == i | EWAS_results2$gene4 == i |  EWAS_results2$gene5 == i| EWAS_results2$gene6 == i |  EWAS_results2$gene7 == i|  EWAS_results2$gene8 == i|
                         EWAS_results2$gene9 == i|EWAS_results2$gene10 == i  |EWAS_results2$gene11 == i  |EWAS_results2$gene12 == i  |EWAS_results2$gene13 == i  |EWAS_results2$gene14 == i  |EWAS_results2$gene15 == i  | EWAS_results2$gene16 == i |EWAS_results2$gene17 == i  |EWAS_results2$gene18 == i  |EWAS_results2$gene19 == i  |EWAS_results2$gene20 == i| EWAS_results2$gene21 == i |EWAS_results2$gene22 == i, ]
  nrow(res)
  mat <- rbind(mat,res)
}

# non-unique genes where multiple CpG sites are available :keep CpG site with strongest significance 

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

write.table(methexpr, file="/results_dir/methexpr_corr.txt", row.names=F,col.names = T, quote=F, sep=";")

library(ggplot2)
customPlot <- list(
  theme_minimal(base_size = 12), 
  scale_fill_manual(values=c("#E43F3F","#268989","gray20","gray60")), 
  scale_colour_manual(values=c("#E43F3F","#268989","gray20","gray60"))
)

library(ggrepel)
pdf("/results_dir/DE_DiffMeth_pabel.pdf", width = 9, height = 6)
ggplot(methexpr , aes(x=BETA, y=log2FoldChange,colour=corr)) + geom_point(size=1) +xlim(-5.5,5.5)+ylim(-5.5,5.5)+ xlab("Differential Methylation (effect size)")+ylab("Differential Expression (log2FC)")+theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+geom_text_repel(aes(label = gene_final), vjust = 0, nudge_y = 0.1, size = 3, max.overlaps=10)+customPlot+geom_hline(aes(yintercept = 0),size=0.1) + geom_vline(aes(xintercept = 0),size=0.1) +theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
dev.off()

#table(methexpr$corr)
# methylation down- expression up   methylation up - expression down methylation down - expression down     methylation up - expression up 
#29                                 63                                 58                                 33 

#### Overlap methylation and expression with Splicing ####
load("/leafviz_dir/leafviz.RData")
clusters <- merge(introns,clusters, by="clusterID")
clusters <- clusters[clusters$FDR <0.05 & abs(clusters$deltapsi) > 0.025,]
clusters$gene.x <- gsub("<i>","",clusters$gene.x)
clusters$gene.x <- gsub("</i>","",clusters$gene.x)
AS_genes <- unique(clusters$gene.x[clusters$gene.x != "." & clusters$FDR < 0.05])
intersect(intersect(EWAS_results, DE_results),AS_genes) # "ZBTB4"  "INPP5E"

library(VennDiagram)
venn.diagram(list(DE = DE_results, AS = AS_genes, Meth = EWAS_results),filename="/results_dir/DE_genes_DiffMeth_genes_AS_genes_overlap.png", imagetype = "png",fill=c("#1f78bf","#b2df8a","#a6cee3"),width=4.5,height=4.5,units="in",resolution=600,category.names = c("","",""))
