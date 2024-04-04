# Author: Eric Zillich, last change: 2024-04-04
## WGCNA - associate networks with traits

library(WGCNA)
library(data.table)
library(dplyr)
library(readr)
library(biomaRt)
library(missMethyl)

DF <- data.frame

setwd(dir0)
options(stringsAsFactors = FALSE)

#read pheno
pheno <- read.csv("/pheno_dir/pheno.txt", sep=";")
pheno$pH[is.na(pheno$pH)] <- 6.386957  #mean impute

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#### EXPRESSION ####

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
  write.table(MEs, file = "/WGCNA_dir/MEs_Expr_filtered.txt", sep = ";", quote = F, row.names = T)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  
  pdf("/WGCNA_dir/Heatmap_modules_Expr_filtered.pdf", width = 6, height = 14)
  par(mar = c(6, 8.5, 3, 3));
  labeledHeatmap(Matrix = moduleTraitCor[,c(1:7)],
                 xLabels = names(datTraits)[c(1:7)],
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix[,c(1:7)],
                 setStdMargins = FALSE,
                 cex.text = 0.8,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  
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
  write.table(mod, "/WGCNA_dir/mod_trait_cor_Expr.txt", sep = ";", quote = F)
  
  # GS-MM plot
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
  
  
  # Create the starting data frame
  probes = colnames(datExpr)
  geneInfo0 = data.frame(genx = probes,
                         geneSymbol =probes,
                         moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  
  # Order modules by their significance for weight
  modOrder = order(-abs(cor(MEs, CUD, use = "p")));
  # Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership)){
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  
  # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.CUD));
  geneInfo = geneInfo0[geneOrder, ]
  
  write.table(geneInfo, file = "/WGCNA_dir/geneInfo_Expr.txt", sep = ";", quote = F, row.names = F)
  write.table(table(geneInfo$moduleColor), file = "/WGCNA_dir/geneInfo_Expr_noMods.txt", sep = ";", quote = F, row.names = F)
  write.table(geneInfo0, file = "/WGCNA_dir/geneInfo0_Expr.txt", sep = ";", quote = F, row.names = F)
  
  # For GO enrichment analysis use clusterProfiler
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  
  GO_gene_list <- geneInfo0$geneSymbol[geneInfo0$moduleColor == "yellow"]
  
  # BP 
  
  module_GO_BP <- enrichGO(GO_gene_list,OrgDb = org.Hs.eg.db,keyType = "SYMBOL",pvalueCutoff=1,ont="BP")
  module_GO_BP <-pairwise_termsim(module_GO_BP)
  
  write.table(module_GO_BP@result[module_GO_BP@result$p.adjust<0.05,],"/WGCNA_dir/GO_statistics_yellow.txt", sep = ";", quote = F, row.names = F)
  
  p1 <- emapplot(module_GO_BP,legend_n=3,cex_line=0.1,cex_label_category=0.75,layout="nicely",cex_category=0.8,showCategory=30)
  ggsave("/WGCNA_dir/GO_modules_BP.pdf",p1,height=9,width=9)
  

# Generate input data for the network plots - WGCNA hub genes
  # select significant modules: yellow
  i = "yellow"
    
    geneInfo_Expr0 <- read.csv("/WGCNA_dir/geneInfo0_Expr.txt", sep=";")
    GS <- geneInfo_Expr0[geneInfo_Expr0$moduleColor == i,]
    
    GS$value <- abs(GS[,paste0("MM.",i)]) * abs(GS$GS.CUD)
    GS_gene <- GS[order(GS$value, decreasing = T),]
  
# 10% hub genes 
    GS_gene_10 <- GS_gene[c(1:round(length(GS_gene$genx)/10)),]
    genes_10 <- unique(GS_gene_10$genx)
    write.table(genes_10, "/WGCNA_dir/yl_hub_genes_10pct.txt",quote=F,row.names = F,sep=";")
    

  