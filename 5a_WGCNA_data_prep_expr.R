# Author: Eric Zillich, last change: 2024-04-04

## Code based on WGCNA tutorials by Langfelder & Horvath

## Data Preparation WGCNA 

library(WGCNA)
library(data.table)
library(dplyr)
library(readr)
library(varhandle)
library(matrixStats)
library(DESeq2)

DF <- data.frame

setwd(dir0)
options(stringsAsFactors = FALSE)

#read pheno
pheno <- read.csv("/pheno_dir/pheno.txt", sep=";")
pheno$pH[is.na(pheno$pH)] <- 6.386957

colnames(pheno)[8] <-"AUD"
colnames(pheno)[9] <-"MDD"

# Add cibersort cell type estimates
ciber <- read.csv("/CIBERSORT_dir/Markers_10fold_Yu_He_BA9.txt", sep=";")
ciber <- ciber[,c(1:7)]
colnames(ciber)[1] <- "rn"

pheno_cc <- merge(pheno,ciber,by="rn")

########### EXPRESSION #####


  phe <- pheno_cc[,-2]
  
  datExpr0 <- get(load(file= "/WGCNA_dir/input_BA9_expression.Rdata"))
  datExpr0 <- t(datExpr0)
  
  #check if genes have too many missings
  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK
  
  if(!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }
  
  setwd("/WGCNA_dir/")
  
  sampleTree = hclust(dist(datExpr0), method = "average");

  sizeGrWindow(12,9)
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  datExpr = datExpr0
  
  #traitData
  allTraits = phe
  
  # Form a data frame analogous to expression data that will hold the clinical traits.
  Samples = rownames(datExpr);
  traitRows = match(Samples, allTraits$rn);
  datTraits = allTraits[traitRows, -1];
  rownames(datTraits) = allTraits$rn[traitRows] 
  
    # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(datTraits, signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits),
                      main = "Sample dendrogram and trait heatmap")
  
  save(datExpr, datTraits, file = "BA9_Expr_dataInput_processed.RData")
  
