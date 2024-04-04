# Author: Eric Zillich, last change: 2024-04-04
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

# Display the current working directory
setwd(dir0)

options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads = 20)


###### EXPRESSION ########

        setwd("WGCNA_dir/")
        lnames = load(file ="BA9_Expr_dataInput_processed.RData")    
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        
        sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
        
        sft$powerEstimate
        
        net = blockwiseModules(datExpr, power = sft$powerEstimate,
                               TOMType = "signed", minModuleSize = 10, maxBlockSize = 36000,
                               reassignThreshold = 0, mergeCutHeight = 0.25,
                               numericLabels = TRUE, pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = "/dir2/BA9_mRNA_TOM",
                               verbose = 3)
        table(net$colors)
        
        setwd("/dir2/")
        
        # Convert labels to colors for plotting
        mergedColors = labels2colors(net$colors)
        # Plot the dendrogram and the module colors underneath
        png("BA9_dendogram_Expr.png", width = 1200, height = 1800)
        plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                            "Module colors",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05)
        dev.off()
        
        moduleLabels = net$colors
        moduleColors = labels2colors(net$colors)
        MEs = net$MEs;
        geneTree = net$dendrograms[[1]];
        save(MEs, moduleLabels, moduleColors, geneTree,
             file = "BA9_networkConstruction_Expr.RData")
        
