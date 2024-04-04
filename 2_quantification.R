# Author: Eric Zillich, last change: 2024-04-04

library(readr)
library(data.table)
library(dplyr)
library(varhandle)
library(Rsubread)

# Read in samples
CUD_mRNA_files <- list.files("/dir/",pattern="*.Coord.out.bam$")# Quantification mRNA Seq

count <- featureCounts(c(paste0("/dir/",CUD_mRNA_files)), annot.ext="/ref_dir/GRCh38.subset.gtf",isGTFAnnotationFile = T,GTF.featureType = "exon",GTF.attrType = "gene_name",isPairedEnd = T, nthreads = 10,countMultiMappingReads = F,reportReads = "CORE")

j <- "mRNA_CUD"
save(count, file = paste0("/dir/",j,"_counts.Rdata"))






