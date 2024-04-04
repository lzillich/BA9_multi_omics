# Author: Eric Zillich, last change: 2024-04-04
# Explore the results of drug signature analysis 
library(ggplot2)
library(forcats)

# Individual perturbagens 

query_result <- read.delim("/CMap_BA9/arfs/TAG/query_result.gct")
query_result <-query_result[-1,] # remove descriptive header
#-log10(0.05) -> 1.30103

query_result_sig <-query_result[abs(as.numeric(query_result$fdr_q_nlog10)) > 1.30103,]

pos_corr_perturbagens <- query_result_sig[order(as.numeric(query_result_sig$norm_cs),decreasing = T),][c(1:20),]
neg_corr_perturbagens <- query_result_sig[order(as.numeric(query_result_sig$norm_cs),decreasing = F),][c(1:20),]

# subset for glucocorticoid inside labels
glucocort_perturbagens <- query_result_sig[grep("Glucocor",query_result_sig$moa),]

# GSEA results

gsea_result <- read.delim("/CMap_BA9/gsea/TAG/arfs/NORM_CS/gsea_result.gct")
gsea_result <- gsea_result[-1,]

gsea_result_sig_MOA <- gsea_result[gsea_result$set_type == "MOA_CLASS",]
gsea_result_sig_MOA <- gsea_result_sig_MOA[order(gsea_result_sig_MOA$fdr_q_nlog10,decreasing=T),]
gsea_result_sig_PATH <- gsea_result[gsea_result$set_type == "PATHWAY_SET",]
gsea_result_sig_PATH <- gsea_result_sig_PATH[order(gsea_result_sig_PATH$fdr_q_nlog10,decreasing=T),]
gsea_result_sig_PCL <- gsea_result[gsea_result$set_type == "PCL",]
gsea_result_sig_PCL <- gsea_result_sig_PCL[order(gsea_result_sig_PCL$fdr_q_nlog10,decreasing=T),]

# MOA
pos_corr_gsea_MOA <- gsea_result_sig_MOA[order(as.numeric(gsea_result_sig_MOA$norm_cs),decreasing = T),][c(1:30),]
pos_corr_gsea_MOA$cell_iname[pos_corr_gsea_MOA$cell_iname=="-666"]<-""
pos_corr_gsea_MOA$name <- paste0(pos_corr_gsea_MOA$src_set_id,"_",pos_corr_gsea_MOA$cell_iname)
pos_corr_gsea_MOA$name <- gsub("*_$","",pos_corr_gsea_MOA$name)

neg_corr_gsea_MOA <- gsea_result_sig_MOA[order(as.numeric(gsea_result_sig_MOA$norm_cs),decreasing = F),][c(1:30),]
neg_corr_gsea_MOA$cell_iname[neg_corr_gsea_MOA$cell_iname=="-666"]<-""
neg_corr_gsea_MOA$name <- paste0(neg_corr_gsea_MOA$src_set_id,"_",neg_corr_gsea_MOA$cell_iname)
neg_corr_gsea_MOA$name <- gsub("*_$","",neg_corr_gsea_MOA$name)

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

#PATHWAY
pos_corr_gsea_PATH <- gsea_result_sig_PATH[order(as.numeric(gsea_result_sig_PATH$norm_cs),decreasing = T),][c(1:30),]
neg_corr_gsea_PATH <- gsea_result_sig_PATH[order(as.numeric(gsea_result_sig_PATH$norm_cs),decreasing = F),][c(1:30),]

# waterfall plots, subset top 20 up and down results for visualization

# perturbagens 

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


# glucocort receptor perturbagens 
glucocort_perturbagens <- glucocort_perturbagens[,c("moa","cell_iname","pert_idose","pert_itime","pert_iname","fdr_q_nlog10","norm_cs")]
glucocort_perturbagens$name <- paste0(glucocort_perturbagens$pert_iname,"_",glucocort_perturbagens$cell_iname,"_",glucocort_perturbagens$pert_idose,"_",glucocort_perturbagens$pert_itime)
glucocort_perturbagens$fdr_q_nlog10 <- as.numeric(glucocort_perturbagens$fdr_q_nlog10)
glucocort_perturbagens$norm_cs <- as.numeric(glucocort_perturbagens$norm_cs)
glucocort_perturbagens$fdr_sig <- ifelse(glucocort_perturbagens$fdr_q_nlog10 > 1.30103, "yes", "no")
glucocort_perturbagens <- glucocort_perturbagens[order(glucocort_perturbagens$norm_cs,decreasing = F),]

p1_1 <- ggplot(glucocort_perturbagens, aes(reorder(name,norm_cs), norm_cs)) +
  geom_col(aes(fill=fdr_sig))+ylim(-2,2) +scale_fill_manual(values=c("#1f78bf","gray60")) +
  coord_flip() +
  labs(x="", y="Normalized connectivity score",
       title="Glucocorticoid receptor targeting drugs",legend="sig") + theme_minimal()
ggsave("/CMap_BA9/glucocort_modulators.pdf",p1_1,height=6,width=8)


# gsea PATH
neg_corr_gsea_PATH <- neg_corr_gsea_PATH[,c("src_set_id","fdr_q_nlog10","norm_cs")]
colnames(neg_corr_gsea_PATH) <- c("pathway","pval","cscore_norm")
neg_corr_gsea_PATH$pval <- as.numeric(neg_corr_gsea_PATH$pval)
neg_corr_gsea_PATH$cscore_norm <- as.numeric(neg_corr_gsea_PATH$cscore_norm)

pos_corr_gsea_PATH <- pos_corr_gsea_PATH[,c("src_set_id","fdr_q_nlog10","norm_cs")]
colnames(pos_corr_gsea_PATH) <- c("pathway","pval","cscore_norm")
pos_corr_gsea_PATH$pval <- as.numeric(pos_corr_gsea_PATH$pval)
pos_corr_gsea_PATH$cscore_norm <- as.numeric(pos_corr_gsea_PATH$cscore_norm)

path <- rbind(pos_corr_gsea_PATH,neg_corr_gsea_PATH)
path$fdr_sig <- ifelse(path$pval > 1.30103, "yes", "no")
path <- path[!duplicated(path$pathway),]
path <- path[order(path$cscore_norm,decreasing = T),]
path <- path[c(1:20,40:59),]

p2 <- ggplot(path, aes(fct_reorder(pathway,cscore_norm), cscore_norm)) +
  geom_col(aes(fill=fdr_sig)) +scale_fill_manual(values=c("gray60","#1f78bf")) +
  coord_flip() + theme_minimal() +
  labs(x="", y="Normalized connectivity score",
       title="GSEA pathway",legend="sig")
ggsave("/CMap_BA9/gsea_PATH.pdf",p2,height=6,width=8)

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

p3 <-ggplot(PCL, aes(fct_reorder(PCL,cscore_norm), cscore_norm)) +
  geom_col(aes(fill=fdr_sig)) +
  coord_flip() + scale_fill_manual(values=c("gray60","#1f78bf")) + theme_minimal() +
  labs(x="", y="Normalized connectivity score",
       title="GSEA perturbagen class",legend="sig")
ggsave("/CMap_BA9/gsea_PCL.pdf",p3,height=6,width=8)


# GSEA MOA

neg_corr_gsea_MOA <- neg_corr_gsea_MOA[,c("id","fdr_q_nlog10","norm_cs")]
colnames(neg_corr_gsea_MOA) <- c("MOA","pval","cscore_norm")
neg_corr_gsea_MOA$pval <- as.numeric(neg_corr_gsea_MOA$pval)
neg_corr_gsea_MOA$cscore_norm <- as.numeric(neg_corr_gsea_MOA$cscore_norm)

pos_corr_gsea_MOA <- pos_corr_gsea_MOA[,c("id","fdr_q_nlog10","norm_cs")]
colnames(pos_corr_gsea_MOA) <- c("MOA","pval","cscore_norm")
pos_corr_gsea_MOA$pval <- as.numeric(pos_corr_gsea_MOA$pval)
pos_corr_gsea_MOA$cscore_norm <- as.numeric(pos_corr_gsea_MOA$cscore_norm)

MOA <- rbind(pos_corr_gsea_MOA,neg_corr_gsea_MOA)
MOA$fdr_sig <- ifelse(MOA$pval > 1.30103, "yes", "no")
MOA <- MOA[!duplicated(MOA$MOA),]
MOA <- MOA[order(MOA$cscore_norm,decreasing = T),]
MOA <- MOA[c(1:20,41:60),]

p4 <- ggplot(MOA, aes(fct_reorder(MOA,cscore_norm), cscore_norm)) +
  geom_col(aes(fill=fdr_sig)) +
  coord_flip() + scale_fill_manual(values=c("gray60","#1f78bf")) + theme_minimal()+
  labs(x="", y="Normalized connectivity score",
       title="GSEA mechanism of action",legend="sig")

ggsave("/CMap_BA9/gsea_MOA.pdf",p4,height=6,width=8)

# Create figure panel
library(cowplot)
pdf("/CMap_BA9/CMap_panel.pdf",width=14,height=12)
plot_grid(p1,p1,p3,p1_1,ncol=2,align="v")
dev.off()

