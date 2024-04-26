library(singscore)

library(ggpubr)

library(pheatmap)
library(viridis)

library(umap)

#--------------------------1. Extracting Regulon-----------------------#

# loading search information analysis results
info_change <- read.csv("fitted_observed_info_change.csv")
row.names(info_change) <- info_change$Target

# loading GeneMANIA results
expanded_TF <- read.table("GeneMANIA_expanded.txt", header=TRUE, sep="\t")

# loading TF-target database
# downloaded from https://github.com/lusystemsbio/NetAct
load("hDB.rdata")

# loading TCGA RNA-seq data
tcga <- read.table("kirc-rsem-fpkm-tcga-t.txt", sep="\t", header=TRUE, row.names = 1)
tcga <- tcga[, -1]
tcga <- log1p(tcga)
tcga_genes <- rownames(tcga)

# select significant TFs
#sig_TFs <- info_change[info_change$padj < 0.1, ]$Target
sig_TFs <- expanded_TF$Symbol

# filter targets using TF-target coexpression
cutoff <- 0.03

TFs <- c()
functional_targets <- c()
PCCs <- c()
sign <- c()

for (TF in sig_TFs){
  if (!(TF %in% names(hDB))){ # significant TF should be annotated by hDB
    next
  }
  
  targets <- hDB[TF][[1]]
  
  if (TF %in% tcga_genes){ # TF should have expression data
    print(TF)
    
    TF_expression_vector <- as.numeric(tcga[TF, ])
    
    for(each_t in targets){
      if (each_t %in% tcga_genes){ # target should have expression data
        target_expression_vector <- as.numeric(tcga[each_t, ])
        
        if(sd(target_expression_vector) == 0){
          next
        }
        
        r <- cor(TF_expression_vector, target_expression_vector)
        
        if (abs(r) >= cutoff){ # remove uncorrelated target
          TFs <- c(TFs, TF)
          functional_targets <- c(functional_targets, each_t)
          PCCs <- c(PCCs, r)
          
          if (r > 0){
            sign <- c(sign, "activation")
          }else{
            sign <- c(sign, "inhibition")
          }
        } 
      }
    }
  }
}

regulons <- data.frame(Sig_TF=TFs, Target=functional_targets, PCC=PCCs, Sign=sign)
reg_sig_TFs <- unique(regulons$Sig_TF)

# calculate regulon size
regulon_size <- c()
for (each_TF in reg_sig_TFs){
  size <- nrow(regulons[regulons$Sig_TF==each_TF, ])
  regulon_size <- c(regulon_size, size)
}

regulon_info <- data.frame(Sig_TF=reg_sig_TFs, Regulon_Size=regulon_size)
regulon_info$Delta_S <- info_change[regulon_info$Sig_TF, ]$Delta_S
regulon_info$padj <- info_change[regulon_info$Sig_TF, ]$padj

# remove below comments if you want to save them into files
# write.csv(regulons, file="regulon.csv", row.names = FALSE, quote = FALSE)
# write.csv(regulon_info, file="regulon_info.csv", row.names = FALSE, quote=FALSE)


#--------------------------2. GSEA analysis--------------------------#
get_activated_genes <- function(regulons, TF){
  gs_up <- regulons[(regulons$Sig_TF==TF) & (regulons$Sign == "activation"), ]$Target
  
  return(gs_up)
}

get_inhibited_genes <- function(regulons, TF){
  gs_dn <- regulons[(regulons$Sig_TF==TF) & (regulons$Sign == "inhibition"), ]$Target
  
  return(gs_dn)
}

reg_sig_TFs <- regulon_info[regulon_info$Regulon_Size >= 10, ]$Sig_TF # select stable regulons

# loading TCGA RNA-seq data
tcga_t <- read.table("kirc-rsem-fpkm-tcga-t.txt", sep="\t", header=TRUE) # tumor samples
tcga_n <- read.table("kirc-rsem-fpkm-tcga.txt", sep="\t", header=TRUE) # normal samples
num_of_tumors <- ncol(tcga_t)-2
num_of_normals <- ncol(tcga_n)-2

kirc_tcga_t <- kirc_tcga_t[, -2]
kirc_tcga_n <- kirc_tcga_n[, -2]

tcga <- merge(tcga_t, tcga_n, by="Hugo_Symbol") # combine tumor and normal samples

# preprocessing tcga data
rownames(tcga) <- tcga[, 1] 
tcga <- tcga[, -1]
tcga <- log1p(tcga)

rm(tcga_n, tcga_t)

# loading genotype data of the cancer driver 
gt <- read.csv("TCGA_VHL_MUT.csv", na.strings="")
corrected_barcodes <- c()
for (i in 1:nrow(gt)){
  records <- strsplit(gt$Barcode[i], "-")
  barcode <- paste(records[[1]], collapse = ".")
  corrected_barcodes <- c(corrected_barcodes, barcode)
}
row.names(gt) <- corrected_barcodes
gt$Corrected_Barcode <- corrected_barcodes

mutated_samples <- gt[gt$MUT_State == "Yes", ]$Corrected_Barcode

genotype <- c() # genotype annotation
samples <- colnames(tcga)
for (i in 1:length(samples)){
  if (i <= num_of_tumors){
    records <- strsplit(samples[i], "\\.")
    barcode <- paste(records[[1]][1:3], collapse = ".")
    if (barcode %in% mutated_samples){
      genotype <- c(genotype, "MUT_tumor")
    }else{
      genotype <- c(genotype, "WT_tumor")
    }
  }else{
    genotype <- c(genotype, "normal")
  }
}

# create first two annotation columns
score_df <- data.frame(row.names = colnames(tcga),
                        Sample_Type=c(rep("tumor", num_of_tumors), rep("normal", num_of_normals)),
                        Genotype=genotype)

# single sample GSEA analysis
ranked_tcga <- rankGenes(tcga)
for (i in 1:length(reg_sig_TFs)){
  TF <- reg_sig_TFs[i]
  gs_up <- get_activated_genes(regulons, TF)
  gs_dn <- get_inhibited_genes(regulons, TF)
  
  if ((length(gs_up) > 0) & (length(gs_dn) > 0)){
    score <- simpleScore(ranked_tcga, upSet = gs_up, downSet = gs_dn)
  }else if ((length(gs_up) > 0) & (length(gs_dn) == 0)){
    score <- simpleScore(ranked_tcga, upSet = gs_up)
  }else{
    score <- simpleScore(ranked_tcga, upSet = gs_dn)
  }
  
  score_df <- cbind(score_df, score$TotalScore)
}

colnames(score_df) <- c("Sample_Type", "Genotype", reg_sig_TFs)

# remove below comments if you want to save them into files
# write.csv(score_df, file="gsea_scores.csv", quote = FALSE, row.names = FALSE)

# PCA analysis and visualization
dat <- as.matrix(score_df[, c(3:ncol(score_df))])

pca_res <- prcomp(dat, center = TRUE,scale. = TRUE)
summ <- summary(pca_res)
pc1_lab <- paste0("PC1 (",round(summ$importance[2,1]*100,2),"%)")
pc2_lab <- paste0("PC2 (",round(summ$importance[2,2]*100,2),"%)")

sample_pos <- data.frame(x=pca_res$x[, 1], y=pca_res$x[, 2],
                         Sample_Type=score_df$Sample_Type, Genotype=score_df$Genotype)

ggscatter(sample_pos, x = "x", y = "y", 
          color = "Genotype", palette = c("#00AFBB", "#E7B800", "#FC4E07") )+
  xlab(pc1_lab)+
  ylab(pc2_lab)

# UMAP analysis and visualization
set.seed(42)
dat_umap <- umap(dat)

sample_pos <- data.frame(x=dat_umap$layout[, 1], y=dat_umap$layout[, 2],
                         Sample_Type=score_df$Sample_Type, Genotype=score_df$Genotype)

ggscatter(sample_pos, x = "x", y = "y",
          color = "Genotype", palette = c("#00AFBB", "#E7B800", "#FC4E07") )+
  xlab("UMAP1")+
  ylab("UMAP2")

pheatmap(score_df[, c(3:ncol(score_df))], scale="column", show_rownames=FALSE, border_color = NA,
         breaks=seq(from=-3, to=3, by=0.06),
         color=viridis(100))

# boxplot analysis and visualization
my_comparisons <- list( c("MUT_tumor", "WT_tumor"), c("MUT_tumor", "normal"), c("WT_tumor", "normal") )
ggboxplot(score_df, x = "Genotype", y = "EPAS1", width = 0.5, order=c("normal", "WT_tumor", "MUT_tumor"), 
          fill = "Genotype", palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(0.4, 0.5, 0.3))+
  xlab(NULL)+guides(fill = "none")

