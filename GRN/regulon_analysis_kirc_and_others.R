library(singscore)

library(ggpubr)
library(cowplot)

library(ggplot2)

library(umap)

library(factoextra)

library(pheatmap)
library(RColorBrewer)
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

# loading KIRC TCGA RNA-seq data
kirc_tcga_t <- read.table("kirc-rsem-fpkm-tcga-t.txt", sep="\t", header=TRUE) # tumor samples
kirc_tcga_n <- read.table("kirc-rsem-fpkm-tcga.txt", sep="\t", header=TRUE) # normal samples
kirc_num_of_tumors <- ncol(kirc_tcga_t)-2
kirc_num_of_normals <- ncol(kirc_tcga_n)-2

kirc_tcga_t <- kirc_tcga_t[, -2]
kirc_tcga_n <- kirc_tcga_n[, -2]

tcga <- merge(kirc_tcga_t, kirc_tcga_n, by="Hugo_Symbol") # combine tumor and normal samples

# loading BRCA TCGA RNA-seq data
brca_tcga_t <- read.table("brca-rsem-fpkm-tcga-t.txt", sep="\t", header=TRUE) # tumor samples
brca_num_of_tumors <- ncol(brca_tcga_t)-2

brca_tcga_t <- brca_tcga_t[, -2]

tcga <- merge(tcga, brca_tcga_t, by="Hugo_Symbol") # add BRCA samples

# loading COAD TCGA RNA-seq data
coad_tcga_t <- read.table("coad-rsem-fpkm-tcga-t.txt", sep="\t", header=TRUE) # tumor samples
coad_num_of_tumors <- ncol(coad_tcga_t)-2

coad_tcga_t <- coad_tcga_t[, -2]

tcga <- merge(tcga, coad_tcga_t, by="Hugo_Symbol") # add COAD samples

# preprocessing tcga data
rownames(tcga) <- tcga[, 1] 
tcga <- tcga[, -1]
tcga <- log1p(tcga)

rm(kirc_tcga_t, kirc_tcga_n, brca_tcga_t, coad_tcga_t)

# create first two annotation columns
score_df <- data.frame(row.names = colnames(tcga),
                        Sample_Type=c(rep("KIRC_tumor", kirc_num_of_tumors), 
                                      rep("KIRC_normal", kirc_num_of_normals),
                                      rep("BRCA_tumor", brca_num_of_tumors),
                                      rep("COAD_tumor", coad_num_of_tumors))
                       )

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

colnames(score_df) <- c("Sample_Type", reg_sig_TFs) # naming columns

# remove below comments if you want to save them into files

# write.csv(score_df, file="gsea_scores.csv", quote = FALSE, row.names = TRUE)

# PCA analysis and visualization
dat <- as.matrix(score_df[, c(2:ncol(score_df))])

pca_res <- prcomp(dat, center = TRUE, scale. = TRUE)
summ <- summary(pca_res)
pc1_lab <- paste0("PC1 (",round(summ$importance[2,1]*100,2),"%)")
pc2_lab <- paste0("PC2 (",round(summ$importance[2,2]*100,2),"%)")

sample_pos <- data.frame(x=pca_res$x[, 1], y=pca_res$x[, 2],
                         Sample_Type=score_df$Sample_Type)

p1 <- ggscatter(sample_pos, x = "x", y = "y", color = "Sample_Type",
          size = 1, alpha=0.7,
          palette = c("#e41a1c", "#377eb8", "#4daf4a","#984ea3"))+
  xlab(pc1_lab)+
  ylab(pc2_lab)+coord_fixed(1)+theme_cowplot(8)

fviz_eig(pca_res)

p2 <- fviz_pca_var(pca_res,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), col.circle="grey", title=,
             repel = TRUE     # Avoid text overlapping
)+xlab(pc1_lab)+ylab(pc2_lab)+theme_bw()+theme_cowplot(8)

ggarrange(p1, p2, ncol = 2, nrow = 1, align="hv")
ggsave("regulon_PCA.pdf", width=174, height=80, units="mm")


contrib1 <- fviz_contrib(pca_res, choice = "var", axes = 1, top = Inf)+theme_pubr()

contrib2 <- fviz_contrib(pca_res, choice = "var", axes = 2, top = Inf)+theme_pubr()

ggarrange(contrib1, contrib2, ncol = 1, nrow = 2)

# ggsave("regulon_PCA.pdf", width=85, height=85, units="mm")

# UMAP analysis and visualization
set.seed(42)
dat_umap <- umap(dat)

sample_pos <- data.frame(x=dat_umap$layout[, 1], y=dat_umap$layout[, 2],
                         Sample_Type=score_df$Sample_Type)

ggscatter(sample_pos, x = "x", y = "y",
          color = "Sample_Type", palette = c("#00AFBB", "#E7B800", "#FC4E07","#66c2a5") )+
  xlab("UMAP1")+
  ylab("UMAP2")


# correlation matrix

dat <- as.matrix(score_df[score_df$Sample_Type =="KIRC_tumor", c(2:ncol(score_df))])


dat <- scale(dat)

M <-cor(dat)
color <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

pdf(file="regulator_cor_KIRC.pdf", width=5.5, height=5)
pheatmap(M, color=color, breaks=seq(from=-1, to=1, by=0.02),
         border_color="black",
         clustering_method="ward.D2",
         cutree_rows=3, cutree_cols=3)

dev.off()



dat <- as.matrix(score_df[, c(2:ncol(score_df))])


dat <- scale(dat)

M <-cor(dat)
color <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

pdf(file="regulator_cor_KIRC_noraml_BRCA_COAD.pdf", width=5.5, height=5)
pheatmap(M, color=color, breaks=seq(from=-1, to=1, by=0.02),
         border_color="black",
         clustering_method="ward.D2",
         cutree_rows=2, cutree_cols=2)

dev.off()