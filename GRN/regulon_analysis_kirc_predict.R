library(singscore)

get_activated_genes <- function(regulons, TF){
  gs_up <- regulons[(regulons$Sig_TF==TF) & (regulons$Sign == "activation"), ]$Target
  
  return(gs_up)
}

get_inhibited_genes <- function(regulons, TF){
  gs_dn <- regulons[(regulons$Sig_TF==TF) & (regulons$Sign == "inhibition"), ]$Target
  
  return(gs_dn)
}

gsea_scoring <- function(filename, tcga_genes, regulon, regulon_info){
  test <- read.csv(filename)
  
  comm_genes <- intersect(tcga_genes$Hugo_Symbol, test$symbol)
  
  test <-test[test$symbol %in% comm_genes, ]
  
  dat <- data.frame(Test=test$fpkm)
  row.names(dat) <- test$symbol
  
  dat <- log1p(dat)
  
  # create annotation column
  score_df <- data.frame(row.names = colnames(dat),
                         Sample_Type=rep("Test", ncol(dat))
  )
  
  # single sample GSEA analysis
  ranked_dat <- rankGenes(dat)
  for (i in 1:length(reg_sig_TFs)){
    TF <- reg_sig_TFs[i]
    gs_up <- get_activated_genes(regulons, TF)
    gs_dn <- get_inhibited_genes(regulons, TF)
    
    if ((length(gs_up) > 0) & (length(gs_dn) > 0)){
      score <- simpleScore(ranked_dat, upSet = gs_up, downSet = gs_dn)
    }else if ((length(gs_up) > 0) & (length(gs_dn) == 0)){
      score <- simpleScore(ranked_dat, upSet = gs_up)
    }else{
      score <- simpleScore(ranked_dat, upSet = gs_dn)
    }
    
    score_df <- cbind(score_df, score$TotalScore)
  }
  
  colnames(score_df) <- c("Sample_Type", reg_sig_TFs) # naming columns
  
  return(score_df)
}

gsea_scoring_tokyo <- function(filename, tcga_genes, regulon, regulon_info){
  test <- read.csv(filename)
  
  comm_genes <- intersect(tcga_genes$Hugo_Symbol, test$Symbol)
  
  test <-test[test$Symbol %in% comm_genes, ]
  
  dat <- data.frame(test[, 2:ncol(test)])
  row.names(dat) <- test$Symbol
  
  # create annotation column
  score_df <- data.frame(row.names = colnames(dat),
                         Sample_Type=rep("Test", ncol(dat))
  )
  
  # single sample GSEA analysis
  ranked_dat <- rankGenes(dat)
  for (i in 1:length(reg_sig_TFs)){
    TF <- reg_sig_TFs[i]
    gs_up <- get_activated_genes(regulons, TF)
    gs_dn <- get_inhibited_genes(regulons, TF)
    
    if ((length(gs_up) > 0) & (length(gs_dn) > 0)){
      score <- simpleScore(ranked_dat, upSet = gs_up, downSet = gs_dn)
    }else if ((length(gs_up) > 0) & (length(gs_dn) == 0)){
      score <- simpleScore(ranked_dat, upSet = gs_up)
    }else{
      score <- simpleScore(ranked_dat, upSet = gs_dn)
    }
    
    score_df <- cbind(score_df, score$TotalScore)
  }
  
  colnames(score_df) <- c("Sample_Type", reg_sig_TFs) # naming columns
  
  return(score_df)
}

gsea_scoring_batch <- function(filename, tcga_genes, regulon, regulon_info){
  test <- read.csv(filename, row.names = 1)
  
  comm_genes <- intersect(tcga_genes$Hugo_Symbol, test$Symbol)
  
  test <-test[test$Symbol %in% comm_genes, ]
  
  dat <- data.frame(test[, 3:ncol(test)])
  row.names(dat) <- test$Symbol
  
  dat <- log1p(dat)
  
  # create annotation column
  score_df <- data.frame(row.names = colnames(dat),
                         Sample_Type=rep("Test", ncol(dat))
  )
  
  # single sample GSEA analysis
  ranked_dat <- rankGenes(dat)
  for (i in 1:length(reg_sig_TFs)){
    TF <- reg_sig_TFs[i]
    gs_up <- get_activated_genes(regulons, TF)
    gs_dn <- get_inhibited_genes(regulons, TF)
    
    if ((length(gs_up) > 0) & (length(gs_dn) > 0)){
      score <- simpleScore(ranked_dat, upSet = gs_up, downSet = gs_dn)
    }else if ((length(gs_up) > 0) & (length(gs_dn) == 0)){
      score <- simpleScore(ranked_dat, upSet = gs_up)
    }else{
      score <- simpleScore(ranked_dat, upSet = gs_dn)
    }
    
    score_df <- cbind(score_df, score$TotalScore)
  }
  
  colnames(score_df) <- c("Sample_Type", reg_sig_TFs) # naming columns
  
  return(score_df)
}


#----------------------1. loading regulons------    --------------------#
regulons <- read.csv("regulon.csv")
regulon_info <- read.csv("regulon_info.csv")

reg_sig_TFs <- regulon_info[regulon_info$Regulon_Size >= 10, ]$Sig_TF # select stable regulons

#----------------------2. loading TCGA-KIRC gene lists------------------#
# kirc_tcga_t <- read.table("kirc-rsem-fpkm-tcga-t.txt", sep="\t", header=TRUE) # tumor samples
# kirc_tcga_genes <- data.frame(Hugo_Symbol=kirc_tcga_t$Hugo_Symbol,
#                               Entrez=kirc_tcga_t$Entrez_Gene_Id)
# 
# write.csv(kirc_tcga_genes, file = "kirc_tcga_gene_list.csv", quote = FALSE)

kirc_tcga_genes <- read.csv("kirc_tcga_gene_list.csv", row.names = 1)

#----------------------3. loading test RNA-seq data and scoring-- ------#
# loading other RNA-seq data
files <- file.path("cell_lines", c("VMRC-RCW_rnaseq.csv", 
                                   "KMRC-20_rnaseq.csv",
                                   "786-0_rnaseq.csv",
                                   "SK-NEP-1.csv"))

# predict VMRC-RCW
score_df <- gsea_scoring(files[1], kirc_tcga_genes, regulon, regulon_info)
row.names(score_df) <- c("VMRC-RCW")

output_file <- file.path("cell_lines", "VMRC_RCW_gsea_scores.csv")
write.csv(score_df, file=output_file, quote = FALSE, row.names = TRUE)


# predict KMRC-20
score_df <- gsea_scoring(files[2], kirc_tcga_genes, regulon, regulon_info)
row.names(score_df) <- c("KMRC-20")

output_file <- file.path("cell_lines", "KMRC_20_gsea_scores.csv")
write.csv(score_df, file=output_file, quote = FALSE, row.names = TRUE)

# predict 786-0
score_df <- gsea_scoring(files[3], kirc_tcga_genes, regulon, regulon_info)
row.names(score_df) <- c("786-0")

output_file <- file.path("cell_lines", "786-0_gsea_scores.csv")
write.csv(score_df, file=output_file, quote = FALSE, row.names = TRUE)


# predict SK-NEP-1
score_df <- gsea_scoring(files[4], kirc_tcga_genes, regulon, regulon_info)
row.names(score_df) <- c("SK-NEP-1")

output_file <- file.path("cell_lines", "SK-NEP-1_gsea_scores.csv")
write.csv(score_df, file=output_file, quote = FALSE, row.names = TRUE)


# predict multiple cell lines
filename = file.path("cell_lines", "RCC_celline_rna_seq_tpm.csv")

score_df <- gsea_scoring_batch(filename, kirc_tcga_genes, regulon, regulon_info)

cell_line_rna_seq <- read.csv(filename, row.names = 1)
row.names(score_df) <- colnames(cell_line_rna_seq)[3:ncol(cell_line_rna_seq)]

output_file <- file.path("cell_lines", "RCC_gsea_scores.csv")
write.csv(score_df, file=output_file, quote = FALSE, row.names = TRUE)


# predict Tokyo cohort
filename = file.path("ccRCC_Tokyo", "ccRCC_Tokyo.csv")

score_df <- gsea_scoring_tokyo(filename, kirc_tcga_genes, regulon, regulon_info)

tokyo_rna_seq <- read.csv(filename)
row.names(score_df) <- colnames(tokyo_rna_seq)[2:ncol(tokyo_rna_seq)]

output_file <- file.path("ccRCC_Tokyo", "ccRCC_gsea_scores.csv")
write.csv(score_df, file=output_file, quote = FALSE, row.names = TRUE)
