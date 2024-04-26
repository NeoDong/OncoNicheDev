library(TCGAbiolinks)
library(survminer)
library(survival)

library(pheatmap)
library(viridis)

library(maftools)

library(reshape2)
library(ggpubr)

library(cowplot)

# load clustering results
clusters <- read.csv("cluster_res.csv")

# add barcode
barcodes <- c()

for (sample in clusters$Sample){
  barcode <- paste0(strsplit(sample, "\\.")[[1]][1:3], collapse="-")
  barcodes <- c(barcodes, barcode)
}

clusters$Barcode <- barcodes

rm(barcode, barcodes, sample)

#----------------------- ----1. Survival analysis------------------------------#
# function to extract survival data from TCGA clinical data
get_TCGA_survival_data<- function(data, clusterCol = NULL){
    if (!all(c("vital_status", "days_to_death", "days_to_last_follow_up") %in% 
             colnames(data))) 
      stop("Columns vital_status, days_to_death and  days_to_last_follow_up should be in data frame")
  
    group <- NULL
    if (is.null(clusterCol)) {
      stop("Please provide the clusterCol argument")
    }
    else if (length(unique(data[, clusterCol])) == 1) {
      stop(paste0("Sorry, but I'm expecting at least two groups\n", 
                  "  Only this group found: ", unique(data[, clusterCol])))
    }
    notDead <- is.na(data$days_to_death)
    if (any(notDead == TRUE)) {
      data[notDead, "days_to_death"] <- data[notDead, "days_to_last_follow_up"]
    }
    if (length(data[which((data[, "days_to_death"] < 0) == T), 
                    "sample"]) > 0 & "sample" %in% colnames(data)) {
      message("Incosistencies in the data were found. The following samples have a negative days_to_death value:")
      message(paste(data[which((data[, "days_to_death"] < 0) == 
                                 T), "sample"], collapse = ", "))
    }
    if (any(is.na(data[, "days_to_death"])) & "sample" %in% colnames(data)) {
      message("Incosistencies in the data were found. The following samples have a NA days_to_death value:")
      message(paste(data[is.na(data[, "days_to_death"]), "sample"], 
                    collapse = ", "))
    }
    data$status <- grepl("dead|deceased", data$vital_status, ignore.case = TRUE)
    data$group <- as.factor(data[, clusterCol])
    data <- data[, c("days_to_death", "status", "group")]
    
    return(data)
}

# downloading TCGA clinical data
clin.gbm <- GDCquery_clinic("TCGA-KIRC", "clinical")
write.csv(clin.gbm, "TCGA_clinical_data.csv", quote=TRUE)

# extract survival data
clin.core.data <- data.frame(Barcode=clin.gbm$submitter_id, 
                             vital_status=clin.gbm$vital_status, 
                             days_to_death=clin.gbm$days_to_death,
                             days_to_last_follow_up=clin.gbm$days_to_last_follow_up,
                             sex=clin.gbm$gender,
                             stage=clin.gbm$ajcc_pathologic_stage,
                             age=clin.gbm$age_at_index)

row.names(clin.core.data) <- clin.core.data$Barcode

# combine clustering results and clinical data
clusters.clin <- merge(clusters, clin.core.data, by="Barcode")

# plot survival curve and do logrank test
surv_df <- get_TCGA_survival_data(clusters.clin, clusterCol="Cluster_Idx")
surv_df$time <- surv_df$days_to_death/30

fit <- survfit(Surv(time, status) ~ group, data = surv_df)

# run this to get exact p-value
# survdiff(Surv(time, status) ~ group, data = surv_df)

g1 <- ggsurvplot(fit, data = surv_df, pval = 2e-08, pval.method=TRUE,
           xlab="Months", ylab="Probability of overall survival",
           palette = c("#CD534CFF", "#0073C2FF","#EFC000FF"))

surv_df$sex <- clusters.clin$sex
surv_df$stage <- clusters.clin$stage
surv_df$age <- clusters.clin$age

model <- coxph(Surv(time, status) ~ group + sex + age + stage, data = surv_df )

g2 <- ggforest(model, main = "Hazard ratio of KIRC", fontsize = 0.9)+
  theme_cowplot()

plot_grid(g1$plot, g2, rel_widths = c(1, 1.3))
ggsave("combined_survival_plot.pdf")

#-----------------------2.Plotting heatmap of clusters-------------------------#

# loading GSEA data and preprocessing
gsea <- read.csv("gsea_scores.csv", row.names = 1)
gsea <- gsea[gsea$Sample_Type=="KIRC_tumor", -1]

gsea_scaled <- scale(gsea)

gsea_scaled_ordered <- gsea_scaled[clusters$Sample, ]

gsea_scaled_ordered <- t(gsea_scaled_ordered)

# loading VHL driver mutation
vhl_mut <- read.csv("TCGA_VHL_MUT.csv")
row.names(vhl_mut) <- vhl_mut$Barcode

# downloading cancer subtypes provided by TCGA
subtypes <- PanCancerAtlas_subtypes()
subtypes <- subtypes[subtypes$cancer.type == "KIRC", ]
# write.csv(subtypes, file="TCGA_KIRC_subtypes.csv", quote=FALSE)

kirc_subtypes <- data.frame(Barcode=subtypes$pan.samplesID, 
                            TCGA_Subtype_mRNA=paste0("m", subtypes$Subtype_mRNA, sep=""))

row.names(kirc_subtypes) <- kirc_subtypes$Barcode

# create cluster annotation
cluster_anno <- data.frame(VHL_MUT=vhl_mut[clusters$Barcode, "MUT_State"],
                           TCGA_Subtype =kirc_subtypes[clusters$Barcode, "TCGA_Subtype_mRNA"],
                           Cluster=paste0('C', clusters$Cluster_Idx, sep=""))

row.names(cluster_anno) <- clusters$Sample

# assigning colors
ann_colors=list(Cluster=c(C0="#CD534CFF",C1="#0073C2FF",C2="#EFC000FF"),
                VHL_MUT=c(Yes="red", No="grey"),
                TCGA_Subtype=c(m1="#E64B35FF", m2="#4DBBD5FF", m3="#00A087FF", m4="#3C5488FF", mNA="white"))

color <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
color <- inferno(100)

group0_num <- sum(clusters$Cluster_Idx == 0)
group1_num <- sum(clusters$Cluster_Idx == 1)
group2_num <- sum(clusters$Cluster_Idx == 2)

pheatmap(gsea_scaled_ordered, cluster_rows=TRUE, cluster_cols=FALSE,
         show_colnames=FALSE, color=color, breaks=seq(from=-3, to=3, by=0.06),
         cutree_rows=3, gaps_col=c(group0_num, group0_num+group1_num),
         annotation_col=cluster_anno, annotation_colors = ann_colors)

#-----------------------3. Mutated gene analysis-------------------------#
query <- GDCquery(
  project = "TCGA-KIRC", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

GDCdownload(query)

maf <- GDCprepare(query)
write.csv(maf, file="TCGA_KIRC_SNV.maf.csv", quote=TRUE)

maf <- read.csv("TCGA_KIRC_SNV.maf.csv", row.name=1)

kirc_genes <- c("VHL", "PBRM1", "SETD2", "KDM5C", "PTEN", "BAP1",
                "MTOR", "TP53")

kirc_maf <- maf[, ]
kirc_maf <- kirc_maf %>% read.maf

# plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', 
#               dashboard = FALSE)

oncoplot(maf = kirc_maf, genes = kirc_genes)

# add barcode
barcodes <- c()

for (sample in maf$Tumor_Sample_Barcode){
  barcode <- paste0(strsplit(sample, "-")[[1]][1:3], collapse="-")
  barcodes <- c(barcodes, barcode)
}

maf$Barcode <- barcodes
rm(barcode, barcodes, sample)


C0_samples <- clusters[clusters$Cluster_Idx == 0, ]$Barcode
C1_samples <- clusters[clusters$Cluster_Idx == 1, ]$Barcode
C2_samples <- clusters[clusters$Cluster_Idx == 2, ]$Barcode

C0_maf <- maf[maf$Barcode %in% C0_samples, ]
C1_maf <- maf[maf$Barcode %in% C1_samples, ]
C2_maf <- maf[maf$Barcode %in% C2_samples, ]

C0_maf <- C0_maf %>% read.maf
C1_maf <- C1_maf %>% read.maf
C2_maf <- C2_maf %>% read.maf

oncoplot(maf = C0_maf, genes = kirc_genes)

oncoplot(maf = C1_maf, genes = kirc_genes, legend_height=4.5)

oncoplot(maf = C2_maf, genes = kirc_genes)

C0_gene_mut <- as.data.frame(getGeneSummary(C0_maf))
C1_gene_mut <- as.data.frame(getGeneSummary(C1_maf))
C2_gene_mut <- as.data.frame(getGeneSummary(C2_maf))

row.names(C0_gene_mut) <- C0_gene_mut$Hugo_Symbol
row.names(C1_gene_mut) <- C1_gene_mut$Hugo_Symbol
row.names(C2_gene_mut) <- C2_gene_mut$Hugo_Symbol

C0_num <- dim(getSampleSummary(C0_maf))[1]
C1_num <- dim(getSampleSummary(C1_maf))[1]
C2_num <- dim(getSampleSummary(C2_maf))[1]

mut_freq <- data.frame(Hugo_Symbol=kirc_genes, 
                       C0=C0_gene_mut[kirc_genes, ]$AlteredSamples/C0_num, 
                       C1=C2_gene_mut[kirc_genes, ]$AlteredSamples/C1_num,
                       C2=C2_gene_mut[kirc_genes, ]$AlteredSamples/C2_num)

mut_freq[is.na(mut_freq)] <- 0

mut_freq <- melt(mut_freq)
colnames(mut_freq) <-c("Gene", "Cluster_Idx", "Frequency")

ggbarplot(mut_freq, x = "Gene", y = "Frequency", 
          fill = "Cluster_Idx", palette = c("#CD534CFF", "#0073C2FF","#EFC000FF"),
          position = position_dodge(0.8))+scale_y_continuous(expand = c(0,0), limits = c(0,0.6))


#-----------------------4. Survival analysis of Tokyo cohort------------------#
# load clustering results
clusters <- read.csv(file.path("ccRCC_Tokyo", "ccRCC_pred.csv"), row.names = 1)
clusters$Sample <- row.names(clusters)

# add barcode
barcodes <- c()

for (sample in clusters$Sample){
  barcode <- paste0(strsplit(sample, "\\.")[[1]][1:2], collapse="-")
  barcodes <- c(barcodes, barcode)
}

clusters$Barcode <- barcodes

rm(barcode, barcodes, sample)

clinical <- read.csv(file.path("ccRCC_Tokyo", "clinicalData.csv"))
row.names(clinical) <- clinical$sample.ID

# combine clustering results and clinical data
selected_clinical <- clinical[clusters$Barcode, ]

# plot survival curve and do logrank test

surv_df <- data.frame(time=selected_clinical$observation.period..month., 
                      group=as.factor(clusters$Predicted_Class))

surv_df$status <- grepl("dead", selected_clinical$outcome, ignore.case = TRUE)


fit <- survfit(Surv(time, status) ~ group, data = surv_df)

# run this to get exact p-value
# survdiff(Surv(time, status) ~ group, data = surv_df)

ggsurvplot(fit, data = surv_df, pval = 2e-08, pval.method=TRUE,
           xlab="Months", ylab="Probability of overall survival",
           palette = c("#CD534CFF", "#0073C2FF","#EFC000FF"))

sum(selected_clinical$outcome == "dead")
