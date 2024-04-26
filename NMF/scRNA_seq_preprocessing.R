library(Seurat)
setwd("/home/xbdong/lab/OncoNicheDev")

# loading data, downloaded from 
# https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/
counts <- readRDS("data/Kidney_gene_count.RDS")
genes <- readRDS("data/df_gene.RDS")
cell_anno <- readRDS("data/df_cell.RDS")

# select kidney cells
clean_cell_anno <- cell_anno [(!is.na(cell_anno$Organ)) & (!is.na(cell_anno$Organ_cell_lineage)),  ]
kidney_anno <- clean_cell_anno[clean_cell_anno$Organ == "Kidney", ]
unique(kidney_anno$Organ_cell_lineage)

# only select kidney cells
kidney_cells <- c("Kidney-Metanephric cells")
kidney_cell_anno <- kidney_anno[kidney_anno$Organ_cell_lineage %in% kidney_cells, ]

counts <- counts[, kidney_cell_anno$sample]

# select protein-coding genes
protein_coding_genes <- genes[genes$gene_type=="protein_coding", ]
counts <- counts[protein_coding_genes$gene_id, ]

# name gene as symbol
row.names(counts) <- protein_coding_genes$gene_short_name

kidney <- CreateSeuratObject(counts = counts, project = "kidney_dev", 
                           min.cells = 10, min.features = 200)

# annotate cells
lineage <- kidney_cell_anno$Organ_cell_lineage
days <- kidney_cell_anno$Development_day

names(lineage) <- kidney_cell_anno$sample
names(days) <- kidney_cell_anno$sample

kidney[["cell.lineage"]] <- lineage
kidney[["developmental.days"]] <- days

# Visualize QC metrics as a violin plot
p <- VlnPlot(kidney, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  theme(aspect.ratio=1)
p
ggsave("feature_and_count_QC.pdf",width = 6, height = 4)

p <- FeatureScatter(kidney, feature1 = "nCount_RNA", 
                    feature2 = "nFeature_RNA")+
  theme(aspect.ratio=1)
p
ggsave("count_feature_QC.pdf", width = 6, height = 6)

# select cells, normalize and log(1+x) transform data
kidney <- subset(kidney, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)
kidney <- NormalizeData(kidney, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify the highly variable genes
kidney <- FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 5000)
top10 <- head(VariableFeatures(kidney), 10)


# plot variable features with and without labels
plot1 <- VariableFeaturePlot(kidney)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
p <- plot1 + plot2
p
ggsave("hv_genes.pdf", width = 10, height = 3)

hv <-VariableFeatures(kidney)

cleanData <- kidney[["RNA"]]@data

hv_clean_data <- cleanData[hv, ]

hv_clean_mat <- as.matrix(hv_clean_data)
write.table(hv_clean_mat, file="write/kidney.tsv", sep="\t", quote=FALSE, row.names = TRUE)

# saveRDS(hv_clean_data, file="write/kidney.rds")
