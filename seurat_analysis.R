setwd('/restricted/projectnb/camplab/projects/Single_Cell_Public/Zhang_SignalTransductTargetTher_2022/Analysis/Seurat/')

library(dplyr)
library(Seurat)
library(patchwork)

raw.data <- readRDS('/restricted/projectnb/camplab/projects/Single_Cell_Public/Zhang_SignalTransductTargetTher_2022/Data/raw_data.rds')
# "dgCMatrix" of 17277 x 293432
# from paper (supplementary materials), samples with less than 500 cells were removed. 
# Cells were required to have more than 1000 UMIs and only genes with more than 1000 UMIs across all cells 
# which left 293,432 cells and 17,277 genes

raw <- CreateSeuratObject(counts = raw.data, project = "zhang_signaltransduct")

# percentage of reads that map to the mitochondrial genome 
raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^MT-")

# qc metrics (nCount_RNA, nFeature_RNA, percent.mt)
# nCount_RNA refers to number of RNA molecules detected in cell
# nFeature_RNA refers to number of genes detected in cell
head(raw@meta.data, 5)

# visualize qc metrics
#VlnPlot(raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = FALSE)

# normalize by library size
# default: normalization.method = "LogNormalize", scale.factor = 10000
raw <- NormalizeData(raw)

# subset of features that exhibit high cell-to-cell variation
# default: selection.method = "vst", nfeatures = 2000
raw <- FindVariableFeatures(raw)
top10 <- head(VariableFeatures(raw), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(raw)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
# top 10 variable features:
# ICLC3, IGHG2, IGHM, IGHA2, IGLL5, IGHA1, CCL21, IGLC2, NTS, IGHG4
# Non-variable count: 15277; Variable count: 2000


# scaling after correcting technical covariates (total cellular read count and mitochondrial read count). 
raw <- ScaleData(raw, vars.to.regress = "nCount_RNA, percent.mt")

raw <- RunPCA(raw, features = VariableFeatures(object = raw))

DimPlot(raw, reduction = "pca") + NoLegend()






