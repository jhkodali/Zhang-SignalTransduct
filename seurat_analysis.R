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
head(raw@meta.data, 5)

# visualize qc metrics
VlnPlot(raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = FALSE)

# normalize by library size
