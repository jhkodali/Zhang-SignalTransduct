setwd('/restricted/projectnb/camplab/projects/Single_Cell_Public/Zhang_SignalTransductTargetTher_2022/Analysis/Seurat/')

library(dplyr)
library(Seurat)
library(patchwork)

raw.data <- readRDS('/restricted/projectnb/camplab/projects/Single_Cell_Public/Zhang_SignalTransductTargetTher_2022/Data/raw_data.rds')
# "dgCMatrix" of 17277 x 293432
# from paper (supplementary materials), counts have already gone through cell ranger quality control 
# which left 293,432 cells and 17,277 genes

raw <- CreateSeuratObject(counts = raw.data, project = "zhang")
