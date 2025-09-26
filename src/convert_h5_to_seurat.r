library(Seurat)
library(SeuratDisk)
library(reticulate)

use_python("/opt/anaconda3/envs/ad/bin/python")

ann_file <- "../WT_v_APP_PS19/results/anndata_objects/all_cells_integrated.h5ad"

Convert(ann_file, dest="h5seurat")

h5s_file <- gsub(".h5ad", ".h5seurat", ann_file)

data <- LoadH5Seurat(h5s_file)
