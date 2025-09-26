library(Seurat)
library(reticulate)

use_python("/opt/anaconda3/envs/ad/bin/python")

ann_file <- "../WT_v_APP_PS19/results/anndata_objects/all_cells_integrated.h5ad"

ad <- import("anndata", convert=F)

ann_data <- ad$read_h5ad(ann_file)

seurat_data <- as.Seurat(ann_data)


