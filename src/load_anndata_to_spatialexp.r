# THIS ONE WORKS!!!!

#library(zellkonverter)
library(SpatialExperiment)
library(reticulate)

use_python("/opt/anaconda3/envs/ad/bin/python")

ann_file <- "../WT_v_APP_PS19/results/anndata_objects/all_cells_integrated.h5ad"

ad <- import("anndata", convert=F)

ann_data <- ad$read_h5ad(ann_file)

sce_obj <- zellkonverter::readH5AD(ann_file, reader="R")
uns_r <- py_to_r(ann_data$uns)

sce_obj@metadata <- uns_r

spatial_coords_matrix <- as.matrix(colData(sce_obj)[, c("X", "Y", "Z")])

# create spatial experiment
se_obj <- SpatialExperiment(
  assays=assays(sce_obj),
  colData=colData(sce_obj),
  rowData=rowData(sce_obj),
  metadata=metadata(sce_obj),
  reducedDims=reducedDims(sce_obj),
  spatialCoords = spatial_coords_matrix
)

# save it?
spatial_out <- gsub(".h5ad", ".spatial_exp.RDS", ann_file)

saveRDS(se_obj, spatial_out)

