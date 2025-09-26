library(sceasy)
library(reticulate)
use_condaenv('ad')

ann_file <- "../WT_v_APP_PS19/results/anndata_objects/all_cells_integrated.h5ad"

sceasy::convertFormat(ann_file, from="anndata", to="seurat",
                      outFile='all_cells_integrated.seurat.rds')
