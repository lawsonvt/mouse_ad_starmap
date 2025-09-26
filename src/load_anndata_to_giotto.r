library(Giotto)
library(GiottoData)
library(reticulate)

ann_file <- "../WT_v_APP_PS19/results/anndata_objects/all_cells_pre_integration.h5ad"

set_giotto_python_path("/opt/anaconda3/envs/ad/bin/python")

from_ann <- anndataToGiotto(ann_file)

