library(SpatialExperiment)
library(Giotto)
library(plyr)
library(dplyr)

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/plaque_proximity_analysis/"
dir.create(out_dir, recursive = T, showWarnings = F)

instrs <- createGiottoInstructions(
  save_dir = out_dir,
  save_plot = TRUE,
  show_plot = FALSE,
  return_plot = FALSE
)

# read in spatial experiment converted file
m3d_se <- readRDS("../WT_v_APP_PS19/results/anndata_objects/all_cells_integrated.spatial_exp.RDS")
# pull out cell metadata
cell_metadata <- as.data.frame(colData(m3d_se))
cell_metadata$cell_id <- rownames(cell_metadata)

cell_metadata$sample_cond <- paste0(cell_metadata$ID,
                                    "_", cell_metadata$Condition)

# colors for consistent plotting
cell_colors <- metadata(m3d_se)$"Cell Type_colors"

colors_xref <- data.frame(Cell.Type=levels(cell_metadata$Cell.Type),
                          cell_color=cell_colors)

cell_color_meta <- merge(cell_metadata[,c("cell_id","Cell.Type")],
                         colors_xref,
                         by="Cell.Type")
rownames(cell_color_meta) <- cell_color_meta$cell_id

# convert to giotto object ...
m3d_g <- spatialExperimentToGiotto(m3d_se)
# add in z axis, since it does not get added
m3d_g@spatial_locs$cell$raw$sdimz <- spatialCoords(m3d_se)[,"Z"]

# split into samples
sample_ids <- unique(cell_metadata$sample_cond)

sample_order <- c("C164B_WT","C158B_APPPS19",
                  "C166A_WT", "C165_APPPS19")

# separate out samples by subsetting
m3d_samples <- lapply(sample_ids, function(id) {
  
  cell_ids <- rownames(cell_metadata[cell_metadata$sample_cond == id,])
  
  subset(m3d_g, cell_ids=cell_ids)
  
})
names(m3d_samples) <- sample_ids

# read in plaque cengtroids
plaque_cents <- list(C165_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_csv_files/C165_APPPS19_plaque_centers_unfiltered_shifted_um.csv"),
                     C164B_WT=read.csv("../WT_v_APP_PS19/data/plaque_csv_files/C164B_WT_plaque_centers_unfiltered_shifted_um.csv"),
                     C166A_WT=read.csv("../WT_v_APP_PS19/data/plaque_csv_files/C166A_WT_plaque_centers_unfiltered_shifted_um.csv"),
                     C158B_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_csv_files/C158B_APPPS19_plaque_centers_unfiltered_shifted_um.csv"))

# add sample ID and make unique plaque IDs
plaque_cents <- lapply(sample_ids, function(id) {
  
  data <- plaque_cents[[id]]
  data$sample_id <- id
  data$plaque_id <- paste0(data$sample_id,
                           "|", data$plaque_id)
  return(data)
  
})
plaque_cents_df <- bind_rows(plaque_cents)

# filter out likely noise
plaque_cents_df <- plaque_cents_df[plaque_cents_df$total_um > 20 &
                                     plaque_cents_df$z_um < 30 &
                                     plaque_cents_df$total_um < 100849,]

table(plaque_cents_df$sample_id)

# pull in boundary data
plaque_bounds <- list(C165_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_boundary_csvs/well10_c165_appps19_plaque_boundary_pixels.csv"),
                      C158B_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_boundary_csvs/well11_c158b_appps19_plaque_boundary_pixels.csv"))

d_models <- names(plaque_bounds)

# add sample ID and make unique plaque IDs
plaque_bounds <- lapply(names(plaque_bounds), function(id) {
  
  data <- plaque_bounds[[id]]
  data$sample_id <- id
  data$plaque_id <- paste0(data$sample_id,
                           "|", data$Label.ID)
  return(data)
  
})
names(plaque_bounds) <- d_models

plaque_bounds_df <- bind_rows(plaque_bounds)

# filter down boundaries
plaque_bounds_f <- lapply(plaque_bounds, function(data) {
  
  return(data[data$plaque_id %in% plaque_cents_df$plaque_id,])
  
})

lapply(plaque_bounds_f, function(data) {length(unique(data$plaque_id))})

