library(SpatialExperiment)
library(Giotto)
library(plyr)
library(dplyr)
library(plotly)

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/plaques_as_sphers_3D_plots/"
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

# add sample ID, make unique plaque IDs, and filter
plaque_cents <- lapply(sample_ids, function(id) {
  
  data <- plaque_cents[[id]]
  data$sample_id <- id
  data$plaque_id <- paste0(data$sample_id,
                           "|", data$plaque_id)
  
  data <- data[data$total_um > 20 &
                 data$z_um < 30 &
                 data$total_um < 100849,]
  
  data$radius <- (3 * data$total_um / (4 * pi))^(1/3)
  
  return(data)
  
})
names(plaque_cents) <- sample_ids
plaque_cents_df <- bind_rows(plaque_cents)

# make a plot for each 

for (sample_id in sample_ids) {
  
  plaque_data <- plaque_cents[[sample_id]]
  plaque_data <- plaque_data[,c("x_um","y_um","z_um","radius")]
  colnames(plaque_data) <- c("x","y","z","radius")
  
  cell_data <- cell_metadata[cell_metadata$sample_cond == sample_id,
                             c("X","Y", "Z", "Cell.Type")]
  colnames(cell_data) <- c("x","y","z","cell_type")
  
  p <- plot_ly(
    data = cell_data,
    x = ~x, y = ~y, z = ~z,
    color = ~cell_type,
    colors = cell_colors,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 2)
  )
  
  # Function to generate sphere mesh
  sphere_mesh <- function(x0, y0, z0, r, n = 30) {
    theta <- seq(0, 2*pi, length.out = n)
    phi <- seq(0, pi, length.out = n)
    theta <- rep(theta, each = n)
    phi <- rep(phi, times = n)
    
    x <- x0 + r * sin(phi) * cos(theta)
    y <- y0 + r * sin(phi) * sin(theta)
    z <- z0 + r * cos(phi)
    
    list(x = matrix(x, n, n), y = matrix(y, n, n), z = matrix(z, n, n))
  }
  
  # Add spheres to the plot
  for (i in 1:nrow(plaque_data)) {
    s <- sphere_mesh(plaque_data$x[i], plaque_data$y[i], plaque_data$z[i], plaque_data$radius[i])
    
    p <- add_surface(
      p,
      x = s$x,
      y = s$y,
      z = s$z,
      opacity = 0.3,
      showscale = FALSE,
      surfacecolor = matrix(plaque_data$radius[i], nrow = nrow(s$x), ncol = ncol(s$x)),
      colorscale = list(c(0, 1), c("gray", "gray"))
    )
  }
  
}


