library(SpatialExperiment)
library(Giotto)
library(plyr)
library(ggplot2)
library(dplyr)
library(plotly)

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/basic_plaque_3d_overlays/"
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

# add sample id to centroids
plaque_cents <- lapply(sample_ids, function(id) {
  
  data <- plaque_cents[[id]]
  data$sample_id <- id
  
  return(data)
  
})
plaque_cents_df <- bind_rows(plaque_cents)

plaque_cents_df$size_scaled <- scales::rescale(log10(plaque_cents_df$total_um + 1), to = c(3, 18))

# determine mins and maxs
volume_min <- min(plaque_cents_df$size_scaled)
volume_max <- max(plaque_cents_df$size_scaled)


# plot out per sample cell color plots

for (sample_id in sample_ids) {
  
  giotto_obj <- m3d_samples[[sample_id]]
  plaque_meta <- plaque_cents_df[plaque_cents_df$sample_id == sample_id,]
  
  # Get your main cell/spot locations
  cell_locs <- getSpatialLocations(giotto_obj, spat_unit = "cell", output = "data.table")
  
  # Get plaque locations and metadata
  # plaque_locs <- giotto_obj@spatial_locs[["cell"]][["plaques"]]
  # plaque_meta <- giotto_obj@cell_metadata[["cell"]][["plaques"]]
  
  # Add cell types and colors to cell_locs
  # Assuming you have a vector of cell_types in the same order as cell_locs
  cell_locs$cell_type <- cell_color_meta[cell_locs$cell_ID,]$Cell.Type  # your cell type labels
  cell_locs$color <- cell_color_meta[cell_locs$cell_ID,]$cell_color     # your cell type colors
  
  # Log scale plaque sizes based on volume
  plaque_meta$size_scaled <- scales::rescale(log10(plaque_meta$total_um + 1), to = c(3, 18))
  
  # Create interactive 3D plot
  fig <- plot_ly()
  
  # Add cells by cell type (one trace per cell type for proper legend)
  unique_cell_types <- sort(unique(cell_locs$cell_type))
  
  for(ct in unique_cell_types) {
    ct_data <- cell_locs[cell_locs$cell_type == ct, ]
    
    fig <- fig %>%
      add_trace(
        data = as.data.frame(ct_data),
        x = ~sdimx, y = ~sdimy, z = ~sdimz,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 2,
          color = unique(ct_data$color),
          opacity = 0.7
        ),
        name = ct,
        visible = TRUE,  # All visible by default
        hoverinfo = "text",
        text = ~paste("Cell ID:", cell_ID, "<br>Cell Type:", cell_type)
      )
  }
  
  # Add plaques with log-scaled sizes
  fig <- fig %>%
    add_trace(
      data = plaque_meta,
      x = ~x_um, y = ~y_um, z = ~z_um,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = ~size_scaled,
        color = ~log10(total_um + 1),
        colorscale = "Greys",
        reversescale=T,
        cmin = volume_min,  # Fixed minimum
        cmax = volume_max,  # Fixed maximum
        showscale = TRUE,
        colorbar = list(
          title = "log10(Volume)<br>(μm³)",
          thickness = 20,
          len = 0.5,
          x = -0.1,
          y = 0.5,
          xanchor = "right",
          yanchor = "middle"
        ),
        line = list(color = "black", width = 1)
      ),
      name = "Plaques",
      visible = TRUE,
      hoverinfo = "text",
      text = ~paste(
        "Plaque ID:", plaque_id,
        "<br>Volume:", round(total_um, 4), "μm³",
        "<br>log10(Volume):", round(log10(total_um + 1), 2),
        "<br>Pixels:", total_pixels
      )
    )
  
  # Update layout with toggle mode for multi-select
  fig <- fig %>%
    layout(
      scene = list(
        xaxis = list(title = "X (μm)"),
        yaxis = list(title = "Y (μm)"),
        zaxis = list(title = "Z (μm)"),
        aspectmode = "data",
        camera = list(
          eye = list(x = 1.5, y = 1.5, z = 1.5)
        )
      ),
      title = paste0(sample_id, "\n3D Spatial Transcriptomics with Plaque Centroids<br><sub>Click legend items to show/hide cell types</sub>"),
      showlegend = TRUE,
      legend = list(
        x = 1.02,
        y = 0.5,
        xanchor = "left",
        yanchor = "middle",
        itemsizing = "constant",
        itemclick = "toggle",      # Click to toggle individual traces on/off
        itemdoubleclick = T,   # Disable double-click behavior
        bgcolor = "rgba(255, 255, 255, 0.8)",
        bordercolor = "gray",
        borderwidth = 1
      )
    )
  
  # Display the plot
  fig
  
  htmlwidgets::saveWidget(fig, paste0(out_dir, sample_id, ".cells_plaque_3d.html"), selfcontained = TRUE)
  
}




