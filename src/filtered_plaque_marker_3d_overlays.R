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

out_dir <- "results/filtered_plaque_marker_3d_overlays/"
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

# add sample id to centroids and filter plaques based on size and location
plaque_cents <- lapply(sample_ids, function(id) {
  
  data <- plaque_cents[[id]]
  data$sample_id <- id
  
  data$plaque_id <- paste0(data$sample_id,
                           "|", data$plaque_id)
  
  data <- data[data$total_um > 20 &
                 data$z_um < 30 &
                 data$total_um < 100849,]
  
  return(data)
  
})
plaque_cents_df <- bind_rows(plaque_cents)
# scale them all together
plaque_cents_df$size_scaled <- scales::rescale(log10(plaque_cents_df$total_um + 1), to = c(3, 18))

# find the markers
endo_markers <- c("Cldn5","Igf2","Rgs5")

for (id in sample_ids) {
  
  g <- m3d_samples[[id]]
  
  # normalize expression values
  g <- processExpression(g, normParam("default"), expression_values = "raw")
  
  expr_values <- getExpression(g,
                               output="matrix",
                               values = "normalized")[endo_markers,]
  
  thresholds <- apply(expr_values, 1, function(x) quantile(x, 0.50, na.rm=T))
  
  # Step 4: Identify cells with high expression of ALL three genes
  high_expr_cells <- colnames(expr_values)[
    expr_values["Cldn5", ] > thresholds["Cldn5"] &
      expr_values["Igf2", ] > thresholds["Igf2"] &
      expr_values["Rgs5", ] > thresholds["Rgs5"]
  ]
  
  high_expr_cells <- high_expr_cells[!is.na(high_expr_cells)]
  
  
  sample_cell_metadata <- pDataDT(g)
  sample_cell_metadata$marker_cluster <- ifelse(
    sample_cell_metadata$cell_ID %in% high_expr_cells,
    "Endothelia Markers",
    "Other"
  )
  
  # Add this metadata back to the Giotto object
  g <- addCellMetadata(g,
                             new_metadata = sample_cell_metadata,
                             by_column = TRUE,
                             column_cell_ID = "cell_ID")
  
  spatPlot2D(g,
             cell_color = "marker_cluster",
             select_cell_groups = "Endothelia Markers",
             point_size = 1,
             other_cell_color = "grey",
             cell_color_code = c("Endothelia Markers" = "red"),
             show_legend = TRUE,
             legend_text = 10)
  
  # code for plotting
  plaque_meta <- plaque_cents_df[plaque_cents_df$sample_id == id,]
  
  sample_cell_metadata$marker_color <- ifelse(sample_cell_metadata$marker_cluster == "Endothelia Markers",
                                              "red",
                                              "lightgrey")
  
  sample_cell_metadata$marker_opacity <- ifelse(sample_cell_metadata$marker_cluster == "Endothelia Markers",
                                                0.8,
                                                0.3)
  
  cell_locs <- sample_cell_metadata[sample_cell_metadata$Z < 30,]
  
  # Create interactive 3D plot
  fig <- plot_ly()
  
  # Add cells by cell type (one trace per cell type for proper legend)
  unique_cell_types <- sort(unique(cell_locs$marker_cluster))
  
  for(ct in unique_cell_types) {
    ct_data <- cell_locs[cell_locs$marker_cluster == ct, ]
    
    fig <- fig %>%
      add_trace(
        data = as.data.frame(ct_data),
        x = ~X, y = ~Y, z = ~Z,
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = 4,
          color = unique(ct_data$marker_color),
          opacity = unique(ct_data$marker_opacity)
        ),
        name = ct,
        visible = TRUE,  # All visible by default
        hoverinfo = "text",
        text = ~paste("Cell ID:", cell_ID, "<br>Cell Type:", marker_cluster)
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
        #        color = ~log10(total_um + 1),
        #        colorscale = "Greys",
        #        reversescale=T,
        color="black",
        #cmin = volume_min,  # Fixed minimum
        #cmax = volume_max,  # Fixed maximum
        #showscale = TRUE,
        #colorbar = list(
        #  title = "log10(Volume)<br>(μm³)",
        #  thickness = 20,
        #  len = 0.5,
        #  x = -0.1,
        #  y = 0.5,
        #  xanchor = "right",
        #  yanchor = "middle"
        #),
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
      title = paste0(id, "\n3D Spatial Transcriptomics with Plaque Centroids<br><sub>Click legend items to show/hide cell types</sub>"),
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
  
  htmlwidgets::saveWidget(fig, paste0(out_dir, id, ".cells_plaque_3d.html"), selfcontained = TRUE)
  
}

# zip it up
curr_wd <- getwd()

setwd(out_dir)

zip(zipfile = "endo_marker_plaque_plots.zip",
    files=list.files(".", pattern="\\.html$"))

setwd(curr_wd)

