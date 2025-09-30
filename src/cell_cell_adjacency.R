
library(Seurat)
library(tidyverse)
library(future)
library(future.apply)
library(patchwork)
source('')

analysis_dir = ''
dir.create(analysis_dir, recursive = T)
setwd(analysis_dir)


# Functions ---------------------------------------------------------------

cluster_dist <- function(data_df) {
  # Ensure required columns exist
  if (!all(c("X_um", "Y_um", "Z_um", "labeled") %in% names(data_df))) {
    stop("Dataframe must contain 'X_um', 'Y_um', 'Z_um', and 'labeled' columns.")
  }
  
  # Get unique cluster labels
  cluster_labels <- unique(data_df$labeled)
  num_clusters <- length(cluster_labels)
  
  # Use future_lapply for parallelizing the outer loop.
  # Each element in the result_list will be a vector representing a row of the distance matrix.
  result_list <- future_lapply(1:num_clusters, function(i) {
    source_cluster_name <- cluster_labels[i]
    source_points <- data_df[data_df$labeled == source_cluster_name, c("X_um", "Y_um", "Z_um")]
    
    # Vector to store the results for the current source cluster (i.e., a row of the final matrix)
    row_results <- numeric(num_clusters)
    
    # Loop through each cluster as the "target"
    for (j in 1:num_clusters) {
      target_cluster_name <- cluster_labels[j]
      target_points <- data_df[data_df$labeled == target_cluster_name, c("X_um", "Y_um", "Z_um")]
      
      # If source and target clusters are the same, calculate average minimum distance
      # to *other* cells within the same cluster.
      if (source_cluster_name == target_cluster_name) {
        # If there's only one cell in the cluster, cannot find a 'minimum distance to another cell'
        if (nrow(source_points) <= 1) {
          row_results[j] <- NA
          next # Skip to the next iteration
        }
        
        min_dists_from_source_to_self_cluster <- numeric(nrow(source_points))
        
        # For each cell in the source cluster, find its minimum distance to any *other* cell in the same cluster
        for (s_idx in 1:nrow(source_points)) {
          current_source_cell_coords <- as.numeric(source_points[s_idx, ])
          
          # Get all other points in the same cluster, excluding the current cell itself
          other_cells_in_same_cluster <- source_points[-s_idx, , drop = FALSE]
          
          # Calculate Euclidean distances from the current source cell to these other points
          temp_coords_matrix <- rbind(current_source_cell_coords, as.matrix(other_cells_in_same_cluster))
          all_pairwise_dists <- as.matrix(dist(temp_coords_matrix))
          
          # The distances from the first point (our current source cell) to all other points
          # (which are the other cells in the same cluster) are in the first row, excluding the first element.
          distances_from_cell_to_other_cells_in_same_cluster <- all_pairwise_dists[1, -1]
          
          # Find the minimum of these distances
          min_dist_for_current_cell <- min(distances_from_cell_to_other_cells_in_same_cluster)
          
          # Store this minimum distance
          min_dists_from_source_to_self_cluster[s_idx] <- min_dist_for_current_cell
        }
        # Calculate the average of these minimum distances
        row_results[j] <- mean(min_dists_from_source_to_self_cluster)
        
      } else {
        # Handle cases where a cluster might be empty
        if (nrow(source_points) == 0 || nrow(target_points) == 0) {
          row_results[j] <- NA
          next # Skip to the next iteration if either cluster is empty
        }
        
        # Vector to store the minimum distance from each source cell to the target cluster
        min_dists_from_source_to_target <- numeric(nrow(source_points))
        
        # For each cell in the source cluster, find its minimum distance to any cell in the target cluster
        for (s_idx in 1:nrow(source_points)) {
          current_source_cell_coords <- as.numeric(source_points[s_idx, ])
          
          # Calculate Euclidean distances from the current source cell to all target points
          # We combine the current cell's coordinates with the target cluster's points
          # into a single matrix. 'dist()' then calculates all pairwise distances.
          temp_coords_matrix <- rbind(current_source_cell_coords, as.matrix(target_points))
          all_pairwise_dists <- as.matrix(dist(temp_coords_matrix))
          
          # The distances from the first point (our current source cell) to all other points
          # (which are the target cluster points) are located in the first row of
          # 'all_pairwise_dists' matrix, excluding the first element (which is the
          # distance from the cell to itself, always 0).
          distances_from_cell_to_target_cluster <- all_pairwise_dists[1, -1]
          
          # Find the minimum of these distances
          min_dist_for_current_cell <- min(distances_from_cell_to_target_cluster)
          
          # Store this minimum distance
          min_dists_from_source_to_target[s_idx] <- min_dist_for_current_cell
        }
        
        # Calculate the average of these minimum distances
        row_results[j] <- mean(min_dists_from_source_to_target)
      }
    }
    return(row_results)
  })
  
  # Combine the list of row results into the final distance matrix
  distance_matrix <- do.call(rbind, result_list)
  colnames(distance_matrix) <- cluster_labels
  rownames(distance_matrix) <- cluster_labels
  
  return(distance_matrix)
}

dist_dotplot <- function(distance_matrix, 
                         plot_title = "Average Minimum Distance Between Clusters", 
                         cluster_order = NULL,
                         size_limits = NULL) {

  # Get the order of cluster labels from the matrix column names
  if(is.null(cluster_order)){
    cluster_order <- colnames(distance_matrix) 
  }
  
  # Convert the distance matrix to a long-format dataframe using melt
  # 'Var1' will be Source_Cluster, 'Var2' will be Target_Cluster, 'value' will be Distance
  df_long <- reshape2::melt(distance_matrix, varnames = c("Source_Cluster", "Target_Cluster"), value.name = "Distance")
  
  # Filter out NA values if any (e.g., if a cluster was empty and resulted in NA distances)
  df_long <- na.omit(df_long)
  
  # Ensure the order of clusters on both axes is the same as in the matrix
  df_long$Source_Cluster <- factor(df_long$Source_Cluster, levels = cluster_order)
  df_long$Target_Cluster <- factor(df_long$Target_Cluster, levels = rev(cluster_order))
  
  # Create the ggplot dot plot
  p <- ggplot(df_long, aes(x = Target_Cluster, y = Source_Cluster, size = Distance, fill = Distance)) +
    geom_point(shape = 21, stroke = .25) + # Use shape 21 for filled circles with no border
    scale_fill_gradientn(colors = pals::magma(10), name = "Avg Min Distance", limits = size_limits) + # Color scale for distance
    scale_size_area(max_size = 20, name = "Avg Min Distance (uM)", limits = size_limits) + # Size scale for distance, max_size for aesthetics
    labs(
      title = plot_title,
      x = "Cell Type",
      y = "Cell Type"
    ) +
    theme_minimal() + # A clean theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # Rotate x-axis labels for readability
      panel.grid.major = element_line(color = "grey90", linetype = "dashed"), # Lighter grid lines
      panel.grid.minor = element_blank(), # Remove minor grid lines
      plot.title = element_text(hjust = 0.5, face = "bold") # Center and bold title
    )
  
  return(p)
}


# Analysis ----------------------------------------------------------------

# Use 100um PVN Dataset
mb = readRDS('S:/PVN_Seurat_Object.Rds')

# 100um
mb_100 = mb

# 5um subset
mb_5 = subset(mb, Z_um > 60 & Z_um < 65)
mb_10 = subset(mb, Z_um > 55 & Z_um < 65)

rm(mb)

# 5um cell density per 100um^2
minx = min(mb_5$X_um)
miny = min(mb_5$Y_um)
maxx = max(mb_5$X_um)
maxy = max(mb_5$Y_um)

cell_locs = mb_5@meta.data %>% dplyr::select(X_um, Y_um)

cells_in = list()
for(x in head(seq(minx, maxx, by = 100), -1)){
  for(y in head(seq(miny, maxy, by = 100), -1)){
    
    cells_in[paste0(round(x),'_',round(y))] = cell_locs %>% 
      filter(X_um > x, X_um < (x+100), Y_um > y, Y_um < (y+100)) %>%
      nrow
    
  }
}
um5_density = cells_in %>% unlist() %>% mean(na.rm = T)
# 4.253247 cells/100um^2

# Make a 2D plot of all clusters in 5um
p = PlotSlice(mb_5, cluster.level = 'RNA_snn_res.2') + ggtitle('All Clusters')
ggsave('All Clusters 5_um.png', p, width = 8, height = 6, bg = 'white', scale = 1.5)

# Do 10um as well. This data is using centroids, but a Xenium run would include
# cells who may only have 1/4 of the cell in frame (the true cell centroid may not be in frame)
p = PlotSlice(mb_10, cluster.level = 'RNA_snn_res.2') + ggtitle('All Clusters')
ggsave('All Clusters 10_um.png', p, width = 8, height = 6, bg = 'white', scale = 1.5)


# 100um cell density per 100um^2
minx = min(mb_100$X_um)
miny = min(mb_100$Y_um)
maxx = max(mb_100$X_um)
maxy = max(mb_100$Y_um)

cell_locs = mb_100@meta.data %>% dplyr::select(X_um, Y_um)

cells_in = list()
for(x in head(seq(minx, maxx, by = 100), -1)){
  for(y in head(seq(miny, maxy, by = 100), -1)){
    
    cells_in[paste0(round(x),'_',round(y))] = cell_locs %>% 
      filter(X_um > x, X_um < (x+100), Y_um > y, Y_um < (y+100)) %>%
      nrow
    
  }
}
um100_density = cells_in %>% unlist() %>% mean(na.rm = T)
# 65.74675 cells/100um^2


# Make a barplot
um_densities = data.frame('um_5' = 4.25,'um_100' = 65.75) %>% reshape2::melt()
p = um_densities %>% ggplot(aes(x = variable, y = value, fill = variable)) + 
  geom_col() + scale_fill_manual(values = c('#123050','#5db7c0')) +
  theme_minimal() + theme(legend.position = 'None') +
  labs(x = 'Thickness', y = 'Cells/100um^2')
ggsave(p, file = 'Cells per 100um2.png', width = 8, height = 6, bg = 'white', scale = .8)


# Major Cell Class Annotation ---------------------------------------------

pvn_markers = FindAllMarkers(mb_100, logfc.threshold = .25)
pvn_markers %>% filter(avg_log2FC > 0, p_val_adj < .05, 
                       gene %in% c('Opalin','Gja1','Hdc','P2ry12','Rgs5'))

mb_100 = RenameIdents(mb_100, c('2' = 'Oligodendrocyte', '6' = 'Astrocyte',
                                '7' = 'Astrocyte','9' = 'Endothelia',
                                '13' = 'Endothelia','15' = 'Microglia',
                                '40', 'Ependymal', '0' = 'I_Neuron',
                                '23' = 'I_Neuron', '37' = 'I_Neuron',
                                '1' = 'E_Neuron', '41' = 'E_Neuron',
                                '4' = 'I_Neuron','17' = 'I_Neuron',
                                '22' = 'I_Neuron', '30' = 'I_Neuron',
                                '29' = 'I_Neuron', '45'= 'I_Neuron',
                                '26' = 'I_Neuron', '47' = 'E_Neuron',
                                '12' = 'E_Neuron','27' = 'E_Neuron',
                                '36' = 'E_Neuron', '43'= 'E_Neuron',
                                '40' = 'E_Neuron', '34'= 'E_Neuron',
                                '24'= 'E_Neuron', '35'= 'E_Neuron'))
mb_100$labeled = Idents(mb_100)

saveRDS(mb_100, file = 'PVN 100um Seurat.Rds')


# Plot specific populations -----------------------------------------------


for(ct in c('Oligodendrocyte', 'Astrocyte','Endothelia','Microglia','I_Neuron','E_Neuron')){
  
  p = PlotSlice(mb_100, cluster = ct, cluster.level = 'labeled', pt.size = 1) + ggtitle(ct)
  ggsave(p, file = paste0(ct, '.png'), width = 10, height = 10, bg = 'white')
  
}

# Nearest Neighbor Analysis -----------------------------------------------

# Get identities and xyz locations
cell_locs = mb_100@meta.data %>%
  filter(labeled %in% c('Oligodendrocyte', 'Astrocyte','Endothelia',
                        'Microglia','I_Neuron','E_Neuron')) %>%
  dplyr::select(X_um, Y_um, Z_um, labeled)

cell_locs_5um = mb_100@meta.data %>% filter(Z_um > 50, Z_um <55)
  filter(labeled %in% c('Oligodendrocyte', 'Astrocyte','Endothelia',
                        'Microglia','I_Neuron','E_Neuron')) %>%
  dplyr::select(X_um, Y_um, Z_um, labeled)


plan(sequential)
gc()
plan(multisession)


dists_100um = cluster_dist(cell_locs)
dists_5um = cluster_dist(cell_locs_5um)

saveRDS(dists_100um, file = 'Major Cell Class Neighborhood Dists in 100uM.Rds')
saveRDS(dists_5um, file = 'Major Cell Class Neighborhood Dists in 5uM.Rds')



# Make dotplot/heatmaps of distance ---------------------------------------

# Get max range of dists so everything is on a common scale between the two plots
min_dist = min(cbind(dists_100um, dists_5um))
max_dist = max(cbind(dists_100um, dists_5um))

p100 = dist_dotplot(dists_100um, plot_title = "Average Minimum Distance Between Clusters (100um)",
                 cluster_order = c('Oligodendrocyte', 'Astrocyte','Endothelia',
                                   'Microglia','I_Neuron','E_Neuron'),
                 size_limits = c(min_dist, max_dist))
ggsave(p100, file = '100um Cell Type Distance Dotplot.png', width = 12, height = 9, bg = 'white', scale = .65)

p5 = dist_dotplot(dists_5um, plot_title = "Average Minimum Distance Between Clusters (5um)",
                 cluster_order = c('Oligodendrocyte', 'Astrocyte','Endothelia',
                                   'Microglia','I_Neuron','E_Neuron'),
                 size_limits = c(min_dist, max_dist))
ggsave(p5, file = '5um Cell Type Distance Dotplot.png', width = 12, height = 9, bg = 'white', scale = .65)


p = p5|p100 + plot_layout(guides = 'collect')
ggsave(p, file = '5 and 100um Cell Type Distance Dotplot.png', width = 16, height = 9, bg = 'white', scale = .65)
