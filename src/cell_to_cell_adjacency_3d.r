library(SpatialExperiment)
library(Giotto)
library(plyr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(future)
library(future.apply)
library(patchwork)

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/cell_to_cell_adjacency_3d/"
dir.create(out_dir, recursive = T, showWarnings = F)

results_folder = out_dir

instrs <- createGiottoInstructions(
  save_dir = results_folder,
  save_plot = TRUE,
  show_plot = FALSE,
  return_plot = FALSE
)

# read in spatil experiment converted file

m3d_se <- readRDS("../WT_v_APP_PS19/results/anndata_objects/all_cells_integrated.spatial_exp.RDS")
# pull out cell metadata
cell_metadata <- as.data.frame(colData(m3d_se))
cell_metadata$cell_id <- rownames(cell_metadata)

# colors for consistent plotting
cell_colors <- metadata(m3d_se)$"Cell Type_colors"

cell_type_counts <- as.data.frame(table(cell_metadata$Cell.Type))
cell_type_counts <- cell_type_counts[order(cell_type_counts$Freq),]
cell_type_counts$Var1 <- factor(as.character(cell_type_counts$Var1),
                                levels=cell_type_counts$Var1)

ggplot(cell_type_counts,
       aes(x=Freq, y=Var1)) +
  geom_bar(stat="identity", color="white", fill="maroon") +
  theme_bw() +
  labs(x="Cell Count", y="Cell Type")
ggsave(paste0(out_dir, "total_cell_type_counts.bar.png"), width=6, height=5)

# convert to giotto object ...
m3d_g <- spatialExperimentToGiotto(m3d_se)
# add in z axis, since it does not get added
m3d_g@spatial_locs$cell$raw$sdimz <- spatialCoords(m3d_se)[,"Z"]

# plot instructions for Giotto
instructions(m3d_g) <- instrs

# make a sample name that has the condition in it included
cell_metadata$sample_cond <- paste0(cell_metadata$ID, " ",cell_metadata$Condition)

sample_ids <- unique(cell_metadata$sample_cond)

sample_order <- c("C164B WT","C158B APPPS19",
                  "C166A WT", "C165 APPPS19")

# separate out samples by subsetting
m3d_samples <- lapply(sample_ids, function(id) {
  
  cell_ids <- rownames(cell_metadata[cell_metadata$sample_cond == id,])
  
  subset(m3d_g, cell_ids=cell_ids)
  
})
names(m3d_samples) <- sample_ids

# connection analysis per sample
cell_connections <- lapply(sample_ids, function(id) {
  
  print(id)
  
  g <- m3d_samples[[id]]
  
  # use Delaundy triangulation to determine spatial network
  g <- createSpatialNetwork(g, method = "Delaunay",
                            delaunay_method = "delaunayn_geometry",
                            name = "spatial_network")
  
  # pull out spatial network
  sp_net <- g@spatial_network$cell$spatial_network
  # convert to dataframe
  sp_graph <- as.data.frame(sp_net@networkDT)
  
  # graph is uni-directional, need to be bi-directional for better counts
  sp_graph_simple <- sp_graph[,c("from","to","distance","weight")]
  sp_graph_simple_rev <- sp_graph_simple
  colnames(sp_graph_simple_rev) <- c("to","from","distance","weight")
  
  sp_graph_simple_both <- rbind(sp_graph_simple,
                                sp_graph_simple_rev)
  
  sp_graph_simple_both <- unique(sp_graph_simple_both)
  
  
  # merge in metadata to determine cell types
  sp_graph_simple_both <- merge(sp_graph_simple_both,
                                cell_metadata[,c("cell_id","Cell.Type","Condition")],
                                by.x="from",
                                by.y="cell_id")
  
  sp_graph_simple_both <- merge(sp_graph_simple_both,
                                cell_metadata[,c("cell_id","Cell.Type", "Condition")],
                                by.x="to",
                                by.y="cell_id",
                                suffixes=c("_from","_to"))
  
  sp_graph_simple_both$sample_id <- id
  
  
  return(sp_graph_simple_both)
})
names(cell_connections) <- sample_ids

# count up totals per cell
cell_connect_totals <- lapply(sample_ids, function(id) {
  
  cell_connect <- cell_connections[[id]]
  
  # sum up connected cell count
  cell_connect_total <- ddply(cell_connect,
                              .(from, Cell.Type_from),
                              summarise,
                              cell_count=length(unique(to)))
  # add in sample metadata
  cell_connect_total$sample <- id
  
  return(cell_connect_total)
})
cell_connect_totals <- bind_rows(cell_connect_totals)

# calc means and medians
cell_connect_stats <- ddply(cell_connect_totals,
                            .(sample),
                            summarise,
                            mean=mean(cell_count),
                            median=median(cell_count),
                            sd=sd(cell_count))

cell_connect_stats$stat <- paste0("Mean: ", round(cell_connect_stats$mean, digits=2),
                                  "\nMedian: ", round(cell_connect_stats$median, digits=2),
                                  "\nStdDev: ", round(cell_connect_stats$sd, digits=2))

cell_connect_totals$sample <- factor(cell_connect_totals$sample,
                                     levels=sample_order)
cell_connect_stats$sample <- factor(cell_connect_stats$sample,
                                    levels=sample_order)

# plot histogram of counts across samples
ggplot(cell_connect_totals,
       aes(cell_count)) +
  geom_histogram(bins=31, color="white", fill="darkblue") +
  theme_bw() +
  labs(x="Number of Connected Cells",
       y="Frequency") +
  facet_wrap(~ sample, ncol=2) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30, 35, 40)) +
  geom_label(data=cell_connect_stats, aes(label=stat), x=30, y=3100, size=3)
ggsave(paste0(out_dir, "sample_cell_connectivity.histogram.png"), width=6, height=5)


# make one big adjacency matrix, to determine cell type order
cell_connections_big <- bind_rows(cell_connections)

cell_type_connect_big <- ddply(cell_connections_big,
                               .(Cell.Type_from,Cell.Type_to),
                               summarise,
                               cell_count=length(from))  

cell_type_connect_big$scaled_count <- scale(log(cell_type_connect_big$cell_count))

connect_big_matrix <- acast(cell_type_connect_big,
                            Cell.Type_from ~ Cell.Type_to,
                            value.var = "scaled_count")

pdf(paste0(out_dir, "all_samples.cell_type_adj_matrix.pdf"), width=8, height=7)
print(Heatmap(connect_big_matrix,
              col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
              name = "Normalized\nEdge Count",
              column_title="All Samples"))
dev.off()

# hard code cell type order across samples
cell_type_order <- c("Astrocyte",
                     "Oligodendrocyte",
                     "E_Homer2_Pou3f1",
                     "I_Npy_Lhx6",
                     "E_Tbr1_Id2",
                     "Endothelia",
                     "E_Rora_Zic1",
                     "I_Chrm2_Pvalb",
                     "Oligodendrocyte progenitor cell",
                     "Microglia",
                     "E_Tshz2_Lamp5",
                     "Hybrid_EI_neurons",
                     "E_Cck_Tcf4",
                     "E_Pou3f1_Cck",
                     "E_Zic1_Lef1",
                     "E_Baiap3_Tmem163",
                     "E_Zic1_Syt9",
                     "Choroid plexus",
                     "Ependymal",
                     "I_Ppp1r1b_Gpr88")

# make per sample heatmaps

for (id in sample_ids) {
  
  print(id)
  
  cell_connect <- cell_connections[[id]]
  
  # sum up cell 2 cell edge counts
  cell_type_connect_total <- ddply(cell_connect,
                                   .(Cell.Type_from,Cell.Type_to),
                                   summarise,
                                   cell_count=length(from))  
  
  # log scale and z-score
  cell_type_connect_total$scaled_count <- scale(log(cell_type_connect_total$cell_count))
  
  connect_matrix <- acast(cell_type_connect_total,
                          Cell.Type_from ~ Cell.Type_to,
                          value.var = "scaled_count")
  
  pdf(paste0(out_dir, id, ".cell_type_adj_matrix.pdf"), width=8, height=7)
  print(Heatmap(connect_matrix,
                col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
                name = "Normalized\nEdge Count",
                column_title=id,
                cluster_rows = F, cluster_columns = F,
                row_order = cell_type_order,
                column_order = cell_type_order))
  dev.off()
  
}

# make per condition heatmaps
conditions <- unique(cell_metadata$Condition)

for (cond in conditions) {
  
  print(cond)
  
  cond_ids <- unique(cell_metadata[cell_metadata$Condition == cond,]$sample_cond)
  
  cell_connect <- bind_rows(cell_connections[cond_ids])
  
  # sum up cell 2 cell edge counts
  cell_type_connect_total <- ddply(cell_connect,
                                   .(Cell.Type_from,Cell.Type_to),
                                   summarise,
                                   cell_count=length(from))  
  
  # log scale and z-score
  cell_type_connect_total$scaled_count <- scale(log(cell_type_connect_total$cell_count))
  
  connect_matrix <- acast(cell_type_connect_total,
                          Cell.Type_from ~ Cell.Type_to,
                          value.var = "scaled_count")
  
  
  pdf(paste0(out_dir, cond, ".cell_type_adj_matrix.pdf"), width=8, height=7)
  print(Heatmap(connect_matrix,
                col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
                name = "Normalized\nEdge Count",
                column_title=cond,
                cluster_rows = F, cluster_columns = F,
                row_order = cell_type_order,
                column_order = cell_type_order))
  dev.off()
  
}


saveRDS(bind_rows(cell_connections), file=paste0(out_dir, "sample_delaunay_connections.RDS"))

for (id in sample_ids) {
  
  g <- m3d_samples[[id]]
  
  g@spatial_locs$cell$raw$sdimy <- 
    g@spatial_locs$cell$raw$sdimy * -1
  
  spatPlot2D(g,
            cell_color="Cell.Type",
            cell_color_code=cell_colors,
            title=id,
            point_size = 1,
            save_param=list(save_name=paste0(id, ".2d_splot")))
  
  
}


# try out stellaromics code
# cluster_dist <- function(data_df) {
#   # Ensure required columns exist
#   if (!all(c("X_um", "Y_um", "Z_um", "labeled") %in% names(data_df))) {
#     stop("Dataframe must contain 'X_um', 'Y_um', 'Z_um', and 'labeled' columns.")
#   }
#   
#   # Get unique cluster labels
#   cluster_labels <- unique(data_df$labeled)
#   num_clusters <- length(cluster_labels)
#   
#   # Use future_lapply for parallelizing the outer loop.
#   # Each element in the result_list will be a vector representing a row of the distance matrix.
#   result_list <- future_lapply(1:num_clusters, function(i) {
#     source_cluster_name <- cluster_labels[i]
#     source_points <- data_df[data_df$labeled == source_cluster_name, c("X_um", "Y_um", "Z_um")]
#     
#     # Vector to store the results for the current source cluster (i.e., a row of the final matrix)
#     row_results <- numeric(num_clusters)
#     
#     # Loop through each cluster as the "target"
#     for (j in 1:num_clusters) {
#       target_cluster_name <- cluster_labels[j]
#       target_points <- data_df[data_df$labeled == target_cluster_name, c("X_um", "Y_um", "Z_um")]
#       
#       # If source and target clusters are the same, calculate average minimum distance
#       # to *other* cells within the same cluster.
#       if (source_cluster_name == target_cluster_name) {
#         # If there's only one cell in the cluster, cannot find a 'minimum distance to another cell'
#         if (nrow(source_points) <= 1) {
#           row_results[j] <- NA
#           next # Skip to the next iteration
#         }
#         
#         min_dists_from_source_to_self_cluster <- numeric(nrow(source_points))
#         
#         # For each cell in the source cluster, find its minimum distance to any *other* cell in the same cluster
#         for (s_idx in 1:nrow(source_points)) {
#           current_source_cell_coords <- as.numeric(source_points[s_idx, ])
#           
#           # Get all other points in the same cluster, excluding the current cell itself
#           other_cells_in_same_cluster <- source_points[-s_idx, , drop = FALSE]
#           
#           # Calculate Euclidean distances from the current source cell to these other points
#           temp_coords_matrix <- rbind(current_source_cell_coords, as.matrix(other_cells_in_same_cluster))
#           all_pairwise_dists <- as.matrix(dist(temp_coords_matrix))
#           
#           # The distances from the first point (our current source cell) to all other points
#           # (which are the other cells in the same cluster) are in the first row, excluding the first element.
#           distances_from_cell_to_other_cells_in_same_cluster <- all_pairwise_dists[1, -1]
#           
#           # Find the minimum of these distances
#           min_dist_for_current_cell <- min(distances_from_cell_to_other_cells_in_same_cluster)
#           
#           # Store this minimum distance
#           min_dists_from_source_to_self_cluster[s_idx] <- min_dist_for_current_cell
#         }
#         # Calculate the average of these minimum distances
#         row_results[j] <- mean(min_dists_from_source_to_self_cluster)
#         
#       } else {
#         # Handle cases where a cluster might be empty
#         if (nrow(source_points) == 0 || nrow(target_points) == 0) {
#           row_results[j] <- NA
#           next # Skip to the next iteration if either cluster is empty
#         }
#         
#         # Vector to store the minimum distance from each source cell to the target cluster
#         min_dists_from_source_to_target <- numeric(nrow(source_points))
#         
#         # For each cell in the source cluster, find its minimum distance to any cell in the target cluster
#         for (s_idx in 1:nrow(source_points)) {
#           current_source_cell_coords <- as.numeric(source_points[s_idx, ])
#           
#           # Calculate Euclidean distances from the current source cell to all target points
#           # We combine the current cell's coordinates with the target cluster's points
#           # into a single matrix. 'dist()' then calculates all pairwise distances.
#           temp_coords_matrix <- rbind(current_source_cell_coords, as.matrix(target_points))
#           all_pairwise_dists <- as.matrix(dist(temp_coords_matrix))
#           
#           # The distances from the first point (our current source cell) to all other points
#           # (which are the target cluster points) are located in the first row of
#           # 'all_pairwise_dists' matrix, excluding the first element (which is the
#           # distance from the cell to itself, always 0).
#           distances_from_cell_to_target_cluster <- all_pairwise_dists[1, -1]
#           
#           # Find the minimum of these distances
#           min_dist_for_current_cell <- min(distances_from_cell_to_target_cluster)
#           
#           # Store this minimum distance
#           min_dists_from_source_to_target[s_idx] <- min_dist_for_current_cell
#         }
#         
#         # Calculate the average of these minimum distances
#         row_results[j] <- mean(min_dists_from_source_to_target)
#       }
#     }
#     return(row_results)
#   })
#   
#   # Combine the list of row results into the final distance matrix
#   distance_matrix <- do.call(rbind, result_list)
#   colnames(distance_matrix) <- cluster_labels
#   rownames(distance_matrix) <- cluster_labels
#   
#   return(distance_matrix)
# }
# 
# plan(multisession, workers = 4)
# 
# 
# sample_mat <- lapply(m3d_samples, function(g) {
#   
#   # metadata
#   cluster_locs <- data.frame(X_um=g@cell_metadata$cell$rna$X,
#                              Y_um=g@cell_metadata$cell$rna$Y,
#                              Z_um=g@cell_metadata$cell$rna$Z,
#                              labeled=g@cell_metadata$cell$rna$Cell.Type)
#   
#   mat <- cluster_dist(cluster_locs)
#   
# })


