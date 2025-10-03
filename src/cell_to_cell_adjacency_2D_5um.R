library(SpatialExperiment)
library(Giotto)
library(plyr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(ggsci)

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/cell_to_cell_adjacency_2D_5um/"
dir.create(out_dir, recursive = T, showWarnings = F)


# read in spatial experiment converted file

m3d_se <- readRDS("../WT_v_APP_PS19/results/anndata_objects/all_cells_integrated.spatial_exp.RDS")
# pull out cell metadata
cell_metadata <- as.data.frame(colData(m3d_se))
cell_metadata$cell_id <- rownames(cell_metadata)
cell_metadata$ID <- paste(cell_metadata$ID, cell_metadata$Condition)

# convert to giotto object ...
m3d_g <- spatialExperimentToGiotto(m3d_se)
# add in z axis, since it does not get added
m3d_g@spatial_locs$cell$raw$sdimz <- spatialCoords(m3d_se)[,"Z"]

# plot instructions for Giotto
instrs <- createGiottoInstructions(
  save_dir = out_dir,
  save_plot = TRUE,
  show_plot = FALSE,
  return_plot = FALSE
)
instructions(m3d_g) <- instrs

cell_metadata$z_int <- floor(cell_metadata$Z)

ggplot(cell_metadata,
       aes(x=z_int)) +
  geom_histogram(fill="maroon", color="white", binwidth = 1) +
  labs(x="Z Coordinate", y="Cell Count") +
  theme_bw()
ggsave(paste0(out_dir, "z_coord_histogram.png"), width=5, height=4)

# add slice information to metadata
cell_metadata$slice <- "unused"

cell_metadata[cell_metadata$Z >= 10 &
                cell_metadata$Z <= 15,]$slice <- "Slice 1"
cell_metadata[cell_metadata$Z >= 20 &
                cell_metadata$Z <= 25,]$slice <- "Slice 2"
cell_metadata[cell_metadata$Z >= 30 &
                cell_metadata$Z <= 35,]$slice <- "Slice 3"

ggplot(cell_metadata,
       aes(x=z_int, fill=slice)) +
  geom_histogram(color="white", binwidth = 1) +
  labs(x="Z Coordinate", y="Cell Count") +
  theme_bw() +
  scale_fill_manual(values=c(pal_futurama()(3), "grey"))
ggsave(paste0(out_dir, "z_coord_histogram.slices.png"), width=6, height=4)

ggplot(cell_metadata[cell_metadata$slice != "unused",],
       aes(x=slice, fill=ID)) +
  geom_bar(color="white") +
  theme_bw() +
  scale_fill_nejm() +
  labs(x=NULL, y="Cell Count")
ggsave(paste0(out_dir, "slice_sample_cell_counts.bar.png"), width=5, height=4)

# colors for consistent plotting
cell_colors <- metadata(m3d_se)$"Cell Type_colors"

# add in slice metadata
m3d_g$slice <- cell_metadata$slice

spatPlot3D(m3d_g,
           cell_color="slice",
           cell_color_code = c(pal_futurama()(3), "grey"),
           point_size = 5,
           save_param = list(save_name="slice_3d_plot"))

# split data into slices

slices <- paste0("Slice ", 1:3)

slice_list <- lapply(slices, function(slice) {
  
  subset(m3d_g, cell_ids=cell_metadata[cell_metadata$slice == slice,]$cell_id)
  
})
names(slice_list) <- slices

# cell to cell adjacency per sample
sample_ids <- unique(cell_metadata$ID)
slices <- names(slice_list)

cell_connections <- lapply(sample_ids, function(id) {
  
  print(id)
  slice_cell_connections <- lapply(slices, function(slice) {
    
    print(slice)
    # pull out slice
    slice_g <- slice_list[[slice]]
    
    # determine cells of current sample
    meta_subset <- cell_metadata[cell_metadata$slice == slice &
                                   cell_metadata$ID == id,]
    # pull out sample
    sample_slice_g <- subset(slice_g, cell_id=meta_subset$cell_id)
    
    # use Delaundy triangulation to determine spatial network
    sample_slice_g <- createSpatialNetwork(sample_slice_g, 
                                           method = "Delaunay",
                                           delaunay_method = "delaunayn_geometry",
                                           name = "spatial_network")
    
    # pull out spatial network
    sp_net <- sample_slice_g@spatial_network$cell$spatial_network
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
    cell_metadata$cell_id <- rownames(cell_metadata)
    
    sp_graph_simple_both <- merge(sp_graph_simple_both,
                                  cell_metadata[,c("cell_id","Cell.Type")],
                                  by.x="from",
                                  by.y="cell_id")
    
    sp_graph_simple_both <- merge(sp_graph_simple_both,
                                  cell_metadata[,c("cell_id","Cell.Type")],
                                  by.x="to",
                                  by.y="cell_id",
                                  suffixes=c("_from","_to"))
    
    sp_graph_simple_both$sample_id <- id
    sp_graph_simple_both$slice <- slice
    
    return(sp_graph_simple_both)
    
  })
  
  return(bind_rows(slice_cell_connections))
  
})
cell_connections <- bind_rows(cell_connections)

cell_connect_totals <- ddply(cell_connections,
                             .(from, sample_id, slice),
                             summarise,
                             cell_count=length(unique(to)))


# calc means and medians
cell_connect_stats <- ddply(cell_connect_totals,
                            .(sample_id, slice),
                            summarise,
                            mean=mean(cell_count),
                            median=median(cell_count),
                            sd=sd(cell_count))

cell_connect_stats$stat <- paste0("Mean: ", round(cell_connect_stats$mean, digits=2),
                                  "\nMedian: ", round(cell_connect_stats$median, digits=2),
                                  "\nStdDev: ", round(cell_connect_stats$sd, digits=2))

ggplot(cell_connect_totals,
       aes(cell_count)) +
  geom_histogram(bins=31, color="white", fill="darkblue") +
  theme_bw() +
  labs(x="Number of Connected Cells",
       y="Frequency") +
  facet_grid(sample_id ~ slice) +
  geom_label(data=cell_connect_stats, aes(label=stat), x=32, y=750, size=2.5) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35,40))
ggsave(paste0(out_dir, "sample_cell_connectivity.sliced.histogram.png"),
       width=8, height=7)


# make heatmaps

# hard code cell type order across samples (from 3D analysis)
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

# heatmap per sample, per slice
for (slice in slices) {
  for (id in sample_ids) {
    
    subset_connections <- cell_connections[cell_connections$sample_id == id &
                                             cell_connections$slice == slice,]
    
    # sum up cell 2 cell edge counts
    cell_type_connect_total <- ddply(subset_connections,
                                     .(Cell.Type_from,Cell.Type_to),
                                     summarise,
                                     cell_count=length(from))  
    
    # log scale and z-score
    cell_type_connect_total$scaled_count <- scale(log(cell_type_connect_total$cell_count))
    
    connect_matrix <- acast(cell_type_connect_total,
                            Cell.Type_from ~ Cell.Type_to,
                            value.var = "scaled_count")
    
    pdf(paste0(out_dir, id, ".", slice, ".cell_type_adj_matrix.pdf"), width=8, height=7)
    print(Heatmap(connect_matrix,
                  col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
                  name = "Normalized\nEdge Count",
                  column_title=paste(id,slice),
                  cluster_rows = F, cluster_columns = F,
                  row_order = cell_type_order,
                  column_order = cell_type_order))
    dev.off()
    
    
  }
}

# heatmap per sample
for (id in sample_ids) {
  
  subset_connections <- cell_connections[cell_connections$sample_id == id,]
  
  # sum up cell 2 cell edge counts
  cell_type_connect_total <- ddply(subset_connections,
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

# condition heatmaps
for (slice in slices) {
  for (cond in unique(cell_metadata$Condition)) {
  
  subset_connections <- cell_connections[cell_connections$sample_id %in% 
                                           unique(cell_metadata[cell_metadata$Condition == cond,]$ID) &
                                           cell_connections$slice == slice,]
  
  # sum up cell 2 cell edge counts
  cell_type_connect_total <- ddply(subset_connections,
                                   .(Cell.Type_from,Cell.Type_to),
                                   summarise,
                                   cell_count=length(from))  
  
  # log scale and z-score
  cell_type_connect_total$scaled_count <- scale(log(cell_type_connect_total$cell_count))
  
  connect_matrix <- acast(cell_type_connect_total,
                          Cell.Type_from ~ Cell.Type_to,
                          value.var = "scaled_count")
  
  pdf(paste0(out_dir, cond, ".",slice,".cell_type_adj_matrix.pdf"), width=8, height=7)
  print(Heatmap(connect_matrix,
                col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
                name = "Normalized\nEdge Count",
                column_title=paste(cond,slice),
                cluster_rows = F, cluster_columns = F,
                row_order = cell_type_order,
                column_order = cell_type_order))
  dev.off()
  
  }
}








# 
# 
# 
# 
# 
# 
# # pull out examples to inspect in 3D
# example_cell <- "S4_A1_7359-0"
# 
# cell_connections_example <- cell_connections[cell_connections$from == example_cell,]
# cell_connections_example_meta <- cell_metadata[cell_metadata$cell_id %in% c(cell_connections_example$to,
#                                                                           cell_connections_example$from),]
# 
# example_subset <- subset(m3d_g, 
#                          cell_ids=cell_metadata[cell_metadata$ID == unique(cell_connections_example$sample_id) &
#                                                   cell_metadata$slice == unique(cell_connections_example$slice),]$cell_id)
# 
# example_meta <- cell_metadata[example_subset@cell_metadata$cell$rna$cell_ID,]
# example_meta$example_points <- "other"
# example_meta[example_meta$cell_id == example_cell,]$example_points <- "center"
# example_meta[example_meta$cell_id %in% cell_connections_example$to,]$example_points <- "connect"
# 
# example_subset$example_points <- example_meta$example_points
# 
# spatPlot3D(example_subset,
#            cell_color="example_points",
#            cell_color_code = c("blue","orange","grey"),
#            point_size = 5,
#            save_param = list(save_name="example_3d_plot"))
# 
# spatPlot3D(example_subset,
#            cell_color="is_example",
#            cell_color_code = c("grey","orange","blue"),
#            point_size = 5,
#            axis_scale = "real",
#            save_param = list(save_name="example_3d_plot"))
# 
# example_subset <- createSpatialNetwork(example_subset, 
#                                        method = "Delaunay",
#                                        delaunay_method = "delaunayn_geometry",
#                                        name = "spatial_network")
# 
# 
# # pull out spatial network
# sp_net <- example_subset@spatial_network$cell$spatial_network
# # convert to dataframe
# sp_graph <- as.data.frame(sp_net@networkDT_before_filter)
# 
# # trying out something from Claude
# # Install required packages if needed
# # install.packages(c("geometry", "dplyr"))
# 
# library(geometry)
# library(dplyr)
# 
# 
# 
# # Function to get nearest neighbors using Delaunay triangulation
# get_delaunay_neighbors <- function(df, cell_col = "cell_id", 
#                                    x_col = "x", y_col = "y", z_col = "z") {
#   
#   # Extract coordinates as matrix
#   coords <- as.matrix(df[, c(x_col, y_col, z_col)])
#   
#   # Perform 3D Delaunay triangulation
#   tet <- delaunayn(coords)
#   
#   # Create adjacency list from tetrahedra
#   # Each row in tet represents a tetrahedron with 4 vertices
#   n <- nrow(df)
#   neighbors_list <- vector("list", n)
#   
#   # For each tetrahedron, connect all 4 vertices
#   for (i in 1:nrow(tet)) {
#     vertices <- tet[i, ]
#     
#     # Add each vertex as neighbor to the others
#     for (j in 1:4) {
#       v <- vertices[j]
#       others <- vertices[-j]
#       neighbors_list[[v]] <- c(neighbors_list[[v]], others)
#     }
#   }
#   
#   # Remove duplicates and convert to dataframe
#   result <- lapply(1:n, function(i) {
#     neighs <- unique(neighbors_list[[i]])
#     if (length(neighs) > 0) {
#       data.frame(
#         cell_id = df[[cell_col]][i],
#         neighbor_id = df[[cell_col]][neighs],
#         stringsAsFactors = FALSE
#       )
#     } else {
#       NULL
#     }
#   })
#   
#   result <- do.call(rbind, result)
#   return(result)
# }
# 
# df <- example_meta[,c("cell_id","X", "Y","Z")]
# colnames(df) <- c("cell_id", "x", "y", "z")
# 
# # Get neighbors
# neighbors_df <- get_delaunay_neighbors(df)
# 
# # View results
# head(neighbors_df, 20)
# 
# # Summary statistics
# neighbor_counts <- neighbors_df %>%
#   group_by(cell_id) %>%
#   summarise(n_neighbors = n()) %>%
#   arrange(desc(n_neighbors))
# 
# print(neighbor_counts)
# 
# # Optional: Add distances to neighbors
# neighbors_with_dist <- neighbors_df %>%
#   left_join(df, by = "cell_id") %>%
#   left_join(df, by = c("neighbor_id" = "cell_id"), suffix = c("", "_neighbor")) %>%
#   mutate(
#     distance = sqrt((x - x_neighbor)^2 + (y - y_neighbor)^2 + (z - z_neighbor)^2)
#   ) %>%
#   select(cell_id, neighbor_id, distance) %>%
#   arrange(cell_id, distance)
# 
# # create a full sample subset, see if there is a huge difference in terms of distances
# 
# sample_meta <- cell_metadata[cell_metadata$ID == sample_ids[1],]
# 
# sample_meta_df <- sample_meta[,c("cell_id","X", "Y","Z")]
# colnames(sample_meta_df) <- c("cell_id", "x", "y", "z")
# 
# # Get neighbors
# sample_neighbors_df <- get_delaunay_neighbors(sample_meta_df)
# 
# # Summary statistics
# sample_neighbor_counts <- sample_neighbors_df %>%
#   group_by(cell_id) %>%
#   summarise(n_neighbors = n()) %>%
#   arrange(desc(n_neighbors))
# 
# summary(sample_neighbor_counts)
# 
# 
# 
