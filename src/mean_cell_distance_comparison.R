library(SpatialExperiment)
library(Giotto)
library(plyr)
library(ggplot2)
library(dplyr)
library(plotly)
library(STDistance) # only works for 2D data
library(data.table)
library(FNN)
library(ComplexHeatmap)
library(circlize)

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/mean_cell_distance_comparison/"
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

# calc mean minimum distances between cell types
cell_types <- unique(cell_metadata$Cell.Type)

sample_dists <- lapply(sample_ids, function(id) {
  
  sample_subset <- cell_metadata[cell_metadata$sample_cond == id,]
  
  sample_dt <- as.data.table(sample_subset)
  
  # Create results matrix
  results <- matrix(NA, nrow = length(cell_types), ncol = length(cell_types),
                    dimnames = list(cell_types, cell_types))
  
  # Extract coordinates
  coords <- as.matrix(sample_dt[, .(X, Y, Z)])
  colnames(coords) <- tolower(colnames(coords))
  
  for (i in seq_along(cell_types)) {
    for (j in seq_along(cell_types)) {
      
      
      # Get indices for each cell type
      idx_from <- which(sample_dt$Cell.Type == cell_types[i])
      idx_to <- which(sample_dt$Cell.Type == cell_types[j])
      
      if (i == j) {
        # Within same cell type: use k=2 to get second nearest neighbor
        if (length(idx_from) < 2) {
          results[i, j] <- NA  # Not enough cells to calculate
        } else {
          min_dists <- knnx.dist(coords[idx_from, , drop = FALSE],
                                 coords[idx_from, , drop = FALSE],
                                 k = 2)[, 2]  # Take 2nd column (2nd nearest)
          results[i, j] <- mean(min_dists)
        }
      } else {
        # Between different cell types: use k=1
        min_dists <- knnx.dist(coords[idx_to, , drop = FALSE],
                               coords[idx_from, , drop = FALSE],
                               k = 1)
        results[i, j] <- mean(min_dists)
      }
      
      # Calculate mean of minimum distances
      results[i, j] <- mean(min_dists)
    }
  }
  
  return(results)
})
names(sample_dists) <- sample_ids

# make slices and determine distances

# add slice information to metadata
cell_metadata$slice <- "unused"

cell_metadata[cell_metadata$Z >= 5 &
                cell_metadata$Z <= 15,]$slice <- "Slice 1"
cell_metadata[cell_metadata$Z >= 17 &
                cell_metadata$Z <= 27,]$slice <- "Slice 2"
cell_metadata[cell_metadata$Z >= 29 &
                cell_metadata$Z <= 39,]$slice <- "Slice 3"

slices <- unique(cell_metadata$slice)
slices <- slices[slices != "unused"]

sample_slice_dists_3d <- lapply(sample_ids, function(id) {
  
  slice_dists <- lapply(slices, function(slice) {
    
    sample_subset <- cell_metadata[cell_metadata$sample_cond == id &
                                     cell_metadata$slice == slice,]
    sample_dt <- as.data.table(sample_subset)
    
    # Create results matrix
    results <- matrix(NA, nrow = length(cell_types), ncol = length(cell_types),
                      dimnames = list(cell_types, cell_types))
    
    # Extract coordinates
    coords <- as.matrix(sample_dt[, .(X, Y, Z)])
    colnames(coords) <- tolower(colnames(coords))
    
    for (i in seq_along(cell_types)) {
      for (j in seq_along(cell_types)) {
        
        
        # Get indices for each cell type
        idx_from <- which(sample_dt$Cell.Type == cell_types[i])
        idx_to <- which(sample_dt$Cell.Type == cell_types[j])
        
        if (i == j) {
          # Within same cell type: use k=2 to get second nearest neighbor
          if (length(idx_from) < 2) {
            results[i, j] <- NA  # Not enough cells to calculate
          } else {
            min_dists <- knnx.dist(coords[idx_from, , drop = FALSE],
                                   coords[idx_from, , drop = FALSE],
                                   k = 2)[, 2]  # Take 2nd column (2nd nearest)
            results[i, j] <- mean(min_dists)
          }
        } else {
          # Between different cell types: use k=1
          min_dists <- knnx.dist(coords[idx_to, , drop = FALSE],
                                 coords[idx_from, , drop = FALSE],
                                 k = 1)
          results[i, j] <- mean(min_dists)
        }
        
        # Calculate mean of minimum distances
        results[i, j] <- mean(min_dists)
      }
    }
    
    return(results)
  })

  # determine the mean results
  mean_results <- Reduce(`+`, slice_dists) / length(slice_dists)
  
  return(mean_results)
  
})
names(sample_slice_dists_3d) <- sample_ids

# do it in 2D as well
sample_slice_dists_2d <- lapply(sample_ids, function(id) {
  
  slice_dists <- lapply(slices, function(slice) {
    
    sample_subset <- cell_metadata[cell_metadata$sample_cond == id &
                                     cell_metadata$slice == slice,]
    sample_dt <- as.data.table(sample_subset)
    
    # Create results matrix
    results <- matrix(NA, nrow = length(cell_types), ncol = length(cell_types),
                      dimnames = list(cell_types, cell_types))
    
    # Extract coordinates
    coords <- as.matrix(sample_dt[, .(X, Y)])
    colnames(coords) <- tolower(colnames(coords))
    
    for (i in seq_along(cell_types)) {
      for (j in seq_along(cell_types)) {
        
        
        # Get indices for each cell type
        idx_from <- which(sample_dt$Cell.Type == cell_types[i])
        idx_to <- which(sample_dt$Cell.Type == cell_types[j])
        
        if (i == j) {
          # Within same cell type: use k=2 to get second nearest neighbor
          if (length(idx_from) < 2) {
            results[i, j] <- NA  # Not enough cells to calculate
          } else {
            min_dists <- knnx.dist(coords[idx_from, , drop = FALSE],
                                   coords[idx_from, , drop = FALSE],
                                   k = 2)[, 2]  # Take 2nd column (2nd nearest)
            results[i, j] <- mean(min_dists)
          }
        } else {
          # Between different cell types: use k=1
          min_dists <- knnx.dist(coords[idx_to, , drop = FALSE],
                                 coords[idx_from, , drop = FALSE],
                                 k = 1)
          results[i, j] <- mean(min_dists)
        }
        
        # Calculate mean of minimum distances
        results[i, j] <- mean(min_dists)
      }
    }
    
    return(results)
  })
  
  # determine the mean results
  mean_results <- Reduce(`+`, slice_dists) / length(slice_dists)
  
  return(mean_results)
  
})
names(sample_slice_dists_2d) <- sample_ids

# average them out
dist_mat_3d <- Reduce(`+`, sample_dists) / length(sample_dists)

dist_mat_slice3d <- Reduce(`+`, sample_slice_dists_3d) / length(sample_slice_dists_3d)
dist_mat_slice2d <- Reduce(`+`, sample_slice_dists_2d) / length(sample_slice_dists_2d)

# consistent order
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

Heatmap(dist_mat_3d,
        column_order = cell_type_order,
        row_order = cell_type_order,
        cluster_rows = F,
        cluster_columns = F)

# not sure heatmap is the way to go, going to try rowmeans per cell

total_means <- data.frame(cell_type=rownames(dist_mat_3d),
                          dist_mean=rowMeans(dist_mat_3d),
                          type="50um 3D")
slice2d_means <- data.frame(cell_type=rownames(dist_mat_slice2d),
                            dist_mean=rowMeans(dist_mat_slice2d),
                            type="10um 2D")

all_means <- rbind(total_means,
                   slice2d_means)

total_means <- total_means[order(total_means$dist_mean,
                                 decreasing = T),]
all_means$cell_type <- factor(as.character(all_means$cell_type),
                              levels=total_means$cell_type)
all_means$type <- factor(as.character(all_means$type),
                         levels=c("10um 2D","50um 3D"))

ggplot(all_means,
       aes(x=dist_mean,
       y=cell_type,
       fill=type)) +
  geom_bar(stat="identity", position="dodge",
           color="white") +
  theme_bw() +
  scale_fill_manual(values=c("grey","maroon")) +
  labs(x="Mean Minimum Distance to All Cell Types",
       y=NULL, fill=NULL) +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "mean_min_dist_all_cells.bar_plot.png"), width=7, height=5)

total_diag <- data.frame(cell_type=rownames(dist_mat_3d),
                         dist_mean=diag(dist_mat_3d),
                         type="50um 3D")

slice2d_diag <- data.frame(cell_type=rownames(dist_mat_slice2d),
                            dist_mean=diag(dist_mat_slice2d),
                            type="5um 2D")

all_diag <- rbind(total_diag,
                   slice2d_diag)

total_diag <- total_diag[order(total_diag$dist_mean,
                                 decreasing = T),]
all_diag$cell_type <- factor(as.character(all_diag$cell_type),
                              levels=total_diag$cell_type)
all_diag$type <- factor(as.character(all_diag$type),
                         levels=c("5um 2D","50um 3D"))

ggplot(all_diag,
       aes(x=dist_mean,
           y=cell_type,
           fill=type)) +
  geom_bar(stat="identity", position="dodge",
           color="white") +
  theme_bw() +
  scale_fill_manual(values=c("grey","maroon")) +
  labs(x="Mean Minimum Distance Within Cell Types",
       y=NULL, fill=NULL) +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "mean_min_dist_within_cells.bar_plot.png"), width=7, height=5)

  








