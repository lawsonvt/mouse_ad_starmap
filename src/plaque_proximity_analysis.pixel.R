library(SpatialExperiment)
library(plyr)
library(dplyr)
library(FNN)
library(ggplot2)
library(ggsci)
library(gtools)
library(scales)

out_dir <- "results/plaque_proximity_analysis/"
dir.create(out_dir, recursive = T, showWarnings = F)

# read in spatial experiment converted file
m3d_se <- readRDS("../WT_v_APP_PS19/results/anndata_objects/all_cells_integrated.spatial_exp.RDS")
# pull out cell metadata
cell_metadata <- as.data.frame(colData(m3d_se))
cell_metadata$cell_id <- rownames(cell_metadata)

cell_metadata$simple_cell_id <- as.numeric(sapply(cell_metadata$cell_id, function(x) {
  
  unlist(strsplit(x, "[-_]"))[3]
  
}))

cell_metadata$sample_cond <- paste0(cell_metadata$ID,
                                    "_", cell_metadata$Condition)

sample_order <- c("C164B_WT","C158B_APPPS19",
                  "C166A_WT", "C165_APPPS19")


# get raw location data for each sample (to extract pixel values)
raw_samples <- list(C165_APPPS19=read.csv("../WT_v_APP_PS19/data/C165_APPPS19/C165_APPPS19_cell_metadata.csv"),
                    C158B_APPPS19=read.csv("../WT_v_APP_PS19/data/C158B_APPPS19/C158B_APPPS19_cell_metadata.csv"))

# clean up files
raw_samples <- lapply(raw_samples, function(data) {
  
  data$X <- NULL
  
  return(data)
})

# read in new plaque centroids
plaque_cents <- list(C165_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_boundary_csvs/well4_c165_appps19_plaque_centroids.csv"),
                     C158B_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_boundary_csvs/well11_c158b_appps19_plaque_centroids.csv"))

sample_ids <- names(raw_samples)

# add sample ID and make unique plaque IDs
plaque_cents <- lapply(sample_ids, function(id) {
  
  data <- plaque_cents[[id]]
  data$sample_id <- id
  data$plaque_id <- paste0(data$sample_id,
                           "|", data$plaque_id)
  return(data)
  
})
names(plaque_cents) <- sample_ids


plaque_cents_df <- bind_rows(plaque_cents)

table(plaque_cents_df$sample_id)

# pull in boundary data
plaque_bounds <- list(C165_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_boundary_csvs/well4_c165_appps19_plaque_boundary_pixels.csv"),
                      C158B_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_boundary_csvs/well11_c158b_appps19_plaque_boundary_pixels.csv"))

# add sample ID and make unique plaque IDs
plaque_bounds <- lapply(sample_ids, function(id) {
  
  data <- plaque_bounds[[id]]
  data$sample_id <- id
  data$plaque_id <- paste0(data$sample_id,
                           "|", data$Label.ID)
  return(data)
  
})
names(plaque_bounds) <- sample_ids

plaque_bounds_df <- bind_rows(plaque_bounds)

# determine average plaque spacing
# Analyze spacing for all plaques
analyze_all_plaque_spacing <- function(plaque_boundaries_df) {
  
  results <- data.frame(
    plaque_id = integer(),
    mean_spacing = numeric(),
    median_spacing = numeric(),
    n_points = integer()
  )
  
  for (id in unique(plaque_boundaries_df$plaque_id)) {
    plaque <- as.matrix(plaque_boundaries_df[plaque_boundaries_df$plaque_id == id,][,c("X","Y","Z")])
    
    if (nrow(plaque) > 1) {
      nn_dist <- knn.dist(plaque, k = 1)
      
      results <- rbind(results, data.frame(
        plaque_id = id,
        mean_spacing = mean(nn_dist),
        median_spacing = median(nn_dist),
        n_points = nrow(plaque)
      ))
    }
  }
  
  return(results)
}

spacing_stats <- analyze_all_plaque_spacing(plaque_bounds_df)

# it seems we can get away with just modeling the bounds as a series of points instead of a 3D surface

# using the data "as is", i.e. no filtering based on size or location

all_boundary_points <- plaque_bounds_df[,c("X","Y","Z")]

sample_min_dists <- lapply(sample_ids, function(id) {
  
  pixel_cell_coords <- raw_samples[[id]]
  sample_cell_metadata <- cell_metadata[cell_metadata$sample_cond == id,]
 
  
  sample_cell_metadata <- merge(sample_cell_metadata,
                                pixel_cell_coords,
                                by.x="simple_cell_id",
                                by.y="Cell_id")
  
  # filter down to used cells
  cell_coords <- sample_cell_metadata
  rownames(cell_coords) <- cell_coords$cell_id
  cell_coords <- as.matrix(cell_coords[,c("X_pixels","Y_pixels","Z_pixels")])
  colnames(cell_coords) <- c("X","Y","Z")
  
  # Find nearest neighbor for each cell to any boundary point
  nn_result <- FNN::knnx.dist(
    data = all_boundary_points,  # reference points (plaque boundaries)
    query = cell_coords,          # query points (cells)
    k = 1                         # find 1 nearest neighbor
  )
  
  sample_cell_metadata$min_plaque_dist <- nn_result
  
  return(sample_cell_metadata)
  
})

sample_min_dists_df <- bind_rows(sample_min_dists)

# make some plots

ggplot(sample_min_dists_df,
       aes(x=min_plaque_dist,
           fill=sample_cond)) +
  geom_histogram( binwidth = 5) +
  facet_wrap(~ sample_cond, ncol=1) +
  scale_fill_manual(values=c("#0072B5FF","#20854EFF")) +
  theme_bw()+
  guides(fill="none") +
  labs(x="Pixel Distance to Plaque", y="Cell Counts")
ggsave(paste0(out_dir, "plaque_dist_bar_plot.png"), width=5, height=5)


ggplot(sample_min_dists_df,
       aes(x=min_plaque_dist,
           color=sample_cond)) +
  geom_density(linewidth = 1) +
  scale_color_manual(values=c("#0072B5FF","#20854EFF")) +
  theme_bw() +
  labs(x="Pixel Distance to Plaque", y="Density of Cells", color=NULL) +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "plaque_dist_density_plot.total.png"), width=5, height=4)

ggplot(sample_min_dists_df,
       aes(x=min_plaque_dist,
           color=sample_cond)) +
  geom_density(linewidth = 1) +
  scale_color_manual(values=c("#0072B5FF","#20854EFF")) +
  #scale_fill_nejm() +
  facet_wrap(~ Cell.Type, ncol=5) +
#  scale_x_log10 () +
  theme_bw() +
  labs(x="Pixel Distance to Plaque", y="Density of Cells", color=NULL) +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "plaque_dist_density_plot.per_celltype.png"), width=10, height=7)

table(sample_min_dists_df[sample_min_dists_df$min_plaque_dist < 100,]$sample_cond)

# for normalization, what are the cell type counts per sample
sample_celltype_counts <- lapply(sample_ids, function(id) {
  
  sample_meta <- cell_metadata[cell_metadata$sample_cond == id,]
  
  celltype_counts <- as.data.frame(table(sample_meta$Cell.Type),
                                   stringsAsFactors=F)
  colnames(celltype_counts) <- c("Cell.Type", "count")
  
  celltype_counts$total_frac <- celltype_counts$count / sum(celltype_counts$count)
  
  celltype_counts$sample_cond <- id
  
  return(celltype_counts)
  
})
sample_celltype_counts <- bind_rows(sample_celltype_counts)

bin_size <- 50
bins <- c(seq(0, 500, by=bin_size), Inf)

sample_min_dists_df$dist_bin <- ""

# bin em
for (i in 1:(length(bins)-1)) {
  
  print(i)
  
  start <- bins[i]
  stop <- bins[i+1]
  
  bin_name <- paste0(start, " - ", stop)
  
  sample_min_dists_df[sample_min_dists_df$min_plaque_dist >= start &
                        sample_min_dists_df$min_plaque_dist < stop,]$dist_bin <- bin_name
  
}

# count up bins per cell, per sample
sample_dist_bin_cell_counts <- ddply(sample_min_dists_df,
                                .(sample_cond, Cell.Type, dist_bin),
                                summarise,
                                bin_celltype_count=length(cell_id))

sample_dist_bin_counts <- ddply(sample_min_dists_df,
                                     .(sample_cond, dist_bin),
                                     summarise,
                                     bin_total_count=length(cell_id))
sample_dist_bin_cell_counts <- merge(sample_dist_bin_cell_counts,
                                     sample_dist_bin_counts,
                                     by=c("sample_cond", "dist_bin"))
sample_dist_bin_cell_counts$bin_frac <- sample_dist_bin_cell_counts$bin_celltype_count / 
  sample_dist_bin_cell_counts$bin_total_count

# merge in total fractions
sample_dist_bin_cell_counts <- merge(sample_dist_bin_cell_counts,
                                     sample_celltype_counts,
                                     by=c("sample_cond","Cell.Type"))


sample_dist_bin_cell_counts$delta_frac <- sample_dist_bin_cell_counts$bin_frac - 
  sample_dist_bin_cell_counts$total_frac

sample_dist_bin_cell_counts$dist_bin <- factor(as.character(sample_dist_bin_cell_counts$dist_bin),
                                               levels=mixedsort(unique(as.character(sample_dist_bin_cell_counts$dist_bin))))
sample_dist_bin_cell_counts$sample_cond <- factor(as.character(sample_dist_bin_cell_counts$sample_cond),
                                                  levels=sample_order)

ggplot(sample_dist_bin_cell_counts,
       aes(x=dist_bin,
           y=delta_frac*100,
           fill=sample_cond)) +
  geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~ Cell.Type, ncol=5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1),
        legend.position = "bottom") +
  scale_fill_manual(values=c("#0072B5FF","#20854EFF")) +
  labs(x="Pixel Distance", y="Cell Percentage Increase", fill=NULL)
ggsave(paste0(out_dir, "plaque_dist.cell_percent_increase.png"), width=12, height=8)

# save the distances
saveRDS(sample_min_dists_df, file=paste0(out_dir, "cell_plaque_min_distances.RDS"))



