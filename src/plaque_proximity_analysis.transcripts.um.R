library(SpatialExperiment)
library(plyr)
library(dplyr)
library(FNN)
library(ggplot2)
library(ggsci)
library(gtools)
library(scales)

out_dir <- "results/plaque_proximity_analysis.transcripts.um/"
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

# load in transcript level data
sample_trans_locs <- list(C165_APPPS19=read.csv("../WT_v_APP_PS19/data/C165_APPPS19/C165_APPPS19_cell_assigned_counts.csv"),
                          C158B_APPPS19=read.csv("../WT_v_APP_PS19/data/C158B_APPPS19/C158B_APPPS19_cell_assigned_counts.csv"))


# read in new plaque centroids
plaque_cents <- list(C165_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_boundary_csvs/well4_c165_appps19_plaque_centroids_um.csv"),
                     C158B_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_boundary_csvs/well11_c158b_appps19_plaque_centroids_um.csv"))

sample_ids <- names(plaque_cents)

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

# pull original plaque centroids to determine um sizes
plaque_cents_old <- list(C165_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_csv_files/C165_APPPS19_plaque_centers_unfiltered_shifted_um.csv"),
                         C158B_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_csv_files/C158B_APPPS19_plaque_centers_unfiltered_shifted_um.csv"))

# add sample ID and make unique plaque IDs
plaque_cents_old <- lapply(sample_ids, function(id) {
  
  data <- plaque_cents_old[[id]]
  data$sample_id <- id
  data$plaque_id <- paste0(data$sample_id,
                           "|", data$plaque_id)
  return(data)
  
})
names(plaque_cents_old) <- sample_ids
plaque_cents_old_df <- bind_rows(plaque_cents_old)

min_plaque_size <- min(plaque_cents_old_df[plaque_cents_old_df$total_um >= 20,]$total_pixels)



# pull in boundary data
plaque_bounds <- list(C165_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_boundary_csvs/well4_c165_appps19_plaque_boundaries_um.csv"),
                      C158B_APPPS19=read.csv("../WT_v_APP_PS19/data/plaque_boundary_csvs/well11_c158b_appps19_plaque_boundaries_um.csv"))

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

length(unique(plaque_bounds_df$plaque_id))


# filter down data to lower 30 Z locations

cell_metadata <- cell_metadata[cell_metadata$Z < 30,]
plaque_cents_df <- plaque_cents_df[plaque_cents_df$z_um < 30 &
                                     plaque_cents_df$area > min_plaque_size,]

plaque_bounds_df <- plaque_bounds_df[plaque_bounds_df$plaque_id %in% plaque_cents_df$plaque_id,]

sample_trans_locs <- lapply(sample_trans_locs, function(data) {
  
  return(data[data$Z_global_um < 30,])
  
})

# reformat for knn
all_boundary_points <- plaque_bounds_df[,c("X_um","Y_um","Z_um")]
colnames(all_boundary_points) <- c("X","Y","Z")



sample_min_dists <- lapply(sample_ids, function(id) {
  
  sample_trans_metadata <- sample_trans_locs[[id]]
  
  # filter down to used cells
  trans_coords <- sample_trans_metadata
  trans_coords <- as.matrix(trans_coords[,c("X_global_um","Y_global_um","Z_global_um")])
  colnames(trans_coords) <- c("X","Y","Z")
  
  # Find nearest neighbor for each cell to any boundary point
  nn_result <- FNN::knnx.dist(
    data = all_boundary_points,  # reference points (plaque boundaries)
    query = trans_coords,          # query points (cells)
    k = 1                         # find 1 nearest neighbor
  )
  
  sample_trans_metadata$min_plaque_dist <- nn_result
  sample_trans_metadata$sample_id <- id
  
  return(sample_trans_metadata)
  
})

sample_min_dists_df <- bind_rows(sample_min_dists)

ggplot(sample_min_dists_df,
       aes(x=min_plaque_dist,
           fill=sample_id)) +
  geom_histogram( binwidth = 5) +
  facet_wrap(~ sample_id, ncol=1) +
  scale_fill_manual(values=c("#0072B5FF","#20854EFF")) +
  theme_bw()+
  guides(fill="none") +
  labs(x="μm Distance to Plaque", y="Transcript Counts")
ggsave(paste0(out_dir, "transcript_plaque_dist_bar_plot.png"), width=5, height=5)

ggplot(sample_min_dists_df,
       aes(x=min_plaque_dist,
           color=sample_id)) +
  geom_density(linewidth = 1) +
  scale_color_manual(values=c("#0072B5FF","#20854EFF")) +
  theme_bw() +
  labs(x="μm Distance to Plaque", y="Density of Transcripts", color=NULL) +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "transcript_plaque_dist_density_plot.png"), width=5, height=4)


# for normalization, what are the gene type counts per sample
sample_genetype_counts <- lapply(sample_ids, function(id) {
  
  sample_meta <- sample_trans_locs[[id]]
  
  genetype_counts <- as.data.frame(table(sample_meta$Gene),
                                   stringsAsFactors=F)
  colnames(genetype_counts) <- c("Gene", "count")
  
  genetype_counts$total_frac <- genetype_counts$count / sum(genetype_counts$count)
  
  genetype_counts$sample_id <- id
  
  return(genetype_counts)
  
})
sample_genetype_counts <- bind_rows(sample_genetype_counts)

bin_size <- 25
bins <- c(seq(0, 300, by=bin_size), Inf)

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

# count up bins per gene, per sample
sample_dist_bin_gene_counts <- ddply(sample_min_dists_df,
                                     .(sample_id, Gene, dist_bin),
                                     summarise,
                                     bin_genetype_count=length(X))

sample_dist_bin_counts <- ddply(sample_min_dists_df,
                                .(sample_id, dist_bin),
                                summarise,
                                bin_total_count=length(X))
sample_dist_bin_gene_counts <- merge(sample_dist_bin_gene_counts,
                                     sample_dist_bin_counts,
                                     by=c("sample_id", "dist_bin"))
sample_dist_bin_gene_counts$bin_frac <- sample_dist_bin_gene_counts$bin_genetype_count / 
  sample_dist_bin_gene_counts$bin_total_count

# merge in total fractions
sample_dist_bin_gene_counts <- merge(sample_dist_bin_gene_counts,
                                     sample_genetype_counts,
                                     by=c("sample_id","Gene"))

sample_dist_bin_gene_counts$delta_frac <- sample_dist_bin_gene_counts$bin_frac - 
  sample_dist_bin_gene_counts$total_frac

sample_dist_bin_gene_counts$dist_bin <- factor(as.character(sample_dist_bin_gene_counts$dist_bin),
                                               levels=mixedsort(unique(as.character(sample_dist_bin_gene_counts$dist_bin))))
sample_dist_bin_gene_counts$sample_id <- factor(as.character(sample_dist_bin_gene_counts$sample_id),
                                                  levels=sample_order)


# find genes with increases
inc_genes <- unique(sample_dist_bin_gene_counts[abs(sample_dist_bin_gene_counts$delta_frac) > 0.01,]$Gene)

ggplot(sample_dist_bin_gene_counts[sample_dist_bin_gene_counts$Gene %in% inc_genes,],
       aes(x=dist_bin,
           y=delta_frac*100,
           fill=sample_id)) +
  geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~ Gene, ncol=5, scales="free_y") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1),
        legend.position = "bottom") +
  scale_fill_manual(values=c("#0072B5FF","#20854EFF")) +
  labs(x="μm Distance To Plaque", y="Gene Transcript Percentage Increase", fill=NULL)
ggsave(paste0(out_dir, "plaque_dist.top_gene_percent_increase.png"), width=13, height=8)


endo_markers <- c("Cldn5","Igf2","Rgs5")


ggplot(sample_min_dists_df[sample_min_dists_df$Gene %in% endo_markers,],
       aes(x=min_plaque_dist,
           color=sample_id)) +
  geom_density(linewidth = 1) +
  facet_wrap(~ Gene, ncol=3) +
  scale_color_manual(values=c("#0072B5FF","#20854EFF")) +
  theme_bw() +
  labs(x="μm Distance to Plaque", y="Density of Transcripts", color=NULL) +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "endothelial_markers.plaque_dist_denisty_plot.png"), width=10, height=5)


ggplot(sample_dist_bin_gene_counts[sample_dist_bin_gene_counts$Gene %in% endo_markers,],
       aes(x=dist_bin,
           y=delta_frac*100,
           fill=sample_id)) +
  geom_bar(stat="identity", position = "dodge") +
  facet_wrap(~ Gene, ncol=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=35, hjust=1),
        legend.position = "bottom") +
  scale_fill_manual(values=c("#0072B5FF","#20854EFF")) +
  labs(x="μm Distance To Plaque", y="Gene Transcript Percentage Increase", fill=NULL)
ggsave(paste0(out_dir, "endothelial_markers.plaque_dist.gene_percent_increase.png"), height=8, width=4)






