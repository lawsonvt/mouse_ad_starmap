library(SpatialExperiment)
library(Giotto)
library(plyr)
library(ggplot2)
library(ggsci)
library(dplyr)
library(plotly)

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/plaque_proximity_analysis.basic_centroid_volume/"
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

# add sample id to centroids and filter plaques based on size and location
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


sample_min_dists <- lapply(sample_ids, function(id) {
  
  print(id)
  
  cells <- cell_metadata[cell_metadata$sample_cond == id,]
  plaques <- plaque_cents[[id]]
  
  cells_min_dist <- lapply(cells$cell_id, function(cell_id) {
    
    #print(cell_id)
    
    cell_data <- cells[cells$cell_id == cell_id,]
    
    cell_pos <- as.matrix(cell_data[,c("X","Y","Z")])
    
    plaque_dists <- lapply(plaques$plaque_id, function(plaque_id) {
      
      plaque_data <- plaques[plaques$plaque_id == plaque_id,]
      
      plaque_pos <- as.matrix(plaque_data[,c("x_um","y_um","z_um")])
      
      center_dist <- sqrt(sum((cell_pos - plaque_pos)^2))
      surface_dist <- center_dist - plaque_data$radius
      
      data.frame(cell_id=cell_id,
                  plaque_id=plaque_id,
                 center_dist=center_dist,
                 surface_dist=surface_dist)
    })
    plaque_dists <- bind_rows(plaque_dists)
    
    plaque_dists <- plaque_dists[order(plaque_dists$surface_dist),]
    
    return(plaque_dists[1,])
  })
  
  cells_min_dist <- bind_rows(cells_min_dist)
  
  cells_min_dist$sample_id <- id
  
  return(cells_min_dist)
})
sample_min_dists_df <- bind_rows(sample_min_dists)

# save to file
saveRDS(sample_min_dists, file=paste0(out_dir, "cell2plaque_min_distance.RDS"))

sample_min_dists_df$sample_id <- factor(as.character(sample_min_dists_df$sample_id),
                                        levels=sample_order)

ggplot(sample_min_dists_df,
       aes(x=surface_dist,
          fill=sample_id)) +
  geom_histogram( binwidth = 5) +
  facet_wrap(~ sample_id, ncol=2) +
  scale_fill_nejm() +
  theme_bw()+
  guides(fill="none") +
  labs(x="Distance to Plaque", y="Cell Counts")
ggsave(paste0(out_dir, "sample_plaque_distances.bar.png"), width=8, height=5)

ggplot(sample_min_dists_df,
       aes(x=surface_dist,
           color=sample_id)) +
  geom_density(linewidth = 2) +
  scale_color_nejm() +
  scale_fill_nejm() +
  theme_bw() +
  labs(x="Distance to Plaque", y="Density of Cells", color=NULL) +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "sample_plaque_distances.density.png"), width=6, height=5)

# merge in cell type
sample_min_dists_df <- merge(sample_min_dists_df,
                             cell_metadata[,c("cell_id","Cell.Type")],
                             by="cell_id")


# for normalization, what are the cell type counts per sample
sample_celltype_counts <- lapply(sample_ids, function(id) {
  
  sample_meta <- cell_metadata[cell_metadata$sample_cond == id,]
  
  celltype_counts <- as.data.frame(table(sample_meta$Cell.Type),
                                   stringsAsFactors=F)
  colnames(celltype_counts) <- c("Cell.Type", "count")
  
  celltype_counts$total_frac <- celltype_counts$count / sum(celltype_counts$count)
  
  celltype_counts$sample_id <- id
  
  return(celltype_counts)
  
})
sample_celltype_counts <- bind_rows(sample_celltype_counts)

# plot out cell type percentages across all samples

sample_celltype_counts <- sample_celltype_counts[order(sample_celltype_counts$total_frac,
                                                       decreasing=T),]

sample_celltype_counts$Cell.Type <- factor(as.character(sample_celltype_counts$Cell.Type),
                                           levels=unique(sample_celltype_counts$Cell.Type))
sample_celltype_counts$sample_id <- factor(as.character(sample_celltype_counts$sample_id),
                                           levels=sample_order)

ggplot(sample_celltype_counts,
       aes(x=Cell.Type,
           y=total_frac*100,
           fill=sample_id)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_nejm() +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1),
        legend.position = "bottom") +
  labs(x=NULL, y="Cell Type Percentage", fill=NULL)
ggsave(paste0(out_dir, "celltype_percentages_total.bar.png"), width=8, height=5)
 
# just focus on bin 1 for now

bin1_cells <- sample_min_dists_df[sample_min_dists_df$surface_dist < 25,]

bin1_celltype_counts <- lapply(sample_ids, function(id) {
  
  sample_bin1 <- bin1_cells[bin1_cells$sample_id == id,]
  
  celltype_counts <- as.data.frame(table(sample_bin1$Cell.Type),
                                   stringsAsFactors=F)
  colnames(celltype_counts) <- c("Cell.Type", "count")
  
  celltype_counts$bin_frac <- celltype_counts$count / sum(celltype_counts$count)
  
  celltype_counts$sample_id <- id
  
  celltype_counts <- merge(celltype_counts,
                           sample_celltype_counts[,c("sample_id","Cell.Type","total_frac")],
                           by=c("sample_id","Cell.Type"))
  celltype_counts$delta_frac <- celltype_counts$bin_frac - celltype_counts$total_frac
  
  celltype_counts <- celltype_counts[order(celltype_counts$delta_frac,
                                           decreasing=T),]
  
  return(celltype_counts)
})
names(bin1_celltype_counts) <- sample_ids

bin1_celltype_counts_df <- bind_rows(bin1_celltype_counts)

for (id in sample_ids) {
  
  data <- bin1_celltype_counts[[id]]
  
  data$Cell.Type <- factor(as.character(data$Cell.Type),
                           levels=rev(data$Cell.Type))
  
  ggplot(data,
         aes(y=Cell.Type,
             x=delta_frac*100)) +
    geom_bar(stat="identity", color="black", fill="grey") +
    theme_bw() +
    labs(x="Percentage Increase in Prevalence\nWithin 25um of Plaque", y=NULL, title = id)
  ggsave(paste0(out_dir, id, ".percent_increase_barplot.close_plaque.png"), width=6, height=5)
  
}

# bin1 cell counts
bin1_sample_counts <- as.data.frame(table(bin1_cells$sample_id),
                                    stringsAsFactors=F)
colnames(bin1_sample_counts) <- c("sample_id","count")
bin1_sample_counts$type <- "Cells within 25um of Plaque"

plaque_sample_counts <- as.data.frame(table(plaque_cents_df$sample_id),
                                      stringsAsFactors=F)
colnames(plaque_sample_counts) <- c("sample_id", "count")
plaque_sample_counts$type <- "Total Plaques"

bin1_sample_counts <- bin1_sample_counts[order(bin1_sample_counts$count,
                                               decreasing=T),]

total_counts <- rbind(plaque_sample_counts,
                      bin1_sample_counts)

total_counts$sample_x <- factor(as.character(total_counts$sample_id),
                                levels=bin1_sample_counts$sample_id)

total_counts$sample_color <- factor(as.character(total_counts$sample_id),
                                levels=sample_order)
total_counts$type <- factor(as.character(total_counts$type),
                            levels=c("Total Plaques","Cells within 25um of Plaque"))

ggplot(total_counts,
       aes(x=sample_x,
           y=count,
           fill=sample_color)) +
  geom_bar(stat="identity") +
  facet_wrap(~ type, nrow=2, scale="free_y") +
  scale_fill_nejm() +
  theme_bw() +
  guides(fill="none") +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y=NULL)
ggsave(paste0(out_dir, "plaque_counts_v_close_cells.bar.png"), width=4, height=5)

# binning!

distance_bins <- list(c(start=-Inf, stop=25),
                      c(start=25, stop=75),
                      c(start=75, stop=125),
                      c(start=125, stop=175),
                      c(start=175, stop=225),
                      c(start=225, stop=Inf))
names(distance_bins) <- c("< 25um",
                          "25um - 75um",
                          "75um - 125um",
                          "125um - 175um",
                          "175um - 225um",
                          ">= 225um")

ggplot(sample_min_dists_df,
       aes(x=surface_dist,
           fill=sample_id)) +
  geom_histogram( binwidth = 5) +
  facet_wrap(~ sample_id, ncol=2) +
  scale_fill_nejm() +
  geom_vline(xintercept = 25, color="black", linetype=2) +
  geom_vline(xintercept = 75, color="black", linetype=2) +
  geom_vline(xintercept = 125, color="black", linetype=2) +
  geom_vline(xintercept = 175, color="black", linetype=2) +
  geom_vline(xintercept = 225, color="black", linetype=2) +
  theme_bw() +
  guides(fill="none") +
  labs(x="Distance to Plaque", y="Cell Counts")
ggsave(paste0(out_dir, "sample_plaque_distances.bar_bin_lines.png"), width=8, height=5)

ggplot(sample_min_dists_df,
       aes(x=surface_dist,
           color=sample_id)) +
  geom_density(linewidth = 2) +
  scale_color_nejm() +
  scale_fill_nejm() +
  geom_vline(xintercept = 25, color="black", linetype=2) +
  geom_vline(xintercept = 75, color="black", linetype=2) +
  geom_vline(xintercept = 125, color="black", linetype=2) +
  geom_vline(xintercept = 175, color="black", linetype=2) +
  geom_vline(xintercept = 225, color="black", linetype=2) +
  theme_bw() +
  labs(x="Distance to Plaque", y="Density of Cells", color=NULL) +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "sample_plaque_distances.density_bin_lines.png"), width=6, height=5)

# put the cells in these bins
sample_min_dists_df$dist_bin <- ""

for (bin_name in names(distance_bins)) {
  
  start <- distance_bins[[bin_name]][1]
  stop <- distance_bins[[bin_name]][2]
  
  sample_min_dists_df[sample_min_dists_df$surface_dist >= start &
                        sample_min_dists_df$surface_dist < stop,]$dist_bin <- bin_name
  
}

# cell fractions across bins

cell_frac_bins <- lapply(names(distance_bins), function(bin_name) {
  
  sample_fracs <- lapply(sample_ids, function(id) {
    
    dists <- sample_min_dists_df[sample_min_dists_df$dist_bin == bin_name &
                                   sample_min_dists_df$sample_id == id,]
    
    celltype_counts <- as.data.frame(table(dists$Cell.Type),
                                     stringsAsFactors=F)
    colnames(celltype_counts) <- c("Cell.Type", "count")
    
    celltype_counts$bin_frac <- celltype_counts$count / sum(celltype_counts$count)
    
    celltype_counts$sample_id <- id
    
    celltype_counts <- merge(celltype_counts,
                             sample_celltype_counts[,c("sample_id","Cell.Type","total_frac")],
                             by=c("sample_id","Cell.Type"))
    celltype_counts$delta_frac <- celltype_counts$bin_frac - celltype_counts$total_frac
    
    return(celltype_counts)
    
  })
  
  sample_fracs <- bind_rows(sample_fracs)
  
  sample_fracs$dist_bin <- bin_name
  
  return(sample_fracs)
})
cell_frac_bins <- bind_rows(cell_frac_bins)

cell_frac_bins$sample_id <- factor(as.character(cell_frac_bins$sample_id),
                                   levels=sample_order)
cell_frac_bins$dist_bin <- factor(as.character(cell_frac_bins$dist_bin),
                                  levels=names(distance_bins))

cell_frac_bins <- cell_frac_bins[order(cell_frac_bins$delta_frac),]

ggplot(cell_frac_bins,
       aes(y=Cell.Type,
           x=delta_frac*100,
           fill=dist_bin)) +
  geom_jitter(pch=21, height=0.2, size=3) +
  theme_bw() +
  facet_wrap(~ sample_id, ncol=2) +
  geom_vline(xintercept = 0, color="blue") +
  scale_fill_brewer(palette="Reds", direction=-1) +
  labs(x="Cell Percentage Increase", y=NULL, fill="Distance\nTo Plaque")
ggsave(paste0(out_dir, "bin_cell_percentage_increase.jitter.png"),
       width=10, height=7)

# plot with just the diseased

ggplot(cell_frac_bins[cell_frac_bins$sample_id %in%
                        c("C158B_APPPS19","C165_APPPS19"),],
       aes(y=Cell.Type,
           x=delta_frac*100,
           fill=dist_bin)) +
  geom_jitter(pch=21, height=0.2, size=3) +
  theme_bw() +
  facet_wrap(~ sample_id, ncol=2) +
  geom_vline(xintercept = 0, color="blue") +
  scale_fill_brewer(palette="Reds", direction=-1) +
  labs(x="Cell Percentage Increase", y=NULL, fill="Distance\nTo Plaque")
ggsave(paste0(out_dir, "bin_cell_percentage_increase.just_disease.jitter.png"),
       width=9, height=5)

# new plot, histogram of cell distances for specific cell types

ggplot(sample_min_dists_df[sample_min_dists_df$Cell.Type == "Microglia",],
       aes(x=surface_dist,
           fill=sample_id)) +
  geom_histogram( binwidth = 5) +
  facet_wrap(~ sample_id, ncol=2) +
  scale_fill_nejm() +
  theme_bw()+
  guides(fill="none") +
  labs(x="Distance to Plaque", y="Cell Counts", title="Microglia Cells")
ggsave(paste0(out_dir, "sample_plaque_distances.microglia.bar.png"), width=8, height=5)

ggplot(sample_min_dists_df[sample_min_dists_df$Cell.Type == "Microglia",],
       aes(x=surface_dist,
           color=sample_id)) +
  geom_density(linewidth = 2) +
  scale_color_nejm() +
  scale_fill_nejm() +
  theme_bw() +
  labs(x="Distance to Plaque", y="Density of Cells", color=NULL,
       title="Microglia Cells") +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "sample_plaque_distances.microglia.density.png"), width=6, height=5)


ggplot(sample_min_dists_df[sample_min_dists_df$Cell.Type == "Choroid plexus",],
       aes(x=surface_dist,
           fill=sample_id)) +
  geom_histogram( binwidth = 5) +
  facet_wrap(~ sample_id, ncol=2) +
  scale_fill_nejm() +
  theme_bw()+
  guides(fill="none") +
  labs(x="Distance to Plaque", y="Cell Counts", title="Choroid plexus Cells")
ggsave(paste0(out_dir, "sample_plaque_distances.choroid_plexus.bar.png"), width=8, height=5)


ggplot(sample_min_dists_df[sample_min_dists_df$Cell.Type == "Choroid plexus",],
       aes(x=surface_dist,
           color=sample_id)) +
  geom_density(linewidth = 2) +
  scale_color_nejm() +
  scale_fill_nejm() +
  theme_bw() +
  labs(x="Distance to Plaque", y="Density of Cells", color=NULL,
       title="Choroid plexus Cells") +
  theme(legend.position = "bottom")
ggsave(paste0(out_dir, "sample_plaque_distances.choroid_plexus.density.png"), width=6, height=5)




