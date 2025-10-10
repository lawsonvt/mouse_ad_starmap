library(SpatialExperiment)
library(Giotto)
library(plyr)
library(ggplot2)
library(ggsci)
library(dplyr)
library(plotly)
library(cowplot)

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/plaque_stats/"
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

# add plaque ID to
plaque_cents <- lapply(sample_ids, function(id) {
  
  data <- plaque_cents[[id]]
  data$sample_id <- id
  return(data)
  
})
plaque_cents_df <- bind_rows(plaque_cents)

# overall cell counts
sample_cell_counts <- as.data.frame(table(cell_metadata$sample_cond))
colnames(sample_cell_counts) <- c("sample_id","count")

ggplot(sample_cell_counts,
       aes(x=reorder(sample_id,count),
           y=count,
           fill=sample_id)) +
  geom_bar(stat="identity") +
  scale_fill_nejm() +
  theme_bw() +
  guides(fill="none") +
  labs(x=NULL, y="Total Cell Count")
ggsave(paste0(out_dir, "sample_cell_counts.png"), width=5, height=4)

# remove small plaques
plaque_cents_df_f <- plaque_cents_df[plaque_cents_df$total_um > 20,]

table(plaque_cents_df_f$sample_id)



plaque_cents_df_f$sample_id <- factor(as.character(plaque_cents_df_f$sample_id),
                                      levels=sample_order)

plaque_counts <- as.data.frame(table(plaque_cents_df_f$sample_id))
colnames(plaque_counts) <- c("sample_id","count")

plaque_counts$sample_id <- factor(as.character(plaque_counts$sample_id),
                                  levels=sample_order)

p2 <- ggplot(plaque_counts,
       aes(x=reorder(sample_id,count),
           y=count,
           fill=sample_id)) +
  geom_bar(stat="identity") +
  scale_fill_nejm() +
  theme_bw() +
  guides(fill="none") +
  labs(x=NULL, y="Plaque Count\nPlaque Volume > 20um^3")

# unfiltered
plaque_counts_all <- as.data.frame(table(plaque_cents_df$sample_id))
colnames(plaque_counts_all) <- c("sample_id","count")

plaque_counts_all$sample_id <- factor(as.character(plaque_counts_all$sample_id),
                                  levels=sample_order)

p1 <- ggplot(plaque_counts_all,
       aes(x=reorder(sample_id,count),
           y=count,
           fill=sample_id)) +
  geom_bar(stat="identity") +
  scale_fill_nejm() +
  theme_bw() +
  guides(fill="none") +
  labs(x=NULL, y="Plaque Count")

plot_grid(p1,p2)
ggsave(paste0(out_dir, "plaque_counts_per_sample.png"), width=10, height=5)

ggplot(plaque_cents_df_f,
       aes(x=total_um,
           color=sample_id)) +
  scale_color_nejm() +
  scale_x_log10() +
  geom_density() +
  theme_bw() +
  labs(x="Plague Volume (in um^3)", y="Density", color=NULL)
ggsave(paste0(out_dir, "plaque_size_per_sample.density.png"), width=6, height=4)

ggplot(plaque_cents_df_f,
       aes(x=total_um,
           fill=sample_id)) +
  facet_wrap(~ sample_id, ncol=2) +
  scale_fill_nejm() +
  scale_x_log10(breaks=c(25,50, 100,250,500,1000,2500,5000, 10000)) +
  geom_histogram( color="white", binwidth = 0.1) +
  guides(fill="none") +
  theme_bw() +
  labs(x="Plague Volume (in um^3)", y="Plaque Count")
ggsave(paste0(out_dir, "plaque_size_per_sample.histograms.png"), width=9, height=7)




  

