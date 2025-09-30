library(SpatialExperiment)
library(Giotto)
library(plyr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/giotto_spatial_experiment/"
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

# NOTE: THIS DOES NOT WORK
# plotUMAP_3D(m3d_g,
#             cell_color = "Cell.Type", 
#             show_center_label = FALSE,
#             dim_reduction_name = "X_umap",
#             save_param = list(save_name = "umap_leiden.cell_types")
# )

# colors for consistent plotting
cell_colors <- metadata(m3d_se)$"Cell Type_colors"

# 2D plot
spatPlot2D(m3d_g,
           cell_color="Cell.Type",
           cell_color_code=cell_colors,
           point_size = 1,
           save_param = list(save_name = "integrated.2d_splot"))

# 3D plot with real axes
spatPlot3D(m3d_g,
           cell_color="Cell.Type",
           axis_scale = "real",
           cell_color_code=cell_colors,
           point_size = 2,
           save_param = list(save_name = "integrated.3d_splot_real"))

# 3D plot with equal sized axes
spatPlot3D(m3d_g,
           cell_color="Cell.Type",
           cell_color_code=cell_colors,
           save_param = list(save_name = "integrated.3d_splot"))


# create 2d plots per sample
sample_ids <- levels(cell_metadata$ID)

# separate out samples by subsetting
m3d_samples <- lapply(sample_ids, function(id) {
  
  cell_ids <- rownames(cell_metadata[cell_metadata$ID == id,])
  
  subset(m3d_g, cell_ids=cell_ids)
  
})
names(m3d_samples) <- sample_ids

# make 2D plots per sample
for (id in sample_ids) {
  
  print(id)
  # 2D plot
  spatPlot2D(m3d_samples[[id]],
             cell_color="Cell.Type",
             cell_color_code=cell_colors,
             title=id,
             point_size = 1,
             save_param=list(save_name=paste0(id, ".2d_splot")))
  # 3D plot with equal sized axes
  spatPlot3D(m3d_samples[[id]],
             cell_color="Cell.Type",
             cell_color_code=cell_colors,
             title=id,
             save_param=list(save_name=paste0(id, ".3d_splot")))
  # 3D plot with real axes
  spatPlot3D(m3d_samples[[id]],
             cell_color="Cell.Type",
             cell_color_code=cell_colors,
             axis_scale="real",
             title=id,
             save_param=list(save_name=paste0(id, ".3d_splot_real")))
  
}

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
  cell_connect_total$sample_id <- id
  cell_connect_total$sample_cond <- paste0(id, " (",
                                       unique(cell_metadata[cell_metadata$ID == id,]$Condition), ")")

  return(cell_connect_total)
})
cell_connect_totals <- bind_rows(cell_connect_totals)

# calc means and medians
cell_connect_stats <- ddply(cell_connect_totals,
                            .(sample_id),
                            summarise,
                            mean=mean(cell_count),
                            median=median(cell_count),
                            sd=sd(cell_count))

cell_connect_stats$stat <- paste0("Mean: ", round(cell_connect_stats$mean, digits=2),
                                  "\nMedian: ", round(cell_connect_stats$median, digits=2),
                                  "\nStdDev: ", round(cell_connect_stats$sd, digits=2))

# plot histogram of counts across samples
ggplot(cell_connect_totals,
       aes(cell_count)) +
  geom_histogram(bins=31, color="white", fill="darkblue") +
  theme_bw() +
  labs(x="Number of Connected Cells",
       y="Frequency") +
  facet_wrap(~ sample_cond, ncol=2) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30, 35, 40)) +
  geom_label(data=cell_connect_stats, aes(label=stat), x=30, y=3100, size=3)
ggsave(paste0(out_dir, "sample_cell_connectivity.histogram.png"), width=6, height=5)

# create adjacency matrices

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
          column_title=id))
  dev.off()
  
}

# one big adjacency matrix
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

# filter out astrocytes and oligodendrocytes

filter_ao <- c("Astrocyte", "Oligodendrocyte")

# create adjacency matrices

for (id in sample_ids) {
  
  print(id)
  
  cell_connect <- cell_connections[[id]]
  
  cell_connect <- cell_connect[!cell_connect$Cell.Type_from %in% filter_ao &
                                 !cell_connect$Cell.Type_to %in% filter_ao,]
  
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
  
  pdf(paste0(out_dir, id, ".no_astro_oligo.cell_type_adj_matrix.pdf"), width=8, height=7)
  print(Heatmap(connect_matrix,
                col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
                name = "Normalized\nEdge Count",
                column_title=id))
  dev.off()
  
}

# one big adjacency matrix
cell_connections_big <- bind_rows(cell_connections)

cell_connections_big <- cell_connections_big[!cell_connections_big$Cell.Type_from %in% filter_ao &
                               !cell_connections_big$Cell.Type_to %in% filter_ao,]

cell_type_connect_big <- ddply(cell_connections_big,
                               .(Cell.Type_from,Cell.Type_to),
                               summarise,
                               cell_count=length(from))  

cell_type_connect_big$scaled_count <- scale(log(cell_type_connect_big$cell_count))

connect_big_matrix <- acast(cell_type_connect_big,
                            Cell.Type_from ~ Cell.Type_to,
                            value.var = "scaled_count")

pdf(paste0(out_dir, "all_samples.no_astro_oligo.cell_type_adj_matrix.pdf"), width=8, height=7)
print(Heatmap(connect_big_matrix,
              col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
              name = "Normalized\nEdge Count",
              column_title="All Samples"))
dev.off()

