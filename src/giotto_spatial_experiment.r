library(SpatialExperiment)
library(Giotto)
library(plyr)
library(ggplot2)
library(dplyr)

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
cell_metadata <- as.data.frame(colData(m3d_se))

# convert to giotto object ...

m3d_g <- spatialExperimentToGiotto(m3d_se)
# add in z axis
m3d_g@spatial_locs$cell$raw$sdimz <- spatialCoords(m3d_se)[,"Z"]

instructions(m3d_g) <- instrs

# plotUMAP_3D(m3d_g,
#             cell_color = "Cell.Type", 
#             show_center_label = FALSE,
#             dim_reduction_name = "X_umap",
#             save_param = list(save_name = "umap_leiden.cell_types")
# )

cell_colors <- metadata(m3d_se)$"Cell Type_colors"

spatPlot2D(m3d_g,
           cell_color="Cell.Type",
           cell_color_code=cell_colors,
           point_size = 1,
           save_param = list(save_name = "integrated.2d_splot"))

spatPlot3D(m3d_g,
           cell_color="Cell.Type",
           axis_scale = "real",
           cell_color_code=cell_colors,
           point_size = 2,
           save_param = list(save_name = "integrated.3d_splot_real"))

spatPlot3D(m3d_g,
           cell_color="Cell.Type",
           cell_color_code=cell_colors,
           save_param = list(save_name = "integrated.3d_splot"))


# create 2d plots per sample
sample_ids <- levels(cell_metadata$ID)

m3d_samples <- lapply(sample_ids, function(id) {
  
  cell_ids <- rownames(cell_metadata[cell_metadata$ID == id,])
  
  subset(m3d_g, cell_ids=cell_ids)
  
})
names(m3d_samples) <- sample_ids

# make 2D plots per sample
for (id in sample_ids) {
  
  print(id)
  spatPlot2D(m3d_samples[[id]],
             cell_color="Cell.Type",
             cell_color_code=cell_colors,
             title=id,
             point_size = 1,
             save_param=list(save_name=paste0(id, ".2d_splot")))
  
  spatPlot3D(m3d_samples[[id]],
             cell_color="Cell.Type",
             cell_color_code=cell_colors,
             title=id,
             save_param=list(save_name=paste0(id, ".3d_splot")))
  
  spatPlot3D(m3d_samples[[id]],
             cell_color="Cell.Type",
             cell_color_code=cell_colors,
             axis_scale="real",
             title=id,
             save_param=list(save_name=paste0(id, ".3d_splot_real")))
  
}

# create spatial network

m3d_g <- createSpatialNetwork(m3d_g, method = "Delaunay",
                              delaunay_method = "delaunayn_geometry",
                              name = "spatial_network")

sp_net <- m3d_g@spatial_network$cell$spatial_network

sp_graph <- as.data.frame(sp_net@networkDT)

sp_graph_simple <- sp_graph[,c("from","to","distance","weight")]
sp_graph_simple_rev <- sp_graph_simple
colnames(sp_graph_simple_rev) <- c("to","from","distance","weight")

sp_graph_simple_both <- rbind(sp_graph_simple,
                              sp_graph_simple_rev)

sp_graph_simple_both <- unique(sp_graph_simple_both)


# merge in metadata
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


sp_graph_total <- ddply(sp_graph_simple_both,
                        .(from, Cell.Type_from),
                        summarise,
                        cell_count=length(unique(to)))

ggplot(sp_graph_total,
       aes(cell_count)) +
  geom_histogram(bins=31, color="white", fill="darkblue") +
  theme_bw() +
  labs(x="Number of Connected Cells",
       y="Frequency") +
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30))

# same analysis per sample
cell_connect_totals <- lapply(sample_ids, function(id) {
  
  print(id)
  
  g <- m3d_samples[[id]]
  
  g <- createSpatialNetwork(g, method = "Delaunay",
                       delaunay_method = "delaunayn_geometry",
                       name = "spatial_network")
  
  sp_net <- g@spatial_network$cell$spatial_network
  
  sp_graph <- as.data.frame(sp_net@networkDT)
  
  sp_graph_simple <- sp_graph[,c("from","to","distance","weight")]
  sp_graph_simple_rev <- sp_graph_simple
  colnames(sp_graph_simple_rev) <- c("to","from","distance","weight")
  
  sp_graph_simple_both <- rbind(sp_graph_simple,
                                sp_graph_simple_rev)
  
  sp_graph_simple_both <- unique(sp_graph_simple_both)
  
  
  # merge in metadata
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
  
  
  sp_graph_total <- ddply(sp_graph_simple_both,
                          .(from, Cell.Type_from),
                          summarise,
                          cell_count=length(unique(to)))
  
  sp_graph_total$sample_id <- id
  sp_graph_total$sample_cond <- paste0(id, " (",
                                       unique(cell_metadata[cell_metadata$ID == id,]$Condition), ")")
  
  return(sp_graph_total)
})
cell_connect_totals <- bind_rows(cell_connect_totals)

ggplot(cell_connect_totals,
       aes(cell_count)) +
  geom_histogram(bins=31, color="white", fill="darkblue") +
  theme_bw() +
  labs(x="Number of Connected Cells",
       y="Frequency") +
  facet_wrap(~ sample_cond, ncol=2) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30))
ggsave(paste0(out_dir, "sample_cell_connectivity.histogram.png"), width=10, height=8)


# spatPlot3D(m3d_g,  
#            show_network = TRUE,
#            network_color = "blue", 
#            spatial_network_name = "spatial_network"
# )