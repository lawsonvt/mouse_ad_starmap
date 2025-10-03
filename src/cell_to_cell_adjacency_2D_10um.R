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

out_dir <- "results/cell_to_cell_adjacency_2D_control/"
dir.create(out_dir, recursive = T, showWarnings = F)


# read in spatial experiment converted file

m3d_se <- readRDS("../WT_v_APP_PS19/results/anndata_objects/all_cells_integrated.spatial_exp.RDS")
# pull out cell metadata
cell_metadata <- as.data.frame(colData(m3d_se))
cell_metadata$cell_id <- rownames(cell_metadata)

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

cell_metadata[cell_metadata$Z >= 5 &
                cell_metadata$Z <= 15,]$slice <- "Slice 1"
cell_metadata[cell_metadata$Z >= 17 &
                cell_metadata$Z <= 27,]$slice <- "Slice 2"
cell_metadata[cell_metadata$Z >= 29 &
                cell_metadata$Z <= 39,]$slice <- "Slice 3"

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

# visualize slices in 3D

# colors for consistent plotting
cell_colors <- metadata(m3d_se)$"Cell Type_colors"

spatPlot3D(m3d_g,
           cell_color="Cell.Type",
           cell_color_code=cell_colors,
           point_size = 2,
           save_plot = F)

# add in slice metadata
m3d_g$slice <- cell_metadata$slice

spatPlot3D(m3d_g,
           cell_color="slice",
           cell_color_code = c(pal_futurama()(3), "grey"),
           point_size = 5,
           save_param = list(save_name="slice_3d_plot"))

spatPlot3D(m3d_g,
           cell_color="slice",
           cell_color_code = c(pal_futurama()(3), "grey"),
           point_size = 2,
           axis_scale = "real",
           save_param = list(save_name="slice_3d_plot.axis_real"))

slices <- paste0("Slice ", 1:3)

# split data into slices

slice_list <- lapply(slices, function(slice) {
  
  subset(m3d_g, cell_ids=cell_metadata[cell_metadata$slice == slice,]$cell_id)
  
})
names(slice_list) <- slices


# cell to cell adjacency per sample
sample_ids <- levels(cell_metadata$ID)
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
  geom_label(data=cell_connect_stats, aes(label=stat), x=35, y=1400, size=3) +
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35,40))
ggsave(paste0(out_dir, "sample_cell_connectivity.sliced.histogram.png"),
       width=8, height=7)







