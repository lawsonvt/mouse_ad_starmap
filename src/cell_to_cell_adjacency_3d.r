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





