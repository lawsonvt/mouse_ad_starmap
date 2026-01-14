library(Giotto)
library(SpatialExperiment)
library(ggplot2)
library(cowplot)
library(Seurat)
library(plyr)
library(dplyr)
library(openxlsx)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(gtools)

out_dir <- "results/find_cluster_markers/"
dir.create(out_dir, recursive = T, showWarnings = F)

# start giotto (for conversion)
instrs <- createGiottoInstructions(
  save_dir = out_dir,
  save_plot = TRUE,
  show_plot = FALSE,
  return_plot = FALSE
)

# read in sample xref
sample_xref <- read.xlsx("sample_xref.xlsx")
# make some factors
sample_xref$type <- factor(sample_xref$type,
                           levels=c("WT", "APP_PS19"))
sample_xref$sex <- factor(sample_xref$sex,
                          levels=c("Female", "Male"))

# read in spatial experiment converted file
m3d_se <- readRDS("../WT_v_APP_PS19/results/anndata_objects/all_cells_integrated.spatial_exp.RDS")

# convert to giotto object ...
m3d_g <- spatialExperimentToGiotto(m3d_se)
# add in z axis, since it does not get added
m3d_g@spatial_locs$cell$raw$sdimz <- spatialCoords(m3d_se)[,"Z"]
# get the actual raw data
m3d_g@expression$cell$rna$raw@exprMat <- m3d_g@expression$cell$rna$counts@exprMat

# remove dimensional reduction
m3d_g@dimension_reduction$cells$cell$rna$`2D_spatial` <- NULL

# convert to Seurat?
m3d_seurat <- giottoToSeuratV5(m3d_g)

# prep data for differential expression
m3d_seurat <- NormalizeData(m3d_seurat)
m3d_seurat <- FindVariableFeatures(m3d_seurat)
m3d_seurat <- ScaleData(m3d_seurat)

harmony_clusters <- unique(m3d_seurat@meta.data$rna_Cluster)

cell_xref <- unique(m3d_seurat@meta.data[,c("rna_Cluster","rna_Cell.Type")])

cell_xref <- cell_xref[order(cell_xref$rna_Cluster),]

Idents(m3d_seurat) <- "rna_Cluster"

cluster_degs <- lapply(harmony_clusters, function(cluster) {
  
  print(cluster)
  
  marker_degs <- FindMarkers(m3d_seurat,
                             ident.1 = cluster,
                             test.use = "MAST")
  marker_degs$gene <- rownames(marker_degs)
  marker_degs$cluster <- cluster
  marker_degs$cell_type <- cell_xref[cell_xref$rna_Cluster == cluster,]$rna_Cell.Type
  marker_degs$delta_pct <- marker_degs$pct.1 - marker_degs$pct.2
  
  return(marker_degs)
  
})
names(cluster_degs) <- harmony_clusters

# write to file
write.xlsx(cluster_degs[mixedorder(names(cluster_degs))], file=paste0(out_dir, "cluster_v_all_degs.xlsx"), 
           colWidths="auto")


cluster_degs_df <- bind_rows(cluster_degs)

# try and filter them down
sig_cluster_degs_df <- cluster_degs_df[cluster_degs_df$avg_log2FC > 0 &
                                         cluster_degs_df$p_val_adj < 0.05,]



# top genes for each

top_n <- 5

cluster_top <- lapply(cluster_degs, function(degs) {
  
  degs <- degs[degs$avg_log2FC > 0 &
                 degs$p_val_adj < 0.05,]
  degs <- degs[order(degs$delta_pct, decreasing = T),]
  
  return(head(degs, n=top_n))
  
})
cluster_top_df <- bind_rows(cluster_top)

cluster_top_df <- cluster_top_df[order(cluster_top_df$cluster),]

top_genes <- cluster_top_df[cluster_top_df$avg_log2FC > 2.5,]$gene

DotPlot(m3d_seurat, features=unique(top_genes)) + RotatedAxis() +
  labs(x=NULL, y="Cluster") + theme(panel.grid.major = element_line(color = "grey", linewidth = 0.2, linetype = "solid"))
ggsave(paste0(out_dir, "marker_dotplot.png"), width=18, height=10, bg="white")


# get cluster xref and output to excel



write.xlsx(cell_xref, paste0(out_dir, "cell_xref.xlsx"), colWidths="auto")



