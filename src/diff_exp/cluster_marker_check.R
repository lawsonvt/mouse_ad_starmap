library(Giotto)
library(SpatialExperiment)
library(ggplot2)
library(Seurat)
library(plyr)
library(dplyr)
library(cowplot)
library(openxlsx)

# cell type markers (from Stellaromics)
cell_markers <- list(Astrocyte=c("Gfap","Gja1"),
                     "Choroid plexus"=c("Igf2","Ppp1r1b","Phactr2"),
                     E_Baiap3_Tmem163=c("Baiap3", "Tmem163"),
                     E_Cck_Tcf4=c("Cck", "Tcf4"),
                     E_Homer2_Pou3f1=c("Homer2", "Pou3f1"),
                     E_Pou3f1_Cck=c("Pou3f1", "Cck"),
                     E_Rora_Zic1=c("Rora", "Zic1"),
                     E_Tbr1_Id2=c("Tbr1", "Id2"),
                     E_Tshz2_Lamp5=c("Tshz2", "Lamp5"),
                     E_Zic1_Lef1=c("Zic1", "Lef1"),
                     E_Zic1_Syt9=c("Zic1", "Syt9"),
                     Endothelia=c("Cldn5", "Rgs5"),
                     Ependymal=c("Gja1","Hdc","Spef2"),
                     Hybrid_EI_neurons=c("Slc17a7", "Gad1"),
                     I_Chrm2_Pvalb=c("Chrm2", "Pvalb"),
                     I_Npy_Lhx6=c("Npy", "Lhx6"),
                     I_Ppp1r1b_Gpr88=c("Ppp1r1b", "Gpr88"),
                     Microglia=c("Epb41l2", "P2ry12", "Selplg"),
                     Oligodendrocyte=c("Mog", "Opalin"),
                     "Oligodendrocyte progenitor cell"=c("Pdgfra", "Ptprz1", "C1ql1"))

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/cluster_marker_check/"
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

# umap plot
DimPlot(m3d_seurat, reduction = "X_umap", group.by = "rna_Cell.Type",
        label=T, cols=cell_colors, label.box=T, label.color = "black", alpha=1, label.size = 3) +
  labs(x="UMAP 1", y="UMAP 2", title=NULL) +
  theme(legend.position = "none")
ggsave(paste0(out_dir, "cell_cluster_umap.png"), width=9, height=8)

p1 <- DimPlot(m3d_seurat, reduction = "X_umap", group.by = "rna_Cluster",
              label=T, label.box=T, label.color = "black", alpha=1, label.size = 2) +
  labs(x="UMAP 1", y="UMAP 2", title="Harmony Clusters") +
  theme(legend.position = "none")

p2 <- DimPlot(m3d_seurat, reduction = "X_umap", group.by = "rna_Cell.Type",
              label=T, cols=cell_colors, label.box=T, label.color = "black", alpha=1, label.size = 2) +
  labs(x="UMAP 1", y="UMAP 2", title="Cell Types") +
  theme(legend.position = "none")

plot_grid(p1,p2)
ggsave(paste0(out_dir, "cluster_vs_cells_umap.png"), width=14, height=6)

DimPlot(m3d_seurat, reduction = "X_umap", group.by = "rna_Cell.Type",
        split.by = "rna_ID",
        label=T, cols=cell_colors, label.box=F, label.color = "black", alpha=1, 
        label.size = 3, raster=F, ncol=2) +
  theme(legend.position = "none")
ggsave(paste0(out_dir, "clusters_per_sample_umap.png"), width=14, height=12)

m3d_seurat@meta.data$rna_Condition <- factor(as.character(m3d_seurat@meta.data$rna_Condition),
                                             levels=c("WT", "APPPS19"))

DimPlot(m3d_seurat, reduction = "X_umap", group.by = "rna_Cell.Type",
        split.by = "rna_Condition",
        label=T, cols=cell_colors, label.box=F, label.color = "black", alpha=1, 
        label.size = 3, raster=F, ncol=2) +
  labs(x="UMAP 1", y="UMAP 2", title="Cell Types")
ggsave(paste0(out_dir, "clusters_per_condition_umap.png"), width=17, height=8)

# cell counts bar plots

sample_cond_levels <- c("C164B_WT","C166A_WT","C158B_APPPS19","C165_APPPS19")

cell_metadata$sample_cond <- factor(as.character(cell_metadata$sample_cond),
                                    levels=sample_cond_levels)

cell_counts <- sort(table(cell_metadata$Cell.Type),
                    decreasing=T)

cell_metadata$Cell.Type <- factor(as.character(cell_metadata$Cell.Type),
                                  levels=names(cell_counts))

ggplot(cell_metadata,
       aes(x=Cell.Type,
           fill=sample_cond)) +
  geom_bar(position="dodge") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="Cell Count", fill="Sample_Condition")
ggsave(paste0(out_dir, "cell_counts_per_sample.png"), width=10, height=6)


# check if any cell markers are missing / typos
assay_genes <- rownames(m3d_seurat)

cell_markers_missing <- lapply(cell_markers, function(markers) {
  
  markers[!markers %in% assay_genes]
  
})

# try doing the one big dot plot for the harmony clusters
cluster_cell_xref <- unique(cell_metadata[,c("Cluster","Cell.Type")])
any(duplicated(cluster_cell_xref$Cluster))
cluster_cell_xref <- cluster_cell_xref[order(cluster_cell_xref$Cluster),]

write.xlsx(cluster_cell_xref, 
           paste0(out_dir, "cluster_cell_xref.xlsx"),
           colWidths="auto")

Idents(m3d_seurat) <- "rna_Cluster"

dotplot_list <- lapply(names(cell_markers), function(cluster_name) {
  
  markers <- cell_markers[[cluster_name]]
  
  DotPlot(m3d_seurat, features = markers) + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none", 
          text = element_text(size=5),
          axis.text = element_text(size=6))
  
})

plot_grid(plotlist = dotplot_list, nrow = 1)
ggsave(paste0(out_dir, "total_cluster_dotplots.png"), width=18, height=7, bg="white")

Idents(m3d_seurat) <- "rna_Cell.Type"

dotplot_list <- lapply(names(cell_markers), function(cluster_name) {
  
  markers <- cell_markers[[cluster_name]]
  
  DotPlot(m3d_seurat, features = markers) + labs(x=NULL, y=NULL, title=cluster_name) +
    theme(legend.position="none", 
          text = element_text(size=5),
          axis.text = element_text(size=6))
  
})

plot_grid(plotlist = dotplot_list, nrow = 4)
ggsave(paste0(out_dir, "total_celltype_dotplots.png"), width=16, height=12, bg="white")



