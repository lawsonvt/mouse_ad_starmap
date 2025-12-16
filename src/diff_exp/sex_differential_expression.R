library(Giotto)
library(SpatialExperiment)
library(ggplot2)
library(cowplot)
library(Seurat)
library(plyr)
library(dplyr)

# Ensure Giotto can access a python env
genv_exists <- checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to
  # install a default Giotto environment
  installGiottoEnvironment()
}

out_dir <- "results/sex_differential_expression/"
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

# split 
seurat_metadata <- m3d_seurat@meta.data

conditions <- unique(seurat_metadata$rna_Condition)
cell_types <- unique(seurat_metadata$rna_Cell.Type)

diff_exp_results <- lapply(conditions, function(condition) {
  
  print(condition)
  
  lapply(cell_types, function(cell_type) {
    
    print(cell_type)
    
    seurat_subset <- subset(m3d_seurat, subset = rna_Condition == condition &
                       rna_Cell.Type == cell_type)
    subset_metadata <- seurat_subset@meta.data
    
    sample_ids <- sort(unique(subset_metadata$rna_ID))
    Idents(seurat_subset) <- "rna_ID"
    
    degs <- FindMarkers(seurat_subset,
                        ident.1=sample_ids[1],
                        ident.2=sample_ids[2],
                        test.use = "MAST")
    degs$gene <- rownames(degs)
    
    degs$condition <- condition
    degs$cell_type <- cell_type
    degs$comparison <- paste0(sample_ids, collapse=" - ")
    
    return(degs)
    
  })
  
})
diff_exp_results <- bind_rows(diff_exp_results)

table(diff_exp_results$cell_type,
      diff_exp_results$comparison)

hist(diff_exp_results$p_val_ad)

sig_diff_exp_results <- diff_exp_results[diff_exp_results$p_val_adj < 0.05,]

table(sig_diff_exp_results$cell_type,
      sig_diff_exp_results$condition)

sig_counts_stats <- lapply(cell_types, function(cell_type) {
  
  sig_cell_subset <- sig_diff_exp_results[sig_diff_exp_results$cell_type == cell_type,]
  
  conditions <- sort(levels(sig_cell_subset$condition))
  
  subset_cond1 <- sig_cell_subset[sig_cell_subset$condition == conditions[1],]
  subset_cond2 <- sig_cell_subset[sig_cell_subset$condition == conditions[2],]
  
  jaccard <- sum(subset_cond1$gene %in% subset_cond2$gene) / 
    (length(subset_cond1$gene) + length(subset_cond2$gene) -
    sum(subset_cond1$gene %in% subset_cond2$gene))
  
  results <- data.frame(cell_type=cell_type,
                        cond1=nrow(subset_cond1),
                        cond2=nrow(subset_cond2),
                        jaccard=jaccard)
  
  colnames(results)[2:3] <- as.character(conditions)
  
  return(results)
  
})
sig_counts_stats <- bind_rows(sig_counts_stats)

ggplot(sig_diff_exp_results,
       aes(y=cell_type,
           fill=condition)) +
  geom_bar(position="dodge") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x="DEGs at FDR < 0.05", y=NULL, title="Sex Differences in Gene Expression",
       fill="Mouse Model") + scale_fill_manual(values=c("maroon","orange"))
ggsave(paste0(out_dir, "sex_degs.bar_plot.png"), width=7, height=5)

# spot checks!

Idents(m3d_seurat) <- "rna_ID"

r1 <- RidgePlot(subset(m3d_seurat, subset = rna_Cell.Type == "E_Tbr1_Id2" &
                   rna_Condition == "APPPS19"), 
          features="Gfap") + theme(legend.position = "none") +
  labs(y=NULL)
v1 <- VlnPlot(subset(m3d_seurat, subset = rna_Cell.Type == "E_Tbr1_Id2" &
                   rna_Condition == "APPPS19"), 
          features="Gfap") + theme(legend.position = "none") +
  labs(x=NULL)

plot_grid(r1, v1)
ggsave(paste0(out_dir, "sig_example.gfap_app_e_tbr1_id2.png"), width=8, height=5)

r2 <- RidgePlot(subset(m3d_seurat, subset = rna_Cell.Type == "Choroid plexus" &
                   rna_Condition == "WT"), 
          features="A2m") + theme(legend.position = "none") +
  labs(y=NULL)
v2 <- VlnPlot(subset(m3d_seurat, subset = rna_Cell.Type == "Choroid plexus" &
                   rna_Condition == "WT"), 
          features="A2m") + theme(legend.position = "none") +
  labs(x=NULL)

plot_grid(r2, v2)
ggsave(paste0(out_dir, "sig_example.a2m_wt_choroid_plexus.png"), width=8, height=5)


r3 <- RidgePlot(subset(m3d_seurat, subset = rna_Cell.Type == "E_Baiap3_Tmem163" &
                   rna_Condition == "WT"), 
          features="Baiap3") + theme(legend.position = "none") +
  labs(y=NULL)
v3 <- VlnPlot(subset(m3d_seurat, subset = rna_Cell.Type == "E_Baiap3_Tmem163" &
                   rna_Condition == "WT"), 
          features="Baiap3") + theme(legend.position = "none") +
  labs(x=NULL)

plot_grid(r3, v3)
ggsave(paste0(out_dir, "sig_example.baiap3_wt_E_Baiap3_Tmem163.png"), width=8, height=5)


