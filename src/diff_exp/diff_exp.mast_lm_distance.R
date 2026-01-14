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
library(MAST)
library(reshape2)

out_dir <- "results/diff_exp.mast_lm_distance/"
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

# read in minimum cell distances to plaques
cell_dists <- readRDS("results/plaque_proximity_analysis.um/cell_plaque_min_distances.RDS")

bins <- mixedsort(unique(cell_dists$dist_bin))

bins_xref <- data.frame(dist_bin=bins,
                        bin_name=paste0("dist_bin", 0:(length(bins)-1)))
bins_xref$bin_name <- factor(bins_xref$bin_name,
                             levels=bins_xref$bin_name)

cell_dists <- merge(cell_dists, 
                    bins_xref,
                    by="dist_bin")

# filter out outliers
cell_dists <- cell_dists[cell_dists$min_plaque_dist < 400,]


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

# merge in cell distances

# subset seurat object
m3d_seurat_ad <- subset(m3d_seurat, cells = cell_dists$cell_id)

rownames(cell_dists) <- cell_dists$cell_id

cell_dists <- cell_dists[rownames(m3d_seurat_ad@meta.data),]

# add metadata
m3d_seurat_ad <- AddMetaData(m3d_seurat_ad, cell_dists$min_plaque_dist[,1],
                             col.name = "min_plaque_dist")

m3d_seurat_ad <- AddMetaData(m3d_seurat_ad, cell_dists$bin_name,
                             col.name = "bin_name")

metadata_ad <- m3d_seurat_ad@meta.data

# MAST based model (thanks to Claude)

# per cell type
cell_types <- unique(m3d_seurat_ad@meta.data$rna_Cell.Type)

zlm_results <- lapply(cell_types, function(cell_type) {
  
  celltype_seurat <- subset(m3d_seurat_ad, subset = rna_Cell.Type == cell_type)
  
  # Get log-normalized counts
  norm_counts <- GetAssayData(celltype_seurat, layer = "data")
  
  # Get metadata
  metadata <- celltype_seurat@meta.data
  
  # Create SingleCellAssay
  sca <- FromMatrix(
    exprsArray = as.matrix(norm_counts),
    cData = metadata,
    fData = data.frame(primerid = rownames(norm_counts))
  )
  
  cdr <- colSums(norm_counts > 0)  # count genes detected per cell
  colData(sca)$cngeneson <- scale(cdr)  # standardize it
  
  pct_detected <- rowSums(norm_counts > 0) / ncol(norm_counts)
  
  pct_detected_df <- data.frame(primerid=names(pct_detected),
                                pct_detected=pct_detected)
  
  # By including cngeneson in your model, You're controlling for the 
  # fact that cells with different detection rates might systematically differ 
  # in measured expression. This helps isolate the biological effect of 
  # distance-to-plaque from technical variation in gene detection.
  
  # 4. Fit the model
  # MAST uses a hurdle model (logistic for detection, linear for expression level)
  zlm_fit <- zlm(~ min_plaque_dist + rna_ID + cngeneson, 
                 sca = sca,
                 method = "glm",
                 ebayes = TRUE)  # empirical Bayes shrinkage
  
  # 5. Likelihood ratio test for distance effect
  summary_fit <- summary(zlm_fit, doLRT = "min_plaque_dist")
  
  # 6. Extract results
  results <- summary_fit$datatable
  results_dt <- merge(
    results[contrast == "min_plaque_dist" & component == "H", 
            .(primerid, `Pr(>Chisq)`)],  # hurdle P-value (combined)
    results[contrast == "min_plaque_dist" & component == "logFC",
            .(primerid, coef, ci.hi, ci.lo)],  # log fold-change
    by = "primerid"
  )
  
  # Adjust for multiple testing
  results_dt$fdr <- p.adjust(results_dt$`Pr(>Chisq)`, method = "BH")
  
  # merge in pct detected
  results_dt <- merge(results_dt,
                      pct_detected_df,
                      by="primerid")
  
  # Sort by significance
  results_dt <- results_dt[order(results_dt$fdr), ]
  
  
  
  results_dt$cell_type <- cell_type
  
  return(results_dt)
  
})
names(zlm_results) <- cell_types


# merge it all together

zlm_results_df <- bind_rows(zlm_results)

zlm_results_df <- zlm_results_df[order(zlm_results_df$coef, decreasing=T),]


# plot function

plot_expression <- function(seurat_obj, gene) {


  expr <- GetAssayData(seurat_obj, slot = "data")[gene, ]
  distance <- seurat_obj@meta.data$min_plaque_dist
  sample <- seurat_obj@meta.data$rna_ID
  
  plot_df <- data.frame(
    expression = expr,
    distance = distance,
    sample = sample,
    detected = expr > 0
  )
  
  plot_df <- plot_df[sample(nrow(plot_df)),]
  
  cell_gene <- paste0(gene, " in ", unique(seurat_obj@meta.data$rna_Cell.Type))
  
  ggplot(plot_df,
         aes(x=distance,
             y=expression,
             color=sample)) +
    geom_point(aes(alpha = detected), size = 0.8) +
    scale_alpha_manual(values = c(0.2, 0.8), name = "Detected") +
    theme_bw() +
    # geom_smooth(method = "loess", color = "black", se = TRUE, linewidth = 1) +
    # Overlay linear fit to compare
    geom_smooth(method = "glm", 
                method.args = list(family = gaussian()),
                color = "black", linetype = "dashed", se = FALSE, linewidth = 0.8) +
    labs(title=cell_gene, x="Minimum Distance from Plaque", y="Log Normalized Gene Expression") +
    annotate("text", x = min(plot_df$distance), y = -0.2, 
            label = paste0(round(mean(plot_df$expression == 0) * 100, 1), 
                           "% zeros"),
            hjust = 0, size = 3, color = "gray30") +
    theme(legend.position = "none") +
    scale_color_manual(values=c("darkblue","goldenrod"))
  
}
# find the biggest changes

sig_zlm_results_df <- zlm_results_df[zlm_results_df$fdr < 0.05 &
                                       zlm_results_df$pct_detected > 0.3,]

top_decreasing <- sig_zlm_results_df[order(sig_zlm_results_df$coef,
                                    decreasing=F),][1:12,]

# get the data

dec_plot_list <- lapply(1:nrow(top_decreasing), function(idx) {
  
  data <- top_decreasing[idx,]
  
  seurat_subset <- subset(m3d_seurat_ad, subset = rna_Cell.Type == data$cell_type)
  
  plot_expression(seurat_subset, data$primerid)
  
})

plot_grid(plotlist=dec_plot_list, nrow=3)
ggsave(paste0(out_dir, "top_decreasing_genes.scatter_plot.png"), width=15, height=10)

top_increasing <- sig_zlm_results_df[order(sig_zlm_results_df$coef,
                                           decreasing=T),][1:12,]

# get the data

inc_plot_list <- lapply(1:nrow(top_increasing), function(idx) {
  
  data <- top_increasing[idx,]
  
  seurat_subset <- subset(m3d_seurat_ad, subset = rna_Cell.Type == data$cell_type)
  
  plot_expression(seurat_subset, data$primerid)
  
})

plot_grid(plotlist=inc_plot_list, nrow=3)
ggsave(paste0(out_dir, "top_increasing_genes.scatter_plot.png"), width=15, height=10)

# bin plots?

# plot function

bin_plot_expression <- function(seurat_obj, gene, bin_size) {
  
  
  expr <- GetAssayData(seurat_obj, slot = "data")[gene, ]
  distance <- seurat_obj@meta.data$min_plaque_dist
  sample <- seurat_obj@meta.data$rna_ID
  
  plot_df <- data.frame(
    expression = expr,
    distance = distance,
    sample = sample,
    detected = expr > 0
  )
  
  # create distance bins
  plot_df$distance_bin <- floor(plot_df$distance / bin_size)
  
  plot_df$distance_bin_name <- sapply(plot_df$distance_bin, function(bin) {
    
    paste0(bin * bin_size, " - ", (bin * bin_size) + bin_size-1)
    
  })
  # make correct order for plotting
  plot_df <- plot_df[order(plot_df$distance),]
  
  plot_df$distance_bin_name <- factor(plot_df$distance_bin_name,
                                      levels=unique(plot_df$distance_bin_name))
  
  
  plot_df <- plot_df[sample(nrow(plot_df)),]
  
  plot_df$sample_bin <- paste0(plot_df$sample, plot_df$distance_bin)
  
  cell_gene <- paste0(gene, " in ", unique(seurat_obj@meta.data$rna_Cell.Type))
  
  ggplot(plot_df,
         aes(x=distance_bin_name,
             y=expression)) +
    geom_boxplot(aes(group=distance_bin)) +
    geom_jitter(aes(alpha = detected, color=sample), size = 0.8, width = 0.2) +
    scale_alpha_manual(values = c(0.2, 0.8), name = "Detected") +
    theme_bw() +
  
    labs(title=cell_gene, x="Minimum Distance from Plaque", y="Log Normalized Gene Expression") +
    annotate("text", x = 1, y = -0.2, 
             label = paste0(round(mean(plot_df$expression == 0) * 100, 1), 
                            "% zeros"),
             hjust = 0, size = 3, color = "gray30") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle=35, hjust=1)) +
    scale_color_manual(values=c("darkblue","goldenrod"))
  
}


bin_size <- 50

dec_bin_plot_list <- lapply(1:nrow(top_decreasing), function(idx) {
  
  data <- top_decreasing[idx,]
  
  seurat_subset <- subset(m3d_seurat_ad, subset = rna_Cell.Type == data$cell_type)
  
  bin_plot_expression(seurat_subset, data$primerid, bin_size)
  
})

plot_grid(plotlist=dec_bin_plot_list, nrow=3)
ggsave(paste0(out_dir, "top_decreasing_genes.bin_plot.png"), width=15, height=10)

inc_bin_plot_list <- lapply(1:nrow(top_increasing), function(idx) {
  
  data <- top_increasing[idx,]
  
  seurat_subset <- subset(m3d_seurat_ad, subset = rna_Cell.Type == data$cell_type)
  
  bin_plot_expression(seurat_subset, data$primerid, bin_size)
  
})

plot_grid(plotlist=inc_bin_plot_list, nrow=3)
ggsave(paste0(out_dir, "top_increasing_genes.bin_plot.png"), width=15, height=10)


