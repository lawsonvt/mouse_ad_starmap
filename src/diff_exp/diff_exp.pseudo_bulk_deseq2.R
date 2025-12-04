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

out_dir <- "results/diff_exp.pseudo_bulk_deseq2/"
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

# add in sample xref data
metadata <- m3d_seurat@meta.data

# merge in Xref data
metadata <- merge(metadata,
                  sample_xref,
                  by.x="rna_ID",
                  by.y="sample")

# fix the order
rownames(metadata) <- metadata$cell_ID
metadata <- metadata[rownames(m3d_seurat@meta.data),]

m3d_seurat <- AddMetaData(m3d_seurat, metadata$sex, col.name="sex")
m3d_seurat <- AddMetaData(m3d_seurat, metadata$type, col.name="type")

# aggregate for pseudobulk analysis
m3d_seurat_agg <- AggregateExpression(m3d_seurat,
                                      group.by=c("rna_Cell.Type","sex","type"),
                                      return.seurat = T)

# per cell type
cell_types <- unique(m3d_seurat_agg@meta.data$rna_Cell.Type)

deg_results <- lapply(cell_types, function(cell_type) {
  
  subset_agg <- subset(m3d_seurat_agg, subset = rna_Cell.Type == cell_type)
  
  # start up DESeq2
  dds <- DESeqDataSetFromMatrix(
    countData = subset_agg[["rna"]]$counts,
    colData = subset_agg@meta.data,
    design = ~ type + sex
  )
  
  # run DESeq
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("type", "WT", "APP-PS19"))
  
  res <- as.data.frame(res)
  res <- res[order(res$pvalue),]
  
  res$gene_name <- rownames(res)
  res$cell_type <- cell_type
  
  # remove low count genes
  res <- res[!is.na(res$padj),]
  
  return(res)
})
names(deg_results) <- cell_types

# combine results
deg_results_df <- bind_rows(deg_results)

# write to excel
write.xlsx(deg_results, file=paste0(out_dir, "degs.wt_minus_app_ps19.xlsx"),
           colWidths="auto")

# save for future analyses
saveRDS(deg_results, file=paste0(out_dir, "degs.wt_minus_app_ps19.RDS"))

# look at sig ones
sig_deg_results_df <- deg_results_df[deg_results_df$padj < 0.1,]

cell_counts <- sort(table(sig_deg_results_df$cell_type))
sig_deg_results_df$cell_type <- factor(as.character(sig_deg_results_df$cell_type),
                                       levels=names(cell_counts))

# make some plots

ggplot(sig_deg_results_df,
       aes(y=cell_type)) +
  geom_bar(fill="maroon", color="black") +
  theme_bw() +
  labs(x="DEGs at FDR < 0.1", y=NULL, title="Differential Expression between WT and APP/PS19")
ggsave(paste0(out_dir, "deg_counts_per_celltype.png"), width=8, height=6)


# heatmaps!

sig_gene_counts <- table(sig_deg_results_df$gene_name)

top_genes <- names(sig_gene_counts[sig_gene_counts > 2])

fc_matrix <- acast(deg_results_df[deg_results_df$gene_name %in% top_genes,],
                   cell_type ~ gene_name,
                   value.var="log2FoldChange")

fdr_matrix <- acast(deg_results_df[deg_results_df$gene_name %in% top_genes,],
                   cell_type ~ gene_name,
                   value.var="padj")

fdr_matrix[is.na(fdr_matrix)] <- 1

fdr_matrix_print <- fdr_matrix
fdr_matrix_print[fdr_matrix <= 0.1] <- "*"
fdr_matrix_print[fdr_matrix > 0.1] <- ""

pdf(paste0(out_dir, "top_degs_heatmap.pdf"), width=8, height=6)
Heatmap(fc_matrix,
        col=colorRamp2(c(-2,0,2), c("blue","white","red")),
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(fdr_matrix_print[i,j], x, y)
        }, name="log2FC")
dev.off()


