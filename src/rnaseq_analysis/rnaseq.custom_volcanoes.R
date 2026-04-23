library(openxlsx)
library(ggplot2)
library(ggrepel)
library(snakecase)

out_dir <- "results/rnaseq.custom_volcanoes/"
dir.create(out_dir, showWarnings = F)

in_file <- "data/pathway_annotations/deg_results_highlight.xlsx"

sheet_names <- getSheetNames(in_file)

data_list <- lapply(sheet_names, function(sheet_name) {
  
  return(read.xlsx(in_file, sheet = sheet_name))
  
})
names(data_list) <- sheet_names

# make volcano plots
for (subset in data_list) {
  
 
  # filter low counts
  subset <- subset[subset$baseMean > 10,]
  
  # cap fold changes at 5
  if (any(subset$log2FoldChange > 5)) {
    subset[subset$log2FoldChange > 5,]$log2FoldChange <- 5
  }
  
  if (any(subset$log2FoldChange < -5)) {
    subset[subset$log2FoldChange < -5,]$log2FoldChange <- -5
  }
  
  subset$log_p <- -log10(subset$pvalue)
  
  subset_sig <- subset[subset$padj < 0.1 &
                         abs(subset$log2FoldChange) > 0.5,]
  
  subset_annot <- subset[!is.na(subset$annotate),]
  
  logp_thresh <- min(subset_sig$log_p)
  
  contrast <- unique(subset$contrast)
  
  p1 <- ggplot(subset,
               aes(x=log2FoldChange,
                   y=log_p)) +
    geom_point(alpha=0.3, color="grey60") +
    geom_hline(yintercept = logp_thresh,
               color="red", linetype=2) +
    geom_vline(xintercept = 0.5,
               color="red", linetype=2) +
    geom_vline(xintercept = -0.5,
               color="red", linetype=2) +
    geom_point(data=subset_sig,
               color="red", alpha=0.5) +
    geom_text_repel(data=subset_annot,
                    aes(label=gene_name),
                    color="black", size=2.5,
                    max.overlaps = 50) +
    theme_bw() +
    labs(x="Log2 Fold Change", y="-log10(P-Value)", title=contrast)
  p1
  ggsave(paste0(out_dir, to_snake_case(contrast), ".volcano_plot.png"), width=5, height=5)
  
}



