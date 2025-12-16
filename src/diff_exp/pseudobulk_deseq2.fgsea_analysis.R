library(msigdbr)  # for gene sets / pathways
library(dplyr) # data manipulation
library(fgsea) # GSEA test
library(openxlsx) # excel output
library(snakecase) # better file names
library(ggplot2) # plotting
library(cowplot) # combining plots
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

# Initial setup ----------------------------------------------------------------

out_dir <- "results/pseudobulk_deseq2.fgsea_analysis/"
dir.create(out_dir, recursive = T, showWarnings = F)

# read in previous results
deg_results <- readRDS("results/diff_exp.pseudo_bulk_deseq2/degs.app_ps19_minus_wt.RDS")

# Load in msigdbr and get gene sets --------------------------------------------
reactome_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
wikipath_gene_sets <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:WIKIPATHWAYS")

gobp_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
gomf_gene_sets <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:MF")

# NOTE: More pathways can be added, check out msigdbr_collections() for more info

# function to convert data frames to lists
list_convert <- function(gene_sets) {
  gene_sets %>%
    split(x = .$gene_symbol, f = .$gs_name)
}

# put into one big list
total_gene_sets <- list(reactome=list_convert(reactome_gene_sets),
                        wikipathways=list_convert(wikipath_gene_sets),
                        gobp=list_convert(gobp_gene_sets),
                        gomf=list_convert(gomf_gene_sets))


# Calculate GSEA results per cell type -----------------------------------------

gsea_results <- lapply(deg_results, function(degs) {
  
  print(unique(degs$cell_type))
  
  ranked_genes <- degs %>%
    filter(!is.na(stat) &
             !is.infinite(stat)) %>%
    arrange(desc(stat)) %>%
    pull(stat, name = gene_name)
  
  print("Run GSEA on ranked DEGs ...")
  
  # for each set of gene sets, run GSEA
  total_gsea_results <- lapply(names(total_gene_sets), function(gs_name) {
    
    print(gs_name)
    
    gene_sets <- total_gene_sets[[gs_name]]
    # run GSEA
    fgsea_results <- fgsea(
      pathways = gene_sets,
      stats = ranked_genes,
      minSize = 15, # Minimal size of a gene set to test
      maxSize = 500 # Maximal size of a gene set to test
    )
    
    return(fgsea_results[order(fgsea_results$pval),])
  })
  names(total_gsea_results) <- names(total_gene_sets)
  
  return(total_gsea_results)
})

# save results
saveRDS(gsea_results, file=paste0(out_dir, "gsea_results.RDS"))

# Output results to excel files ------------------------------------------------

# for each cell type, pull out results and write to excel
for (cell_type in names(gsea_results)) {
  
  cell_results <- gsea_results[[cell_type]]
  
  write.xlsx(cell_results,
             file = paste0(out_dir, 
                           to_snake_case(cell_type), 
                           ".gsea_results.xlsx"),
             colWidths="auto")
}

# combine them to see where pathways are actually being affected ---------------

gsea_results_df <- lapply(names(gsea_results), function(cell_type) {
  
  cell_results <- bind_rows(gsea_results[[cell_type]])
  cell_results$cell_type <- cell_type

  return(cell_results)  
})
gsea_results_df <- bind_rows(gsea_results_df)

# make the pathway name prettier for plotting
gsea_results_df$pathway_pretty <- sapply(gsea_results_df$pathway, function(x) {
  
  str_wrap(paste0(unlist(str_split(x, "_"))[-1], collapse=" "), width=40)
           
})


sig_gsea_results_df <- gsea_results_df[gsea_results_df$padj < 0.1,]

sig_pathway_counts <- table(sig_gsea_results_df$pathway)

top_pathways <- names(sig_pathway_counts[sig_pathway_counts > 1])

nes_matrix <- acast(gsea_results_df[gsea_results_df$pathway %in% top_pathways,],
                   cell_type ~ pathway_pretty,
                   value.var="NES")

fdr_matrix <- acast(gsea_results_df[gsea_results_df$pathway %in% top_pathways,],
                    cell_type ~ pathway_pretty,
                    value.var="padj")

fdr_matrix[is.na(fdr_matrix)] <- 1

fdr_matrix_print <- fdr_matrix
fdr_matrix_print[fdr_matrix <= 0.1] <- "*"
fdr_matrix_print[fdr_matrix > 0.1] <- ""

pdf(paste0(out_dir, "top_pathways_heatmap.pdf"), width=10, height=9)
Heatmap(nes_matrix,
        col=colorRamp2(c(-2,0,2), c("blue","white","red")),
        cell_fun = function(j, i, x, y, w, h, col) {
          grid.text(fdr_matrix_print[i,j], x, y)
        }, name="NES", column_names_max_height = unit(10, "cm"),
        column_names_gp = gpar(fontsize=8))
dev.off()













