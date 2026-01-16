library(DESeq2)
library(openxlsx)
library(gtools)
library(RColorBrewer)
library(pheatmap)
library(pbayes)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggrepel)

out_dir <- "results/deseq2_workflow/"

dir.create(out_dir, showWarnings = F)

exp_matrix <- read.delim("../WT_v_APP_PS19/rnaseq/3SK7PT_results 1/3SK7PT-expression-matrix.tsv")

metadata <- read.xlsx("../WT_v_APP_PS19/rnaseq/APP-PS19 RNAseq samples list_QW.xlsx", 
                      startRow = 2)

# reformat counts
rownames(exp_matrix) <- exp_matrix$gene_id
any(duplicated(exp_matrix$gene_name))

gene_xref <- unique(exp_matrix[,c("gene_id","gene_name","gene_biotype")])

# if no gene name, just use ensembl id
gene_xref[gene_xref$gene_name == "",]$gene_name <-
  gene_xref[gene_xref$gene_name == "",]$gene_id

# counts
counts_mat <- exp_matrix[,mixedsort(grep("count", colnames(exp_matrix), value=T))]
cpm_mat <- exp_matrix[,mixedsort(grep("cpm", colnames(exp_matrix), value=T))]

# fix column names
colnames(counts_mat) <- paste0("sample", 
                               sapply(colnames(counts_mat), 
                                      function(x) {unlist(strsplit(x, "\\_"))[2]}))

colnames(cpm_mat) <- paste0("sample", 
                               sapply(colnames(cpm_mat), 
                                      function(x) {unlist(strsplit(x, "\\_"))[2]}))

# convert counts to integers
counts_mat <- round(counts_mat)

# fix metadata
metadata$sample_id <- paste0("sample", metadata$`samples.#`)
metadata$condition <- paste0(metadata$Strain, "_",
                        metadata$age)
metadata$condition <- factor(metadata$condition,
                             levels=c("WT_3m",
                                      "WT_9m",
                                      "APPPS19_3m",
                                      "APPPS19_9m",
                                      "APPPS19_12m"))
metadata$sex <- ifelse(grepl("M", metadata$ID), 
                       "M", "F")

# make coldata for DESeq
coldata <- metadata[,c("sample_id", "age", "Strain", "condition","sex")]
rownames(coldata) <- coldata$sample_id



# make DDS object
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = coldata,
                              design = ~0 + condition)

# pre filter counts data (for plotting purposes)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds <- DESeq(dds)

# QC plots
vsd <- vst(dds)

# sample distance matrix
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$sample_id, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

png(paste0(out_dir, "sample_dist.heatmap.png"), 
    width=6, height=4, units="in", res=500)
print(pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors))
dev.off()


plotPCA(vsd, intgroup=c("condition")) +
  scale_color_manual(values=c(brewer.pal(3, "Blues")[2:3],
                              brewer.pal(4, "Oranges")[2:4])) +
  theme_bw()
ggsave(paste0(out_dir, "pca_conditions.png"), width=6, height=4)

plotPCA(vsd, intgroup=c("sex")) +
  theme_bw()
ggsave(paste0(out_dir, "pca_sex.png"), width=6, height=4)

# contrast list
contrasts <- list("APPPS19_3m-WT_3m"=c("condition","APPPS19_3m","WT_3m"),
                  "APPPS19_9m-WT_9m"=c("condition","APPPS19_9m","WT_9m"),
                  "WT_9m-WT_3m"=c("condition","WT_9m","WT_3m"),
                  "APPPS19_9m-APPPS19_3m"=c("condition","APPPS19_9m","APPPS19_3m"),
                  "APPPS19_12m-APPPS19_9m"=c("condition","APPPS19_12m","APPPS19_9m"),
                  "APPPS19_12m-APPPS19_3m"=c("condition","APPPS19_12m","APPPS19_3m"))

results_list <- lapply(contrasts, function(contrast) {
  
  # calculate results
  res <- results(dds, contrast) 
  
  res <- as.data.frame(res)
  # remove NAs
  res <- res[!is.na(res$padj),]
  
  # posterior probability
  res$post_p <- pbayes(res$pvalue, n_cores=2)$posterior_prob
  
  # merge in gene names
  res$gene_id <- rownames(res)
  
  res <- merge(gene_xref,
               res, by="gene_id")
  
  # sort the data
  res <- res[order(res$pvalue),]
  
  res$contrast <- paste0(contrast[2], "-",
                         contrast[3])
  
  return(res)
  
})

# save the results
saveRDS(results_list, file=paste0(out_dir, "deg_results.list.RDS"))

results_df <- bind_rows(results_list)

contrast_order <- c("APPPS19_3m-WT_3m",
                    "APPPS19_9m-WT_9m",
                    "WT_9m-WT_3m",
                    "APPPS19_9m-APPPS19_3m",
                    "APPPS19_12m-APPPS19_9m",
                    "APPPS19_12m-APPPS19_3m")
results_df$contrast <- factor(results_df$contrast,
                              levels=contrast_order)

results_df$log_p <- -log(results_df$pvalue)

sig_results_df <- results_df[results_df$padj < 0.05 &
                               abs(results_df$log2FoldChange) > 0.5,]
table(sig_results_df$contrast)

# volcano plots

logp_threshs <- data.frame(contrast=unique(sig_results_df$contrast),
                           thresh=sapply(unique(unique(sig_results_df$contrast)), function(contrast) {
                             -log(max(sig_results_df[sig_results_df$contrast == contrast,]$pvalue))
                           }))
logp_threshs$contrast <- factor(logp_threshs$contrast,
                                levels=contrast_order)


group <- c("APPPS19_3m-WT_3m",
           "APPPS19_9m-WT_9m")



ggplot(results_df[results_df$contrast %in% group,],
       aes(x=log2FoldChange,
           y=log_p)) +
  geom_point(alpha=0.4, color="black") +
  facet_wrap(~ contrast, scales="free", nrow=1) +
  geom_hline(data=logp_threshs[logp_threshs$contrast %in% group,],
             aes(yintercept = thresh),
             color="red", linetype=2) +
  geom_point(data=sig_results_df[sig_results_df$contrast %in% group,],
             color="red") +
  theme_bw() + labs(x="Log2 Fold Change", y="-log(P-Value)")
ggsave(paste0(out_dir, "appps19_m_wt.volcano_plots.png"), width=7, height=4)


group <- c("WT_9m-WT_3m")

ggplot(results_df[results_df$contrast %in% group,],
       aes(x=log2FoldChange,
           y=log_p)) +
  geom_point(alpha=0.4, color="black") +
  facet_wrap(~ contrast, scales="free", nrow=1) +
  geom_hline(data=logp_threshs[logp_threshs$contrast %in% group,],
             aes(yintercept = thresh),
             color="red", linetype=2) +
  geom_point(data=sig_results_df[sig_results_df$contrast %in% group,],
             color="red") +
  theme_bw() + labs(x="Log2 Fold Change", y="-log(P-Value)")
ggsave(paste0(out_dir, "wt_v_age.volcano_plots.png"), width=4, height=4)

group <- c("APPPS19_9m-APPPS19_3m",
           "APPPS19_12m-APPPS19_9m",
           "APPPS19_12m-APPPS19_3m")

ggplot(results_df[results_df$contrast %in% group,],
       aes(x=log2FoldChange,
           y=log_p)) +
  geom_point(alpha=0.4, color="black") +
  facet_wrap(~ contrast, scales="free", nrow=1) +
  geom_hline(data=logp_threshs[logp_threshs$contrast %in% group,],
             aes(yintercept = thresh),
             color="red", linetype=2) +
  geom_point(data=sig_results_df[sig_results_df$contrast %in% group,],
             color="red") +
  theme_bw() + labs(x="Log2 Fold Change", y="-log(P-Value)")
ggsave(paste0(out_dir, "appps19_v_age.volcano_plots.png"), width=10, height=4)

# counts

ggplot(sig_results_df,
       aes(y=contrast)) +
  geom_bar(color="black", fill="grey") +
  theme_bw() +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(breaks=seq(0, 350, 50)) +
  labs(y="Contrast", x="DEGs at FDR < 0.05")
ggsave(paste0(out_dir, "deg_counts.png"), width=6, height=4)

# pull out top contrast results, plot expression

for (contrast in names(contrasts)) {
  
  contrast_results <- results_list[[contrast]]
  
  condition1 <- unlist(strsplit(contrast, "\\-"))[1]
  condition2 <- unlist(strsplit(contrast, "\\-"))[2]
  
  samples_meta <- metadata[metadata$condition %in% c(condition1,condition2),]

  samples_vsd <- assay(vsd)[,samples_meta$sample_id]  
  
  samples_vsd_long <- melt(samples_vsd)
  colnames(samples_vsd_long) <- c("gene_id","sample_id","value")
  
  # merge in metadata 
  samples_vsd_long <- merge(samples_vsd_long,
                            samples_meta,
                            by="sample_id")
  
  samples_vsd_long <- merge(samples_vsd_long,
                            contrast_results,
                            by="gene_id")
  samples_vsd_long$print_gene <- paste0(samples_vsd_long$gene_name,
                                        "\nlog2FC:  ",round(samples_vsd_long$log2FoldChange, digits=2),
                                        "\nFDR: ", formatC(samples_vsd_long$padj, format = "e", digits=2))
  
  top_genes <- contrast_results$gene_id[0:9]
  
  # refactor for plotting
  samples_vsd_long$condition <- factor(as.character(samples_vsd_long$condition), 
                                       levels=c(condition1,condition2))
  
  ggplot(samples_vsd_long[samples_vsd_long$gene_id %in% top_genes,],
         aes(x=condition,
             y=value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width=0.1) +
    facet_wrap(~ print_gene, ncol=3, scales="free_y") +
    theme_bw() +
    labs(x=NULL, y="VST Expression", title=contrast)
  ggsave(paste0(out_dir, "top_genes.", contrast, ".vst_box_plot.png"), width=8, height=6)
  
  
  
}









