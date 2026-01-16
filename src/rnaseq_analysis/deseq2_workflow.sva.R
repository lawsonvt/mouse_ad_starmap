library(DESeq2)
library(sva)
library(openxlsx)
library(gtools)
library(RColorBrewer)
library(pheatmap)
library(pbayes)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(limma)

out_dir <- "results/deseq2_workflow.sva/"

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
metadata$sex <- factor(ifelse(grepl("M", metadata$ID), 
                       "M", "F"))

# make coldata for DESeq
coldata <- metadata[,c("sample_id", "age", "Strain", "condition","sex")]
rownames(coldata) <- coldata$sample_id

# make DDS object
dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData = coldata,
                              design = ~condition)


# pre filter counts data (for plotting purposes)
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Get normalized counts for SVA
dds_norm <- estimateSizeFactors(dds)
norm_counts <- counts(dds_norm, normalized = TRUE)

# Apply variance stabilizing transformation for SVA input
dat <- vst(dds_norm, blind = FALSE)
dat_matrix <- assay(dat)

# Create model matrices
# Full model: includes your biological variables of interest
# Example: if comparing condition (e.g., treatment vs control)
mod <- model.matrix(~ sex + condition, colData(dds))

# Null model: intercept only (no biological variables)
mod0 <- model.matrix(~ 1, colData(dds))

# Estimate number of surrogate variables
# This step can take several minutes
n_sv <- num.sv(dat_matrix, mod, method = "be")

# Calculate surrogate variables
svobj <- sva(dat_matrix, mod, mod0, n.sv = n_sv)

# Add SVs to colData
for (i in 1:n_sv) {
  colData(dds)[, paste0("SV", i)] <- svobj$sv[, i]
}

# Create new design formula including SVs
# Build the SV terms dynamically
sv_terms <- paste0("SV", 1:n_sv, collapse = " + ")
design_formula <- as.formula(paste("~", sv_terms, "+ sex + condition"))

# Update the design
design(dds) <- design_formula

dds <- DESeq(dds)

# compare PCA plots!

# PCA before correction
vsd_before <- vst(dds_norm, blind = FALSE)
pcaData_before <- plotPCA(vsd_before, intgroup = "condition", returnData = TRUE)
percentVar_before <- round(100 * attr(pcaData_before, "percentVar"))

p1 <- ggplot(pcaData_before, aes(x = PC1, y = PC2, color = condition, shape=sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_before[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_before[2], "% variance")) +
  ggtitle("PCA - Before SVA Correction") +
  scale_color_manual(values=c(brewer.pal(3, "Blues")[2:3],
                              brewer.pal(4, "Oranges")[2:4])) +
  theme_bw()

# PCA after correction
vsd_corrected <- vst(dds, blind = FALSE)
assay(vsd_corrected) <- limma::removeBatchEffect(
  assay(vsd_corrected),
  covariates = svobj$sv,
  design = mod
)

pcaData_after <- plotPCA(vsd_corrected, intgroup = "condition", returnData = TRUE)
percentVar_after <- round(100 * attr(pcaData_after, "percentVar"))

p2 <- ggplot(pcaData_after, aes(x = PC1, y = PC2, color = condition, shape=sex)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_after[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_after[2], "% variance")) +
  ggtitle("PCA - After SVA Correction") +
  scale_color_manual(values=c(brewer.pal(3, "Blues")[2:3],
                              brewer.pal(4, "Oranges")[2:4])) +
  theme_bw()

plot_grid(p1,p2)
ggsave(paste0(out_dir, "pca.sva_before_and_after.png"), width=10, height=4)

# proceed with the analysis!


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

# save to excel
write.xlsx(results_list, file=paste0(out_dir, "deg_results.xlsx"),
           colWidths="auto")

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
  geom_text_repel(data=sig_results_df[sig_results_df$contrast %in% group,],
                  aes(label=gene_name),
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
  geom_text_repel(data=sig_results_df[sig_results_df$contrast %in% group,],
                  aes(label=gene_name),
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
  geom_text_repel(data=sig_results_df[sig_results_df$contrast %in% group,],
                  aes(label=gene_name),
                  color="red") +
  theme_bw() + labs(x="Log2 Fold Change", y="-log(P-Value)")
ggsave(paste0(out_dir, "appps19_v_age.volcano_plots.png"), width=10, height=4)

# counts

ggplot(sig_results_df,
       aes(y=contrast)) +
  geom_bar(color="black", fill="grey") +
  theme_bw() +
  scale_y_discrete(limits=rev) +
  scale_x_continuous(breaks=seq(0, 350, 25)) +
  labs(y="Contrast", x="DEGs at FDR < 0.05")
ggsave(paste0(out_dir, "deg_counts.png"), width=6, height=4)

# pull out top contrast results, plot expression

for (contrast in names(contrasts)) {
  
  contrast_results <- results_list[[contrast]]
  
  # put in print gene for plotting
  contrast_results$print_gene <- paste0(contrast_results$gene_name,
                                        "\nlog2FC:  ",round(contrast_results$log2FoldChange, digits=2),
                                        "\nFDR: ", formatC(contrast_results$padj, format = "e", digits=2))
  
  contrast_results$print_gene <- factor(contrast_results$print_gene,
                                        levels=unique(contrast_results$print_gene))
  
  condition1 <- unlist(strsplit(contrast, "\\-"))[1]
  condition2 <- unlist(strsplit(contrast, "\\-"))[2]
  
  samples_meta <- metadata[metadata$condition %in% c(condition1,condition2),]
  
  samples_vsd <- assay(vsd_corrected)[,samples_meta$sample_id]  
  
  samples_vsd_long <- melt(samples_vsd)
  colnames(samples_vsd_long) <- c("gene_id","sample_id","value")
  
  # merge in metadata 
  samples_vsd_long <- merge(samples_vsd_long,
                            samples_meta,
                            by="sample_id")
  
  samples_vsd_long <- merge(samples_vsd_long,
                            contrast_results,
                            by="gene_id")
  
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

# example of using posteriors
app_res_9_3 <- results_list[["APPPS19_9m-APPPS19_3m"]]
app_res_12_9 <- results_list[["APPPS19_12m-APPPS19_9m"]]
app_res_12_3 <- results_list[["APPPS19_12m-APPPS19_3m"]]

# add in column suffixes
colnames(app_res_9_3)[4:ncol(app_res_9_3)] <- paste0(colnames(app_res_9_3)[4:ncol(app_res_9_3)],
                                                     ".9_3")
colnames(app_res_12_9)[4:ncol(app_res_12_9)] <- paste0(colnames(app_res_12_9)[4:ncol(app_res_12_9)],
                                                     ".12_9")
colnames(app_res_12_3)[4:ncol(app_res_12_3)] <- paste0(colnames(app_res_12_3)[4:ncol(app_res_12_3)],
                                                       ".12_3")


# merge em up
app_res <- merge(app_res_9_3,
                 app_res_12_9,
                 by=c("gene_id","gene_name","gene_biotype"))
app_res <- merge(app_res,
                 app_res_12_3,
                 by=c("gene_id","gene_name","gene_biotype"))


# joint p for change in both
app_res$joint_p_12_3 <- app_res$post_p.9_3 * app_res$post_p.12_9
app_res$same_dir_12_3 <- sign(app_res$log2FoldChange.9_3) == sign(app_res$log2FoldChange.12_9)


app_res$joint_p_late <- (1-app_res$post_p.9_3) * app_res$post_p.12_9

app_res$joint_p_early <- app_res$post_p.9_3 * (1-app_res$post_p.12_9)


# make some plots
app_metadata <- coldata[coldata$Strain == "APPPS19",]


samples_vsd <- assay(vsd_corrected)[,app_metadata$sample_id]  

samples_vsd_long <- melt(samples_vsd)
colnames(samples_vsd_long) <- c("gene_id","sample_id","value")

# merge in metadata 
samples_vsd_long <- merge(samples_vsd_long,
                          app_metadata,
                          by="sample_id")

samples_vsd_long <- merge(samples_vsd_long,
                          gene_xref,
                          by="gene_id")


late_genes <- app_res[order(app_res$joint_p_late, decreasing=T),c("gene_id","gene_name","joint_p_late"),]

ggplot(samples_vsd_long[samples_vsd_long$gene_id %in% late_genes$gene_id[1:5],],
       aes(x=condition,
           y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1) +
  facet_wrap(~ gene_name, ncol=3, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="VST Expression")
ggsave(paste0(out_dir, "late_change_genes.appps19.boxplot.png"), width=7, height=5)


early_genes <- app_res[order(app_res$joint_p_early, decreasing=T),c("gene_id","gene_name","joint_p_early"),]

ggplot(samples_vsd_long[samples_vsd_long$gene_id %in% early_genes$gene_id[1:9],],
       aes(x=condition,
           y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1) +
  facet_wrap(~ gene_name, ncol=3, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="VST Expression")
ggsave(paste0(out_dir, "early_change_genes.appps19.boxplot.png"), width=7, height=7)


constant <- app_res[order(app_res$joint_p_12_3, decreasing=T),
                             c("gene_id","gene_name","joint_p_12_3","same_dir_12_3"),]


ggplot(samples_vsd_long[samples_vsd_long$gene_id %in% constant$gene_id[1:6],],
       aes(x=condition,
           y=value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1) +
  facet_wrap(~ gene_name, ncol=3, scales="free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(x=NULL, y="VST Expression", title=contrast)
ggsave(paste0(out_dir, "consistent_change_genes.appps19.boxplot.png"), width=7, height=5)


