library(openxlsx)
library(ggplot2)
library(snakecase)
library(stringr)

out_dir <- "results/rnaseq.pathway_plots/"

dir.create(out_dir, showWarnings = F)

# function to clean up pathway names
pathway_pretty <- function(p) {

  # drop DB name and replace underscores
  p <- paste0(unlist(strsplit(p, "_"))[-1], collapse=" ")
  
  # make title
  p <- str_to_title(p)
  
  # drop WP IDs
  p <- gsub("Wp[0-9]+", "", p)
  p <- trimws(p)
  
  # wrap text
  p <- str_wrap(p, width=40)
  
  return(p)
  
}

# get excel files

gsea_files <- list.files("data/pathway_annotations/", ".xlsx", full.names = T)

# filter down files
gsea_files <- gsea_files[!grepl("(deg)|\\~", gsea_files)]

# make plots for each file
for (gsea_file in gsea_files) {
  
  data <- read.xlsx(gsea_file, sheet = "annotate")
  
  data$pathway_pretty <- sapply(data$pathway, pathway_pretty)
  data$logp <- -log10(data$pval)
  # ensure proper order
  data <- data[order(data$pval),]
  data$pathway_pretty <- factor(as.character(data$pathway_pretty),
                                levels=rev(as.character(data$pathway_pretty)))

  ggplot(data,
         aes(x=logp, y=pathway_pretty,fill=NES)) +
    geom_point(pch=21, size=5) +
    theme_bw() +
    xlim(0, max(data$logp)) +
    scale_fill_gradient2(
      low  = "blue",
      mid  = "white",
      high = "red",
      midpoint = 0
    ) +
    labs(x="-log10(PValue)", y=NULL)  
  
  out_file <- paste0(out_dir, gsub(".xlsx", ".png", basename(gsea_file)))
  
  ggsave(out_file, width=7, height=6)
  
}