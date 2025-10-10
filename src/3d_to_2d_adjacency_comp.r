library(ComplexHeatmap)
library(circlize)
library(plyr)
library(reshape2)
library(ggplot2)
library(ggsci)

out_dir <- "results/3d_to_2d_adjacency_comp/"

dir.create(out_dir, showWarnings = F)

connect_3d <- readRDS("results/cell_to_cell_adjacency_3d/sample_delaunay_connections.RDS")
connect_2d <- readRDS("results/cell_to_cell_adjacency_2D_10um_flatten/sample_slice_delaunay_connections.RDS")

connect_2d$condition <- sapply(connect_2d$sample_id, function(x) {unlist(strsplit(x, " "))[2]})
connect_3d$condition <- sapply(connect_3d$sample_id, function(x) {unlist(strsplit(x, " "))[2]})

cell_type_order <- c("Astrocyte",
                     "Oligodendrocyte",
                     "E_Homer2_Pou3f1",
                     "I_Npy_Lhx6",
                     "E_Tbr1_Id2",
                     "Endothelia",
                     "E_Rora_Zic1",
                     "I_Chrm2_Pvalb",
                     "Oligodendrocyte progenitor cell",
                     "Microglia",
                     "E_Tshz2_Lamp5",
                     "Hybrid_EI_neurons",
                     "E_Cck_Tcf4",
                     "E_Pou3f1_Cck",
                     "E_Zic1_Lef1",
                     "E_Baiap3_Tmem163",
                     "E_Zic1_Syt9",
                     "Choroid plexus",
                     "Ependymal",
                     "I_Ppp1r1b_Gpr88")

sample_ids <- unique(connect_2d$sample_id)
conditions <- unique(connect_2d$condition)

# make heatmap per sample

for (id in sample_ids) {
  
  sample_3d <- connect_3d[connect_3d$sample_id == id,]
  sample_2d <- connect_2d[connect_2d$sample_id == id,]
  
  total_3d <- ddply(sample_3d,
                                   .(Cell.Type_from,Cell.Type_to),
                                   summarise,
                                   cell_count=length(from)) 
  total_3d$type <- "3D"
  
  total_2d_slice <- ddply(sample_2d,
                          .(Cell.Type_from,Cell.Type_to, slice),
                          summarise,
                          slice_cell_count=length(from)) 
  
  total_2d <- ddply(total_2d_slice, 
                    .(Cell.Type_from,Cell.Type_to),
                    summarise,
                    cell_count=mean(slice_cell_count))
  total_2d$type <- "2D"
  
  total_all <- rbind(total_3d,
                     total_2d)
  
  total_all$scaled_count <- scale(log(total_all$cell_count))
  
  connect_matrix_3d <- acast(total_all[total_all$type == "3D",],
                          Cell.Type_from ~ Cell.Type_to,
                          value.var = "scaled_count")
  
  connect_matrix_2d <- acast(total_all[total_all$type == "2D",],
                             Cell.Type_from ~ Cell.Type_to,
                             value.var = "scaled_count")
  
  h1 <- Heatmap(connect_matrix_3d,
                col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
                name = "Normalized\nEdge Count",
                column_title=paste(id, "3D"),
                cluster_rows = F, cluster_columns = F,
                row_order = cell_type_order,
                column_order = cell_type_order)
  
  h2 <- Heatmap(connect_matrix_2d,
                col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
                name = "Normalized\nEdge Count",
                column_title=paste(id, "2D"),
                cluster_rows = F, cluster_columns = F,
                row_order = cell_type_order,
                column_order = cell_type_order)
  
  pdf(paste0(out_dir, id, ".compared_cell_type_adj_matrix.pdf"), width=14, height=7)
  print(h1+h2)
  dev.off()
  
}

# make heatmap per condition

for (cond in conditions) {
  
  sample_3d <- connect_3d[connect_3d$condition == cond,]
  sample_2d <- connect_2d[connect_2d$condition == cond,]
  
  total_3d <- ddply(sample_3d,
                    .(Cell.Type_from,Cell.Type_to),
                    summarise,
                    cell_count=length(from)) 
  total_3d$type <- "3D"
  
  total_2d_slice <- ddply(sample_2d,
                          .(Cell.Type_from,Cell.Type_to, slice),
                          summarise,
                          slice_cell_count=length(from)) 
  
  total_2d <- ddply(total_2d_slice, 
                    .(Cell.Type_from,Cell.Type_to),
                    summarise,
                    cell_count=mean(slice_cell_count))
  total_2d$type <- "2D"
  
  total_all <- rbind(total_3d,
                     total_2d)
  
  total_all$scaled_count <- scale(log(total_all$cell_count))
  
  connect_matrix_3d <- acast(total_all[total_all$type == "3D",],
                             Cell.Type_from ~ Cell.Type_to,
                             value.var = "scaled_count")
  
  connect_matrix_2d <- acast(total_all[total_all$type == "2D",],
                             Cell.Type_from ~ Cell.Type_to,
                             value.var = "scaled_count")
  
  h1 <- Heatmap(connect_matrix_3d,
                col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
                name = "Normalized\nEdge Count",
                column_title=paste(cond, "3D"),
                cluster_rows = F, cluster_columns = F,
                row_order = cell_type_order,
                column_order = cell_type_order)
  
  h2 <- Heatmap(connect_matrix_2d,
                col=colorRamp2(c(-1.5,0,1.5), c("blue","white","red")),
                name = "Normalized\nEdge Count",
                column_title=paste(cond, "2D"),
                cluster_rows = F, cluster_columns = F,
                row_order = cell_type_order,
                column_order = cell_type_order)
  
  pdf(paste0(out_dir, cond, ".compared_cell_type_adj_matrix.pdf"), width=14, height=7)
  print(h1+h2)
  dev.off()
}



# do overall adjacency comparison

sample_order <- c("C164B WT","C158B APPPS19",
                  "C166A WT", "C165 APPPS19")

connect_3d_sums <- ddply(connect_3d,
                         .(sample_id, from),
                         summarise,
                         cell_count=length(unique(to)))
connect_3d_sums$type <- "3D"

connect_2d_sums <- ddply(connect_2d,
                         .(sample_id, from),
                         summarise,
                         cell_count=length(unique(to)))
connect_2d_sums$type <- "2D"

connect_sums <- rbind(connect_3d_sums,
                      connect_2d_sums)

connect_sums$sample_id <- factor(connect_sums$sample_id,
                                 levels=sample_order)


ggplot(connect_sums,
       aes(x=cell_count,
           fill=type)) +
  geom_density(adjust=3) +
  facet_wrap(~ sample_id, ncol=2) +
  theme_bw() +
  scale_fill_manual(values=c("grey","maroon")) +
  scale_x_continuous(breaks=seq(0,40,by=5)) +
  theme(legend.position = "bottom") +
  labs(x="Number of Connected Cells",
       y="Density of Cells",
       fill=NULL)
ggsave(paste0(out_dir, "mean_connected_cells.density_per_sample.png"),
       width=6, height=5)






