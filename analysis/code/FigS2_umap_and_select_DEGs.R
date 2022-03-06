library(monocle3)

# Spleen data downloaded from https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/
counts <- readRDS("../data/input/human_spleen/spleen_counts.RDS")
cell_meta <- readRDS("../data/input/human_spleen/df_cell.RDS")
cell_meta <- cell_meta[!is.na(cell_meta$Main_cluster_name),]
cell_meta <- cell_meta[cell_meta$Organ == "Spleen",]
rownames(cell_meta) <- cell_meta$sample
gene_meta <- readRDS("../data/input/human_spleen/df_gene.RDS")
gene_meta$gene_id <- gsub("[.].*$","",gene_meta$gene_id) # Strip off the '.' and the version number in gene id
rownames(gene_meta) <- gene_meta$gene_id
rownames(counts) <- gsub("[.].*$","",rownames(counts))

counts <- counts[,cell_meta$sample]

# Generate cds object for monocle3
cds <- new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)

# Perform PCA 100 dims followed by UMAP
cds <- monocle3::preprocess_cds(cds, num_dim = 100)
cds <- monocle3::reduce_dimension(cds, reduction_method = 'UMAP')

# Now select top 3000 DEGs
res <- graph_test(cds)
degs <- res[order(res$q_value),][1:3000,]$gene_id

# Subset cds object using the top 3000 DEGs and redo analysis
sub_cds <- cds[degs,]
sub_cds <- monocle3::preprocess_cds(sub_cds, num_dim = 100)
sub_cds <- monocle3::reduce_dimension(sub_cds, reduction_method = 'UMAP')
sub_cds <- monocle3::cluster_cells(sub_cds)

library(ggplot2)

# Visualize by UMAP
p <- monocle3::plot_cells(
  sub_cds,
  color_cells_by = "Main_cluster_name",
  show_trajectory_graph = F,
  label_cell_groups = T,
  group_label_size = 6,
  label_leaves=TRUE,
  cell_size = 0.8,
  label_branch_points=TRUE,
  graph_label_size=3,
  label_groups_by_cluster=FALSE
) + 
  theme(axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18))

ggsave("../results/figS2A_umap.svg",
        width = 7,
        height = 7)

# Generate inputs for scSTEM
colnames(cell_meta) <- gsub("sample", "cell_id", colnames(cell_meta))
colnames(cell_meta) <- gsub("Development_day", "time_point", colnames(cell_meta))
cell_meta <- cell_meta[order(cell_meta$time_point),]
cell_meta$time_point <- factor(cell_meta$time_point, labels = 1:length(unique(cell_meta$time_point)))
gene_meta <- gene_meta[degs,]
sub_cds <- sub_cds[gene_meta$gene_id, cell_meta$cell_id]

# Save input files for scSTEM. These can be processed by scSTEM GUI.
library(Matrix)
writeMM(exprs(sub_cds),
        "../data/input/human_spleen/counts_3kgenes_final.mtx")
write.csv(cell_meta,
          "../data/input/human_spleen/cell_meta_3kgenes_final.csv")
write.csv(gene_meta[degs,],
          "../data/input/human_spleen/gene_meta_3kgenes_final.csv")
