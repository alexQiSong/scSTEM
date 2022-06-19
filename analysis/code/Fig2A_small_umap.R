library(Matrix)
library(monocle3)

# Read data
counts <- readMM("../data/input/mouse_blood/counts_5k_final.mtx")
cell_meta <- read.csv("../data/input/mouse_blood/cell_meta_5kgenes_final.csv")
gene_meta <- read.csv("../data/input/mouse_blood/gene_meta_5kgenes_final.csv")

gene_names <- gene_meta$gene_id
rownames(counts) <- gene_meta$gene_id
colnames(counts) <- cell_meta$cell_id
rownames(gene_meta) <- gene_meta$gene_id
rownames(cell_meta) <- cell_meta$cell_id

# Preprocess and run UMAP
cds <- monocle3::new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)
cds <- monocle3::preprocess_cds(cds, num_dim = 100)
cds <- monocle3::reduce_dimension(cds, reduction_method = 'UMAP')
cds <- monocle3::cluster_cells(cds)

# Select the cells clusters presented in the paper for analysis. The results might look different on different OS.
p <- plot_cells(
  cds,
  color_cells_by = "cluster",
  show_trajectory_graph = F,
  label_cell_groups = T,
  group_label_size = 15,
  label_leaves=TRUE,
  cell_size = 0.8,
  label_branch_points=TRUE,
  graph_label_size=3,
  label_groups_by_cluster=FALSE
)
print(p)

# Get cell ids
selected_cells = names(clusters(cds))[clusters(cds) %in% c(7,8,14,15)]

# Get UMAP coordinates
dimred <- SingleCellExperiment::reducedDim(cds, "UMAP") %>%
  magrittr::set_colnames(c("comp_1","comp_2"))
dimred <- dimred[rownames(dimred) %in% selected_cells,]

# Get UMAP coordinates for cells
dimred <- SingleCellExperiment::reducedDim(cds, "UMAP") %>%
  magrittr::set_colnames(c("comp_1","comp_2"))
dimred <- dimred[rownames(dimred) %in% selected_cells,]

# Replot UMAP of this subset using customized color assignments
cell_plot_data <- as_tibble(dimred, rownames = NA) %>%
  dplyr::mutate(cell_type = cell_meta[rownames(dimred),]$Main_cell_type)

# Label cells
rows = !cell_plot_data$cell_type %in% c("Definitive erythroid lineage", "Primitive erythroid lineage","White blood cells")
cell_plot_data$cell_type[rows] = "Other"

# Plot UMAP results
plt <- ggplot2::ggplot() +
  ggplot2::geom_point( # Plot definitive erythroid lineage
    ggplot2::aes(
      x = comp_1,
      y = comp_2,
    ),
    data = cell_plot_data[cell_plot_data$cell_type == "Definitive erythroid lineage",],
    size = 3,
    colour = "#73a7fa"
  ) +
  ggplot2::geom_point( # Plot Primitive erythroid lineage
    ggplot2::aes(
      x = comp_1,
      y = comp_2,
    ),
    data = cell_plot_data[cell_plot_data$cell_type == "Primitive erythroid lineage",],
    size = 3,
    colour = "#f78c2f"
  ) +
  ggplot2::geom_point( # Plot White blood cells
    ggplot2::aes(
      x = comp_1,
      y = comp_2,
    ),
    data = cell_plot_data[cell_plot_data$cell_type == "White blood cells",],
    size = 3,
    colour = "#ed88fc"
  ) +
  ggplot2::geom_point( # Plot all other cells
    ggplot2::aes(
      x = comp_1,
      y = comp_2,
    ),
    data = cell_plot_data[cell_plot_data$cell_type == "Other",],
    size = 3,
    colour = "#FDDBD0"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.line = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank()
  )
print(plt)
