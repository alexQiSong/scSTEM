library("monocle3")
library("Matrix")
library("dplyr")
library("dynamicTreeCut")
library("ggplot2")

###############################################
# 1. read data
###############################################
counts <- readMM("../data/input/human_immune/counts_3kgenes_final.mtx")
cell_meta <- read.csv("../data/input/human_immune/cell_meta_3kgenes_final.csv")
gene_meta <- read.csv("../data/input/human_immune/gene_meta_3kgenes_final.csv")
gene_meta$gene_id <- gsub("[.].*$","",gene_meta$gene_id) # Strip off the '.' and the version number in gene id
gene_names <- gene_meta$gene_id
rownames(counts) <- gene_meta$gene_id
colnames(counts) <- cell_meta$cell_id
rownames(gene_meta) <- gene_meta$gene_id
rownames(cell_meta) <- cell_meta$cell_id

cds <- monocle3::new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)
cds <- monocle3::preprocess_cds(cds, num_dim = 100)
cds <- monocle3::reduce_dimension(cds, reduction_method = 'UMAP')
cds <- monocle3::cluster_cells(cds)

# View cells
p <- plot_cells(
  cds,
  color_cells_by = "partition",
  group_cells_by = "partition",
  show_trajectory_graph = F,
  label_cell_groups = T,
  group_label_size = 15,
  label_leaves=TRUE,
  cell_size = 0.8,
  label_branch_points=TRUE,
  graph_label_size=3,
  label_groups_by_cluster=FALSE
)
ggsave("../results/human_immune_umap.svg",width = 7, height = 7, dpi = 300)
###############################################
# 2. Build trajectory using monocle3
#    and get cells mapped to trajectory paths
###############################################
cds <- monocle3::learn_graph(cds,
                             use_partition = T,
                             close_loop = F)

# Use this cell cluster, the human immune cell cluster.
# The cluster number may be different on different OS.
pars <- monocle3::partitions(cds)
selected_cells <- names(pars[pars == 4])

# Cell meta file contains time point information
tp_table <- cell_meta[selected_cells,]

# Get principal graph
gr <- monocle3::principal_graph(cds)[["UMAP"]]

# Get closest cells for all milestone nodes in the principal graph
closest_nodes <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex[selected_cells,] %>%
  as.character() %>%
  paste0("Y_",.)

# Use the connected components that have nodes mapped to the selected cells
selected_comp <- igraph::components(gr)$membership %>%
  .[unique(closest_nodes)] %>%
  unique()

# Get milestone nodes in the selected connected component
selected_nodes <- igraph::components(gr)$membership %in% selected_comp %>%
  igraph::components(gr)$membership[.] %>%
  names()

# Get UMAP coordinates for cells and milestone nodes
# Use only the nodes and cells for current selected partitions/components
dimred <- SingleCellExperiment::reducedDim(cds, "UMAP") %>%
  magrittr::set_colnames(c("comp_1","comp_2"))
dimred <- dimred[rownames(dimred) %in% selected_cells,]

dimred_milestones <- t(cds@principal_graph_aux$UMAP$dp_mst) %>%
  magrittr::set_colnames(colnames(dimred))
dimred_milestones <- dimred_milestones[rownames(dimred_milestones) %in% selected_nodes,]

# Get milestone network
milestone_network <-
  igraph::induced_subgraph(gr, v = selected_nodes) %>%
  igraph::as_data_frame() %>%
  dplyr::transmute(
    from,
    to,
    length = sqrt(rowSums((dimred_milestones[from, ] - dimred_milestones[to, ])^2)),
    directed = FALSE
  )

# Milestone percentage as 1 for all cell-node pairs.
milestone_percentages <- data.frame(cell_id = selected_cells,
                                    milestone_id = closest_nodes,
                                    percentage = 1,
                                    stringsAsFactors = F) %>%
  as_tibble()

dimred_segment_progressions <-
  milestone_network %>%
  dplyr::select(from, to) %>%
  dplyr::mutate(percentage = purrr::map(seq_len(dplyr::n()), ~ c(0, 1))) %>%
  tidyr::unnest(percentage) %>%
  as_tibble(rownames = NA)

dsp_names <-
  dimred_segment_progressions %>%
  {ifelse(.$percentage == 0, .$from, .$to)}
dimred_segment_points <- dimred_milestones[dsp_names, , drop = FALSE]

# Wrap up trajectory as a dynverse trajectory object
traj <- dynwrap::wrap_data(cell_ids = selected_cells) %>%
  dynwrap::add_trajectory(
    milestone_ids = selected_nodes,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages
  ) %>%
  dynwrap::add_dimred(
    dimred = dimred,
    dimred_milestones = dimred_milestones,
    dimred_segment_progressions = dimred_segment_progressions,
    dimred_segment_points = dimred_segment_points,
    project_trajectory = F
  )

# Function to get root node from a trajectory object and a time point table (the earliest time point as root)
get_root_node <- function(traj, tp_table){
  
  # Assign cells to closest milestone node (highest percentage)
  mp <- traj$milestone_percentages
  mp <- mp[with(mp, order(cell_id,percentage,decreasing = T)),]
  closest_ms <- mp[match(unique(mp$cell_id), mp$cell_id),]
  
  # Iterate over time points from the earliest to the last
  for(tp in unique(sort(tp_table$time_point))){
    tp_cell_ids <- tp_table[tp_table$time_point == tp,]$cell_id
    if(any(tp_cell_ids %in% closest_ms$cell_id)){
      root_id <- closest_ms[closest_ms$cell_id %in% tp_cell_ids,]$milestone_id %>%
        table() %>%
        which.max() %>%
        names()
      
      # Return the node id of the root node
      return(root_id)
    }
  }
}

# add pseudotime to the model
root_id <- get_root_node(traj, tp_table)
cds <- monocle3::order_cells(cds, root_pr_nodes = root_id)
traj <- dynwrap::add_pseudotime(
  traj,
  monocle3::pseudotime(cds)[selected_cells]
) %>%
  dynwrap::add_root(
    root_milestone_id = root_id
  )

# Function for getting all paths from root node to all leave nodes
get_all_paths <- function(traj, root_id){
  
  root_id <- as.character(root_id)
  # if network is a directed network
  if(any(traj$milestone_network$directed)){
    
    # Get all leaf nodes and exclude the root node
    gr <- igraph::graph_from_edgelist(as.matrix(traj$milestone_network[,c("from","to")]))
    deg_out <- igraph::degree(gr, mode = "out")
    leaf_ids <- names(deg_out[deg_out == 0])
    leaf_ids <- leaf_ids[!(leaf_ids %in% root_id)]
    
    # Get all paths from root node to all leaf nodes
    paths <- igraph::shortest_paths(gr, root_id, leaf_ids, mode = 'out')$vpath
    non_empty_paths <- list()
    
    # If network is an undirected network
  }else{
    # Get all leaf nodes and exclude the root node
    gr <- igraph::graph_from_edgelist(as.matrix(traj$milestone_network[,c("from","to")]))
    deg<- igraph::degree(gr, mode = "all")
    leaf_ids <- names(deg[deg == 1])
    leaf_ids <- leaf_ids[!(leaf_ids %in% root_id)]
    
    # Get all paths from root node to all leaf nodes
    paths <- igraph::shortest_paths(gr, root_id, leaf_ids, mode = 'all')$vpath
    non_empty_paths <- list()
  }
  
  # Remove empty paths (paths going to unreachable leaf nodes)
  for(i in 1:length(paths)){
    if(length(paths[[i]]) > 0 & igraph::as_ids(paths[[i]])[1] == root_id){
      non_empty_paths[[leaf_ids[i]]] <- paths[[i]]$name
    }
  }
  names(non_empty_paths) <- paste0("path",1:length(non_empty_paths))
  return(non_empty_paths)
}

# Get node IDs for each path
all_paths <- get_all_paths(traj = traj, root_id = traj$root_milestone_id)

# Assign cells to closest milestone node (highest percentage)
mp <- traj$milestone_percentages
mp <- mp[with(mp, order(cell_id,percentage,decreasing = T)),]
closest_ms <- mp[match(unique(mp$cell_id), mp$cell_id),]

# Get cells mapped to each path
all_path_cells <- list()
for(path_node_ids in all_paths){
  all_path_cells[[length(all_path_cells)+1]] <- closest_ms$cell_id[closest_ms$milestone_id %in% path_node_ids]
}
names(all_path_cells) <- names(all_paths)

# Get a simplified trajectory tree
traj_simp <- dynwrap::simplify_trajectory(traj)

# Get coordinates of trajectory projections
traj_proj <- dynwrap::project_trajectory(
  trajectory = traj_simp,
  dimred = traj_simp$dimred
)

# Get coordinates of milestones
ms_plot_data <- as.data.frame(traj_proj$dimred_milestones)

# Remove cells not in the partition.
# The cluster number may be different on different OS.
pars <- monocle3::partitions(cds)
cell_coord <- SingleCellExperiment::reducedDim(cds,'UMAP')
cell_coord <- cell_coord[pars %in% 4,]

# Function for plotting the path and cells mapped to the path
plot_cell_and_path <- function(path){
  
  # Plot the paths in the main figure and cells mapped to these paths.
  path_cells <- all_path_cells[[path]]
  path_nodes <- all_paths[[path]]
  
  # Get groupings of the trajectory edges
  idx <- grep("BEGIN|END",rownames(traj_proj$dimred_segment_points))
  lens <- idx[seq(2,length(idx),2)] - idx[seq(1,length(idx),2)] + 1
  edge_group <- seq(1:length(lens))
  edge_group <- rep(edge_group, times = lens)
  
  # Get segments mapped to currently selected path
  idx_path <- lapply(
    1:length(idx),
    function(i){
      begin_end <- sub(
        "MILESTONE_BEGIN_W|MILESTONE_END_W",
        "",
        rownames(traj_proj$dimred_segment_points)[idx[i]]
      )
      begin_node <- substr(begin_end,1,(nchar(begin_end)-1)/2)
      end_node <- substr(begin_end,(nchar(begin_end)-1)/2+2,nchar(begin_end))
      if(begin_node %in% path_nodes &
         end_node %in% path_nodes){
        return(idx[i])
      }
    }
  ) %>% unlist()
  
  idx_path <- idx_path[!is.null(idx_path)]
  idx_path <- lapply(
    seq(1,length(idx_path),2),
    function(i){
      idx_path[i]:idx_path[i+1]
    }
  ) %>% unlist()
  
  # Grouping for trajectory path color
  path_group <- rep("other",length(edge_group))
  path_group[idx_path] <- "path"
  
  # Make data for plotting trajectory
  traj_plot_data <- data.frame(
    traj_proj$dimred_segment_points,
    edge_group = edge_group,
    path_group = path_group
  ) %>%
    as_tibble()
  
  # Get cell coordinates and cells mapped to currently selected path
  cell_plot_data <- as_tibble(traj_simp$dimred, rownames = NA)
  mask <- rownames(cell_plot_data) %in% path_cells
  cell_groups <- rep("other",length(mask))
  cell_groups[mask] <- "path"
  cell_plot_data <- dplyr::mutate(cell_plot_data, cell_groups = cell_groups)
  
  plt <- ggplot2::ggplot() +
    ggplot2::geom_point( # Plot other cells not mapped to the current path
      ggplot2::aes(
        x = comp_1,
        y = comp_2,
      ),
      data = cell_plot_data[cell_plot_data$cell_groups == "other",],
      size = 1.5,
      colour = "#9C9C9C"
    ) +
    ggplot2::geom_point( # Plot cells mapped to the current path
      ggplot2::aes(
        x = comp_1,
        y = comp_2,
      ),
      data = cell_plot_data[cell_plot_data$cell_groups == "path",],
      size = 1.5,
      colour = "#FFE120"
    ) +
    ggplot2::geom_path( # Plot other edges not mapped to the current path
      ggplot2::aes(
        x = comp_1,
        y = comp_2,
        group = edge_group
      ),
      data = traj_plot_data[traj_plot_data$path_group == "other",],
      size = 3,
      color = "#333333"
    ) +
    ggplot2::geom_path( # Plot edges mapped to the current path
      ggplot2::aes(
        x = comp_1,
        y = comp_2,
        group = edge_group
      ),
      data = traj_plot_data[traj_plot_data$path_group == "path",],
      size = 3,
      color = "#FF0000"
    ) +
    ggplot2::geom_point( # Plot milestones
      ggplot2::aes(
        x = comp_1,
        y = comp_2,
      ),
      data = ms_plot_data,
      size = 3,
      colour = "#9C9C9C"
    ) +
    ggplot2::geom_point( # Plot milestone borders
      ggplot2::aes(
        x = comp_1,
        y = comp_2,
      ),
      data = ms_plot_data,
      size = 3,
      shape = 1,
      stroke = 2,
      colour = "#333333"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "black"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )
  return(plt)
}

# Plot paths
for(i in 1:9){
  path_to_plot <- paste0("path",i)
  p <- plot_cell_and_path(path_to_plot)
  out_file <- paste0("../results/","human_immune_path",i,".svg")
  ggsave(out_file, width = 9, height = 11, dpi = 300)
}

# Use cells in this path (NK cell path)
# This may be slightly different for different OS
path_cells <- all_path_cells[["path1"]]

gene_names = gene_meta$gene_id
###############################################
# 3. Perform clustering analysis
###############################################
###############################################
# 3.1 scHOT clustering analysis
library("scHOT")
# Use this to keep clustering results, scSTEM results were obtained from
# running scSTEM GUI.
cluster_res <- data.frame(method = c("scHot","tradeSeq","scSTEM"),
                          time = c(0,0,60),
                          clusters = c(0,0,2),
                          ratio = c(0,0,7/117)
)

# Select the top 300 genes and compute weighted spearman correlation
start_time <- proc.time()
cor = matrix(0, ncol = 300, nrow = 300)
for(i in 1:300){
  for(j in i:300){
    cor[i,j] = weightedSpearman(counts[i,],counts[j,])
  }
}
cor[lower.tri(cor)] = t(cor)[lower.tri(t(cor))]

# Hierarchical clustering based on correlation matrix
dendro <- hclust(dist(cor, method = "maximum"))
cluster_labels <- cutreeDynamic(dendro)
end_time <- proc.time()

# Record the running time for scHOT
cluster_res[cluster_res$method == "scHot","time"] <- (end_time - start_time)["elapsed"]

# Perform GO analysis for each cluster using clusterprofiler
# and get top 10 go terms.
library(clusterProfiler)
library(org.Hs.eg.db)
ratios <- c()
for(cluster in unique(cluster_labels)){
  cluster_genes <- gene_meta$gene_id[1:300][cluster_labels %in% cluster]
  ego <- enrichGO(gene          = cluster_genes,
                  universe      = gene_meta$gene_id[1:300],
                  keyType       = "ENSEMBL",
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  res <- ego@result
  gene_counts <- strsplit(res$GeneRatio,"/") %>%
              unlist(.) %>%
              .[c(T,F)] %>%
              as.numeric(.)
  cluster_size <- strsplit(res$GeneRatio,"/") %>%
            unlist(.) %>%
            .[c(F,T)] %>%
            as.numeric(.)
  background_counts <- strsplit(res$BgRatio,"/") %>%
            unlist(.) %>%
            .[c(T,F)] %>%
            as.numeric(.)
  background_size <- strsplit(res$BgRatio,"/") %>%
            unlist(.) %>%
            .[c(F,T)] %>%
            as.numeric(.)
  enrich_fold <- gene_counts/(cluster_size/background_size * background_counts)
  res$enrich_fold <- enrich_fold
  res$gene_counts <- gene_counts
  
  if(sum(res[res$p.adjust < 0.05,]$Description %>% grepl(".*natural killer.*",.)) > 0){
    schot_nk_genes <- res[res$p.adjust < 0.05,]$Description %>%
      grepl(".*natural killer.*",.) %>%
      res[.,"core_enrichment"] %>%
      strsplit("/") %>%
      unlist() %>%
      unique()
  }else{
    schot_nk_genes = c()
  }
  
  # Compute NK gene ratio
  ratios<-c(ratios,length(tradeseq_nk_genes)/length(cluster_genes))
  
  # View Top 10
  res <- res[order(res$enrich_fold,decreasing = T),][1:10,]
  res$Description <- stringr::str_wrap(res$Description, width = 35) # Format GO names.
  res$Description <- factor(res$Description, levels = res$Description) # Format GO names.
  
  print(res)
}

# Check NK gene ratios. None of the clusters were enriched for NK genes
print(ratios)

# Record NK gene ratios for scHOT
cluster_res[cluster_res$method == "scHot","ratio"] <- 0

# Record number of clusters for scHOT
cluster_res[cluster_res$method == "scHot","clusters"] <- length(unique(cluster_labels))

###############################################
# 3.2 tradeSeq clustering analysis
library(tradeSeq)
tmp_counts <- counts[,path_cells]
tmp_pseudotime <- monocle3::pseudotime(cds)
tmp_pseudotime <- tmp_pseudotime[path_cells]
tmp_pseudotime <- matrix(tmp_pseudotime, ncol = 1, nrow = length(tmp_pseudotime))
cellWeights <- matrix(1, ncol = 1, nrow = length(tmp_pseudotime))

start_time <- proc.time()
sce <- fitGAM(counts = as.matrix(tmp_counts),
              pseudotime = tmp_pseudotime,
              cellWeights = cellWeights)

# Run clustering
library(clusterExperiment)
library(ggplot2)
nPointsClus <- 20
clusPat <- clusterExpressionPatterns(sce, nPoints = 20,
                                     genes = rownames(counts))
end_time <- proc.time()

# Record the running time for tradeSeq
cluster_res[cluster_res$method == "tradeSeq","time"] <- (end_time - start_time)["elapsed"]

clusterLabels <- primaryCluster(clusPat$rsec)
cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq == -1] # remove unclustered genes

# Record number of clusters for scHOT
cluster_res[cluster_res$method == "tradeSeq","clusters"] <- length(unique(cUniq))

# View cluster sizes for the top five clusters
table(clusterLabels[clusterLabels %in% cUniq[1:5]])

# Perform GO analysis and plot the GO enrichment for the top five clusters.
library(clusterProfiler)
library(org.Hs.eg.db)
for(cluster in cUniq[1:5]){
  cluster_genes <- gene_names[clusterLabels == cluster]
  ego <- enrichGO(gene          = cluster_genes,
                  universe      = gene_names,
                  keyType       = "ENSEMBL",
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  res <- ego@result
  gene_counts <- strsplit(res$GeneRatio,"/") %>%
    unlist(.) %>%
    .[c(T,F)] %>%
    as.numeric(.)
  cluster_size <- strsplit(res$GeneRatio,"/") %>%
    unlist(.) %>%
    .[c(F,T)] %>%
    as.numeric(.)
  background_counts <- strsplit(res$BgRatio,"/") %>%
    unlist(.) %>%
    .[c(T,F)] %>%
    as.numeric(.)
  background_size <- strsplit(res$BgRatio,"/") %>%
    unlist(.) %>%
    .[c(F,T)] %>%
    as.numeric(.)
  enrich_fold <- gene_counts/(cluster_size/background_size * background_counts)
  res$enrich_fold <- enrich_fold
  res$gene_counts <- gene_counts
  
  if(sum(res[res$p.adjust < 0.05,]$Description %>% grepl(".*natural killer.*",.)) > 0){
    tradeseq_nk_genes <- res[res$p.adjust < 0.05,]$Description %>%
      grepl(".*natural killer.*",.) %>%
      res[.,"core_enrichment"] %>%
      strsplit("/") %>%
      unlist() %>%
      unique()
  }else{
    tradeseq_nk_genes = c()
  }
  
  # View top 10
  res <- res[order(res$enrich_fold,decreasing = T),][1:10,]
  res$Description <- stringr::str_wrap(res$Description, width = 35) # Format GO names.
  res$Description <- factor(res$Description, levels = res$Description) # Format GO names.
  print(res)
  
  # Compute NK gene ratio
  ratios<- length(tradeseq_nk_genes)/length(cluster_genes)
}

# Check NK gene ratios. None of the clusters were enriched for NK genes
print(ratios)

# Record NK gene ratios for tradeSeq
cluster_res[cluster_res$method == "tradeSeq","ratio"] <- 0

###################################################
# 4. Differential Expression (DE) Analsysis
###################################################
###############################################
# 4.1 tradeSeq DE analysis

# Use this to keep DE results, scSTEM results were obtained from
# running scSTEM GUI.
de_res <- data.frame(method = c("Monocle3","tradeSeq","scSTEM"),
                          time = c(0,0,60),
                          clusters = c(0,0,2),
                          ratio = c(0,0,(7/117)*100)
)

library(tradeSeq)

start_time <- proc.time()
tradeseq_de <- startVsEndTest(sce, l2fc = log2(2))
tradeseq_de$pvalue <- p.adjust(tradeseq_de$pvalue, method = "fdr") # Correct pvalues
end_time <- proc.time()

# Running time includes GAM model fit time in the previous steps
de_res[de_res$method == "tradeSeq","time"] <- (end_time - start_time)["elapsed"] + cluster_res[cluster_res$method == "tradeSeq","time"]

# Perform GO analysis using pvalues and significant genes
library(clusterProfiler)
library(org.Hs.eg.db)
geneList <- tradeseq_de$pvalue
names(geneList) <- gene_names
geneList <- 1-geneList[order(geneList)]
geneList <- rank(geneList) # Replace pvalues with ranking
tradeseq_ego <- gseGO(geneList = geneList,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENSEMBL",
                      ont = "BP",
                      readable = T
                      )

# Rank by enrichment scpre
res <- tradeseq_ego@result

# Get NK related genes
tradeseq_nk_genes <- res[res$p.adjust < 0.05,]$Description %>%
  grepl(".*natural killer.*",.) %>%
  res[.,"core_enrichment"] %>%
  strsplit("/") %>%
  unlist() %>%
  unique()

# Get tradeseq DEGs
tradeseq_de_genes <- rownames(tradeseq_de[tradeseq_de$pvalue < 0.05,])

# Compute NK gene ratio among all DE genes
de_res[de_res$method == "tradeSeq","ratio"] <- length(intersect(tradeseq_nk_genes, tradeseq_de_genes))/length(tradeseq_de_genes)

# Top 10 for plotting
res <- res[order(res$enrichmentScore,decreasing = T),][1:10,]
res$Description <- stringr::str_wrap(res$Description, width = 35) # Format GO names.
res$Description <- factor(res$Description, levels = res$Description) # Format GO names.

library(ggplot2)
# This is figS1B bottom left
out_file <- "../result/figS1B_tradeseq_path1_GO.svg"
p <- ggplot(data = res, aes(x = Description, y = setSize, fill = p.adjust)) +
  geom_bar(stat = "identity") + 
  scale_fill_continuous(name = "P-value", limits = c(0,0.05), breaks = c(0, 0.05)) +
  theme_minimal(base_size = 22) + 
  theme(legend.position="bottom") +
  coord_flip() +
  xlab("GO category") +
  ylab("Counts")
ggsave(out_file, width = 9, height = 11, dpi = 300)

###############################################
# 4.2 Monocle3 DE analysis

# Subset the cds object. Get the subset for path one cells
path1_cds <- cds[,colData(cds)$cell_id %in% path_cells]

start_time <- proc.time()
monocle3_de <- graph_test(path1_cds, neighbor_graph="knn")
end_time <- proc.time()

# Record running time
de_res[de_res$method == "Monocle3","time"] <- (end_time - start_time)["elapsed"]

geneList <- monocle3_de$p_value
geneList[is.na(geneList)] <- 1
names(geneList) <- gene_names
geneList <- 1-geneList
geneList <- rank(geneList) # Replace pvalues with ranking
geneList <- geneList[order(geneList,decreasing = T)]
monocle3_ego <- gseGO(geneList = geneList,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENSEMBL",
                      ont = "BP"
                      )

# Rank by enrichment scores
res <- monocle3_ego@result

# Get NK related genes
monocle3_nk_genes <- res[res$p.adjust < 0.05,]$Description %>%
  grepl(".*natural killer.*",.) %>%
  res[.,"core_enrichment"] %>%
  strsplit("/") %>%
  unlist() %>%
  unique()

# Get monocle3 DEGs
monocle3_de_genes <- monocle3_de[monocle3_de$q_value < 0.05,]$gene_id

# Compute NK gene ratio
de_res[de_res$method == "Monocle3","ratio"] <- length(intersect(monocle3_nk_genes, monocle3_de_genes))/length(monocle3_de_genes)

# Top 10 for plotting
res <- res[order(res$enrichmentScore,decreasing = T),][1:10,]
res$Description <- stringr::str_wrap(res$Description, width = 35) # Format GO names.
res$Description <- factor(res$Description, levels = res$Description) # Format GO names.

# This is figS1B bottom right
out_file <- "../result/figS1B_monocle3_path1_GO.svg"
p <- ggplot(data = res, aes(x = Description, y = setSize, fill = p.adjust)) +
  geom_bar(stat = "identity") + 
  scale_fill_continuous(name = "P-value", limits = c(0,0.05), breaks = c(0, 0.05)) +
  theme_minimal(base_size = 22) + 
  theme(legend.position="bottom") +
  coord_flip() +
  xlab("GO category") +
  ylab("Counts")
ggsave(out_file, width = 9, height = 11, dpi = 300)

###################################################
# 5. Plot benchmark results
###################################################
###################################################
# 5.1 Clustering results

p1 <- ggplot(data = cluster_res[,c("method","time")],
             aes(x = method,
                 y = time,
                 fill = method)
) +
  scale_fill_manual(values = c(scHot = "bisque1",
                               tradeSeq = "coral",
                               scSTEM = "skyblue")) +
  geom_bar(stat = "identity") + 
  theme_minimal(base_size = 25) +
  theme(legend.position="top") +
  labs(fill = NULL) +
  xlab(NULL) +
  ylab("Running time (seconds)") +
  ylim(c(0,50000)) + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

p2 <- ggplot(data = cluster_res[,c("method","clusters")],
             aes(x = method, y = clusters, fill = method)) +
  scale_fill_manual(values = c(scHot = "bisque1",
                               tradeSeq = "coral",
                               scSTEM = "skyblue")) +
  geom_bar(stat = "identity") + 
  theme_minimal(base_size = 25) +
  theme(legend.position="top") +
  labs(fill = NULL) +
  xlab(NULL) +
  ylab("Number of filtered clusters") + 
  ylim(c(0,70))

p3 <- ggplot(data = cluster_res[,c("method","ratio")],
             aes(x = method, y = ratio, fill = method)) +
  scale_fill_manual(values = c(scHot = "bisque1",
                               tradeSeq = "coral",
                               scSTEM = "skyblue")) +
  geom_bar(stat = "identity") + 
  theme_minimal(base_size = 25) +
  theme(legend.position="top") +
  labs(fill = NULL) +
  xlab(NULL) +
  ylab("NK gene ratio") +
  scale_y_continuous(labels = scales::percent)

p <- gridExtra::grid.arrange(p1, p2, p3,
                             nrow = 1,
                             top = grid::textGrob("Clustering Analysis",gp=grid::gpar(fontsize=27)))
ggsave("../result/figS1A_clusteringComparison.svg",
       p,
       width = 16,
       height = 5,
       dpi = 300)

###################################################
# 5.2 DE results
p1 <- ggplot(data = de_res[,c("method","time")],
             aes(x = method,
                 y = time,
                 fill = method)
) +
  scale_fill_manual(values = c(Monocle3 = "darkgoldenrod1",
                               tradeSeq = "coral",
                               scSTEM = "skyblue")) +
  geom_bar(stat = "identity") + 
  theme_minimal(base_size = 25) +
  theme(legend.position="top") +
  labs(fill = NULL) +
  xlab(NULL) +
  ylab("Running time (seconds)") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

p2 <- ggplot(data = de_res[,c("method","ratio")],
             aes(x = method, y = ratio, fill = method)) +
  scale_fill_manual(values = c(Monocle3 = "darkgoldenrod1",
                               tradeSeq = "coral",
                               scSTEM = "skyblue")) +
  geom_bar(stat = "identity") + 
  theme_minimal(base_size = 25) +
  theme(legend.position="top") +
  labs(fill = NULL) +
  xlab(NULL) +
  ylab("NK gene ratio") +
  scale_y_continuous(labels = scales::percent, limits = c(0,0.1) )

p <- gridExtra::grid.arrange(p1, p2,
                             nrow = 1,
                             top = grid::textGrob("Differential Expression Analysis",gp=grid::gpar(fontsize=27)))
ggsave("../result/figS1B_DEComparison.svg",
       p,
       width = 12,
       height = 5,
       dpi = 300)
