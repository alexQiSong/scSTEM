library("monocle3")
library("Matrix")
library("dplyr")
library("dynamicTreeCut")
library("ggplot2")

# Read data from human immune cell data set.
counts <- readMM("../data/input/human_immune/counts_3kgenes_final.mtx")
cell_meta <- read.csv("../data/input/human_immune/cell_meta_3kgenes_final.csv")
gene_meta <- read.csv("../data/input/human_immune/gene_meta_3kgenes_final.csv")

gene_meta$gene_id <- gsub("[.].*$","",gene_meta$gene_id) # Strip off the '.' and the version number in gene id
gene_names <- gene_meta$gene_id
rownames(counts) <- gene_meta$gene_id
colnames(counts) <- cell_meta$cell_id
rownames(gene_meta) <- gene_meta$gene_id
rownames(cell_meta) <- cell_meta$cell_id

#############################################
# Normalization, PCA and UMAP
#############################################
cds <- new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)

# Perform PCA followed by UMAP
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, reduction_method = 'UMAP')

# Cluster cells
cds <- cluster_cells(cds)

#############################################
# Get immune cell cluster and its trajectory
#############################################
# Infer trajectories using monocle3
cds <- monocle3::learn_graph(cds,
                             use_partition = T,
                             close_loop = F)

# Use cells in cluster 1, the immune cell cluster.
pars <- monocle3::partitions(cds)
selected_cells <- names(pars[pars == 1])

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

#############################################
# Format trajectory as dynverse trajectory model
#############################################
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

#############################################
# Map cells to trajectory paths
#############################################

# Get node IDs for each path
all_paths <- get_all_paths(traj = traj, root_id = traj$root_milestone_id)

# Assign cells to closest milestone node (highest percentage)
mp <- traj$milestone_percentages
mp <- mp[with(mp, order(cell_id,percentage,decreasing = T)),]
closest_ms <- mp[match(unique(mp$cell_id), mp$cell_id),]

# Get cells mapped to each segment.
all_path_cells <- list()
for(i in 1:length(all_paths)){
  
  path_nodes <- all_paths[[i]]
  path_cells <- list()
  all_path_cells[[i]] <- list()
  
  # Find branching nodes by degrees
  deg <- igraph::degree(gr, v=path_nodes[2:length(path_nodes)],mode = "all")
  key_node_pos <- c(1,which(deg > 2)+1,length(deg)+1) # index of root -> branching nodes ... ->leaf node
  
  # Assign cells to current path based on different scenario
  if(length(key_node_pos) <= 1 ){
    print(paste0("path",i," is too short for analysis, skipped..."))
    next
  }else if(length(key_node_pos) == 2){
    
    nodes <- path_nodes[key_node_pos[1]:key_node_pos[2]]
    cells <- closest_ms$cell_id[closest_ms$milestone_id %in% nodes]
    pt <- sort(traj$pseudotime[cells])
    
    key_pos <- quantile(
      1:length(pt),
      probs = c(0,1/3,2/3,1)
    ) %>%
      as.integer()
    
    for(k in 1:(length(key_pos)-1)){
      path_cells[[k]] <-  names(pt[key_pos[k]:key_pos[k+1]])
    }
    
  }else if(length(key_node_pos) == 3){
    for(k in 1:(length(key_node_pos)-1)){
      nodes <- path_nodes[key_node_pos[k]:key_node_pos[k+1]]
      cells <- closest_ms$cell_id[closest_ms$milestone_id %in% nodes]
      pt <- sort(traj$pseudotime[cells])
      med_pos <- as.integer(median(1:length(pt)))
      path_cells[[length(path_cells)+1]] <- names(pt[1:med_pos])
      path_cells[[length(path_cells)+1]] <- names(pt[(med_pos+1):length(pt)])
    }
  }else{
    for(k in 1:(length(key_node_pos)-1)){
      nodes <- path_nodes[key_node_pos[k]:key_node_pos[k+1]]
      path_cells[[k]] <- closest_ms$cell_id[closest_ms$milestone_id %in% nodes]
    }
  }
  all_path_cells[[i]] <- path_cells
  names(all_path_cells[[i]]) <- c(paste0("segment",1:length(all_path_cells[[i]])))
}

# Get uniq path cells
is_terminal <- c() # Use this to indicate whether the given segment is at the terminal
all_path_cells_uniq <- list()
for(i in 1:length(all_path_cells)){
  for(j in 1:length(all_path_cells[[i]])){
    if(length(all_path_cells_uniq) == 0){
      all_path_cells_uniq[[length(all_path_cells_uniq)+1]] <- all_path_cells[[i]][[j]]
      is_terminal <- c(is_terminal,F)
    }else{
      found = F
      for(k in 1:length(all_path_cells_uniq)){
        n_cells <- intersect(all_path_cells[[i]][[j]], all_path_cells_uniq[[k]]) %>%
          length()
        if(n_cells == length(all_path_cells_uniq[[k]])){
          found = T
          break
        }
      }
      if(!found){
        all_path_cells_uniq[[length(all_path_cells_uniq)+1]] <- all_path_cells[[i]][[j]]
        if(j == length(all_path_cells[[i]])){
          is_terminal <- c(is_terminal,T)
        }else{
          is_terminal <- c(is_terminal,F)
        }
      }
    }
  }
}

#############################################
# Fit linear model to path cells
#############################################
all_pvals <- list()
for(k in 1:length(all_path_cells_uniq)){
  
  # Perform filtering for x and Y
  x <- traj$pseudotime[all_path_cells_uniq[[k]]]
  x <- x[!is.infinite(x)]
  Y <- normalized_counts(cds)[gene_names, names(x)]
  
  # Use genes expressed in >5% cells
  filtered_Y <- Y[rowSums(Y != 0)/ncol(Y) > 0.05,]
  
  # Perform linear regressions
  res<-lapply(1:nrow(filtered_Y), function(j){
    lmModel <- lm(filtered_Y[j,] ~ x)
    f <- summary(lmModel)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    r2 <- summary(lmModel)$r.squared
    return(c(p,r2))
  }) %>%
    unlist() %>%
    matrix(ncol = 2, byrow = T)
  pvals <- p.adjust(res[,1],method = "fdr")
  r2s <- res[,2]
  all_pvals[[k]] <- pvals
}

# Check
rates <- c()
num_genes <- c()
for(k in 1:length(all_pvals)){
  rates <- c(rates,sum(all_pvals[[k]] < 0.05)/length(all_pvals[[k]]))
  num_genes <- c(num_genes, length(all_pvals[[k]]))
}

ord <- order(-rates)
res <- tibble(rates = rates[ord],
              is_terminal = is_terminal[ord],
              num_genes = num_genes[ord],
              index = 1:length(rates)
)

#############################################
# Generate final plot
#############################################
library(ggplot2)
p<-ggplot(data=res, aes(x=index, y=rates, fill = is_terminal)) +
  geom_bar(stat="identity") + 
  geom_text(aes(label=num_genes), vjust=-0.3, size=7) +
  theme_minimal(base_size = 25) +
  theme(axis.text.x = element_blank()) +
  xlab(NULL) +
  ylab("Percentage of significant genes") + 
  scale_y_continuous(labels = scales::percent) +
  scale_fill_discrete(name = "Terminal segment")

ggsave("../results/figS11_linearModel.svg",
       p,
       width = 14,
       height = 7,
       dpi = 300)
