library(dyntoy)
library(magrittr)
library(tibble)
library(dplyr)

# Specify root folder for analysis
root_folder <- "scSTEM/analysis/"

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

# dataset: dynverse dataset object
# traj: dynverse trajectory object
# all_paths: paths from root to all leave nodes
# tmp_folder: temporary folder to save inputs for STEM.
get_stem_input<-function(norm_expressions, gene_ids, traj, all_paths, path_names, top_milestone){
  
  gr <- igraph::graph_from_edgelist(as.matrix(traj$milestone_network[,c("from","to")]))
  df <- list()
  
  ###############################################################################################
  # Assign cells to each path.
  # For paths having too few branching points, we will split the segments based on pseudotime.
  ###############################################################################################
  for(i in 1:length(all_paths)){
    
    path_nodes <- all_paths[[i]]
    path_cells <- list()
    df[[i]] <- list()
    deg <- igraph::degree(gr, v=path_nodes[2:length(path_nodes)],mode = "out")
    key_node_pos <- c(1,which(deg >= 2)+1,length(deg)+1) # index of root -> branching nodes ... ->leaf node
    
    # Assign cells to current path based on different scenario
    if(length(key_node_pos) <= 2){
      
      nodes <- path_nodes[key_node_pos[1]:key_node_pos[2]]
      cells <- top_milestone$cell_id[top_milestone$milestone_id %in% nodes]
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
        cells <- top_milestone$cell_id[top_milestone$milestone_id %in% nodes]
        pt <- sort(traj$pseudotime[cells])
        med_pos <- as.integer(median(1:length(pt)))
        path_cells[[length(path_cells)+1]] <- names(pt[1:med_pos])
        path_cells[[length(path_cells)+1]] <- names(pt[(med_pos+1):length(pt)])
      }
    }else{
      for(k in 1:(length(key_node_pos)-1)){
        nodes <- path_nodes[key_node_pos[k]:key_node_pos[k+1]]
        path_cells[[k]] <- top_milestone$cell_id[top_milestone$milestone_id %in% nodes]
      }
    }
    
    # Compute metrics using cells assigned to path
    for(k in 1:length(path_cells)){
      df[[i]][[k]] <- Matrix::rowMeans(Matrix::t(norm_expressions)[gene_ids, path_cells[[k]]])
    }
    names(df[[i]]) <- c(paste0("segment",1:length(df[[i]])))
  }
  
  names(df) <- path_names
  return(df)
}

stem_analysis<-function(
  input_name,
  tmp_folder,
  stem_path,
  setting_template_path,
  species
){
  
  # Read setting file template
  settings <- readLines(setting_template_path)
  
  # Specify species for GO annotations
  settings[3] <- paste0("Gene_Annotation_Source\t",species)
  
  # Run stem in command line mode
  # Decide what operating system current R is running on.
  osname = Sys.info()['sysname']
  
  # Specify input file path
  settings[2] <- paste0("Data_File\t",input_name)
  
  # Write current setting file
  writeLines(settings, file.path(tmp_folder, "setting"))
  # Enclose the file paths so that the spaces would not affect.
  if(osname == "Windows"){
    cmd = paste("java", "-mx1024M", "-jar",
                            paste0('"',stem_path,'"'),
                            "-b", paste0('"',paste0(tmp_folder, "/setting"),'"'),
                             '"/"')
  }else{
    
    # on linux/MAC we need to escape the spaces.
    stem_path = gsub(" ","\\\ ", stem_path, fixed = TRUE)
    cmd = paste("java", "-mx1024M", "-jar",
                stem_path,
                "-b", paste0(tmp_folder, "/setting"),
                '"/"')
  }
  
  # Run STEM clustering using commandline
  system(cmd)
  
  # remove current setting file
  unlink(paste0(tmp_folder, "setting"))
}

traj_types = c("linear","bifurcating","multifurcating")
de_rates = c(0.1,0.3,0.5)
for(traj_type in traj_types){
    for(de_rate in de_rates){
      # Run STEM
      cat("=========================================\n")
      cat(sprintf("Running simulations for %s %s\n",traj_type, de_rate))
      cat("=========================================\n")
      
      # Generate trajectory structure
      traj <- dyntoy::generate_trajectory(id = "aa",
                                          model = traj_type,
                                          num_cells = 10000
      )
      
      # Generate counts by zinb model
      counts <- dyntoy:::generate_counts( traj,
                                          num_features = 3000,
                                          dropout_rate = 0.9,
                                          dropout_probability_facto = 15000,
                                          differentially_expressed_rate = de_rate
      )

      # Get root milestone from the trajectory
      deg <- igraph::graph_from_edgelist(el = as.matrix(traj$milestone_network[,c("from","to")])) %>%
        igraph::degree( mode = "in")
      root_milestone <- names(deg[deg == 0])
      
      # Get closest milestone for each cell
      top_milestone <- traj$milestone_percentages %>%
        dplyr::group_by(cell_id) %>%
        dplyr::arrange(dplyr::desc(percentage)) %>%
        dplyr::slice_head(n = 1)
      
      # Get root cells
      root_cells <- top_milestone$cell_id[top_milestone$milestone_id == root_milestone]
      
      # Add root and pseudotime to the trajectory object
      traj <- dynwrap::add_root(traj,root_milestone = root_milestone, flip_edges = F) %>%
        dynwrap::add_pseudotime()
      
      # Get all paths from the trajectory tree
      all_paths <- get_all_paths(traj = traj, root_id = root_milestone)
      
      # Get cells mapped to each path
      all_path_cells <- list()
      for(path_node_ids in all_paths){
        all_path_cells[[length(all_path_cells)+1]] <- top_milestone$cell_id[top_milestone$milestone_id %in% path_node_ids]
      }
      names(all_path_cells) <- names(all_paths)
      
      # Get cell positisons in 2D space
      dimred <- dynwrap::calculate_trajectory_dimred(traj, adjust_weights = FALSE)
      pos <- as.matrix(dimred$cell_positions[,c(2,3)])
      rownames(pos) <- dimred$cell_positions$cell_id
      
      # Plot trajectory and cells using the 2D coordinates
      dynplot::plot_dimred(traj,
                           color_cells = "pseudotime",
                           dimred = pos,
                           size_milestones = 3,
                           size_cells = 4,
                           size_transitions = 1.5,
                           border_radius_percentage = 0
                          )
      # Save trajectory plot
      outfile <- file.path(root_folder,"results",paste0("trajPlot","_",traj_type,"_","DE",de_rate,".svg"))
      ggplot2::ggsave(outfile, width = 9, height = 11, dpi = 300)
      norm_expressions <- dynnormaliser::normalise_filter_counts(counts$counts)$expression
        
      # Get summarized expression data for STEM
      df <- get_stem_input(
        norm_expressions = norm_expressions,
        gene_ids = counts$tde_overall$feature_id,
        traj = traj,
        all_paths = all_paths,
        path_names = names(all_paths),
        top_milestone = top_milestone
      )
        
      # save summarized expressions as tsv files
      file_names <- file.path(root_folder, "data", "input", paste0(names(df),".tsv"))
      
      for(i in 1:length(df)){
        data <- do.call(cbind, df[[i]])
        data <- data.frame(gene_id = rownames(data), data)
        write.table(data,
                    file_names[i],
                    row.names = F,
                    col.names = T,
                    quote=F,
                    sep = "\t",
                    na = "")
      }

      outdir_name = file.path(root_folder, "results")
      cluster_res <- list()
      for(i in 1:length(file_names)){
        # Run STEM and get output table
        stem_analysis(
          input_name = file_names[i],
          tmp_folder = outdir_name,
          stem_path = file.path(root_folder,"code","stem.jar"),
          setting_template_path = file.path(root_folder,"code","stem_setting_template"),
          species = "Human (EBI)" # Ignore this for simulated data
        )
        
        # Read STEM clustering output file. If they do not exist, this means all genes were filtered out by STEM for current path
        prof_file <- file.path(outdir_name,"setting_profiletable.txt")
        gene_file <- file.path(outdir_name,"setting_genetable.txt")
        
        if(file.exists(prof_file) & file.exists(gene_file)){
          prof_tab <- read.csv(prof_file, sep = "\t") %>% as_tibble()
          gene_tab <- read.csv(gene_file, sep = "\t") %>% as_tibble()
          
          # Get genes assigned to significant clusters 
          tmp_tab <- left_join(gene_tab, prof_tab, by = c("Profile" = "Profile.ID"))
          tmp_tab <- tmp_tab[tmp_tab$Cluster...1.non.significant. != -1,]
          tmp_tab <- tmp_tab[,c("gene_id", "Cluster...1.non.significant.")]
          tmp_tab$path <- names(df)[i]
          colnames(tmp_tab) <- c("gene_id", "cluster", "path")
          
          cluster_res[[length(cluster_res) + 1]] <- tmp_tab
          
          # Remove output files.
          unlink(file.path(outdir_name,"setting_genetable.txt"))
          unlink(file.path(outdir_name,"setting_profiletable.txt"))
        }
        
        # Remove setting file and STEM input file.
        unlink(file.path(outdir_name,gsub(".tsv$","",file_names[i])))
        unlink(file.path(outdir_name, file_names[i]))
      }
      cluster_res <- do.call(rbind, cluster_res)
      
      # Save cluster result
      filename <- file.path(root_folder,"data","input","simulation",paste0("cluster_",traj_type,"_DE",de_rate,".csv"))
      write.csv(cluster_res,
                filename,
                quote = F,
                row.names = F)
      
      # Save trajectory object
      filename <- file.path(root_folder,"data","input","simulation",paste0("traj_",traj_type,"_DE",de_rate,".RDS"))
      saveRDS(traj, filename)
      
      # Save count matrix
      filename <- file.path(root_folder,"data","input","simulation",paste0("counts_",traj_type,"_DE",de_rate,"_.mtx"))
      Matrix::Matrix(counts$counts, sparse = T) %>%
      Matrix::writeMM(filename)
      
      # Save DEG info
      filename <- file.path(root_folder,"data","input","simulation",paste0("geneMeta_",traj_type,"_DE",de_rate,"_.csv"))
      gene_meta <- counts$tde_overall
      colnames(gene_meta) <- c("gene_id","DEG")
      write.csv(gene_meta,
                filename,
                quote = F,
                row.names = F)
      
      # Save cell info, including root cell info and pseudotime info.
      cell_meta <- data.frame(cell_id = rownames(counts$counts), time_point = 1)
      cell_meta$time_point[cell_meta$cell_id %in% root_cells] <- 0
      cell_meta$pseudotime <- traj$pseudotime[cell_meta$cell_id] # Add pseudotime to cell info table
      filename <- file.path(root_folder,"data","input","simulation",paste0("cellMeta_",traj_type,"_DE",de_rate,"_.csv"))
      write.csv(cell_meta,
                filename,
                quote = F,
                row.names = F)
      
      cat("Done\n")
    }
} 
