# library(monocle3)
# library(magrittr)
# library(tibble)
# library(SingleCellExperiment)
# library(magrittr)
# library(hdf5r), this package is needed for infer_trajectory()
# package 'hexbin' is needed
# Other dependencies: ROGUE; shinyFiles; shinyWidgets, R.utils
# dynplot package has some bugs, to fix this, do: devtools::install_github("dynverse/dynplot@devel", force = TRUE)

# traj: a trajectory object.
# dataset: a dataset object.
# tp_table: time point information table
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

# Rewire the milestone network by specified root milestone
rewire_by_root <- function(traj, root_id){
  return(traj %>% dynwrap::add_root(root_milestone_id = root_id))
}

# Extract all paths from root node to all leave nodes
# traj: dynverse trajectory object.
# root_id: root node id.
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


# Compute the entropy scores
get_ds <- function(exp_mat, gene_ids, cell_ids){
  ds <- rep(NA,length(gene_ids))
  names(ds) <- gene_ids
  #dss <- tryCatch({
  #  ROGUE::SE_fun(exp_mat[gene_ids, cell_ids])$ds
  #},error = function(x){
  #  ds
  #})
  dss <- ROGUE::SE_fun(exp_mat[gene_ids, cell_ids])$ds
  ds[names(dss)] <- dss
  return(abs(ds))
}

# Compute change rate
get_cr <- function(ptime, norm_exp_mat, gene_ids, cell_ids){

  x <- ptime[cell_ids]
  Y <- norm_exp_mat[gene_ids, cell_ids]

  # Compute z scores for expressions > 0.
  # and filter expressions using z scores
  mask <- Y>0
  means <- Matrix::rowSums(Y)/Matrix::rowSums(mask)
  sd <- sqrt(Matrix::rowSums(((Y-means)*mask)^2)/Matrix::rowSums(mask))
  z <- (Y-means)/sd
  z_filter <- z < sd * 1.96 & z > sd * -1.96 & mask

  # Use entries of which  sd * -1.96 < z-score < sd * 1.96 and expression > 0 as y
  # for linear model and pseudotime as x for linear model.
  rates <- unlist(lapply(1:nrow(Y), function(i, z_filter, x){
    z_select = !is.na(z_filter[i,]) & (z_filter[i,] == T)
    if(any(z_select)){
      return(lm(Y[i, z_select] ~ order(x[z_select]))$coefficients[2])
    }else{
      return(NA)
    }
  },z_filter, x))
  names(rates) <- rownames(Y)
  return(rates)
}

# Get cells closest to a path
# traj: dynverse trajectory object
# path_node_ids: node ids for a path
get_path_cells <- function(traj, closest_ms, path_node_ids){

  # Assign cells to closest milestone node (highest percentage)
  mp <- traj$milestone_percentages
  mp <- mp[with(mp, order(cell_id,percentage,decreasing = T)),]
  closest_ms <- mp[match(unique(mp$cell_id), mp$cell_id),]

  # Assign cells to path nodes
  return(closest_ms$cell_id[closest_ms$milestone_id %in% path_node_ids])
}

# dataset: dynverse dataset object
# traj: dynverse trajectory object
# all_paths: paths from root to all leave nodes
# tmp_folder: temporary folder to save inputs for STEM.
get_stem_input<-function(dataset, traj, all_paths, path_names, closest_ms, metric){

  # Add pseudotime if metric = 'cr2'
  if(metric == 'cr2'){
    traj <- traj %>% dynwrap::add_pseudotime(pseudotime = dynwrap::calculate_pseudotime(traj))
  }


  gr <- igraph::graph_from_edgelist(as.matrix(traj$milestone_network[,c("from","to")]))
  df <- list()
  gene_ids <- dataset$feature_ids

  shiny::incProgress(0, detail = sprintf("%g%% done",0))

  ###############################################################################################
  # Assign cells to each path.
  # For paths having two few branching points, we will split the segments based on pseudotime.
  ###############################################################################################
  for(i in 1:length(all_paths)){

    path_nodes <- all_paths[[i]]
    path_cells <- list()
    df[[i]] <- list()
    if(traj$directed){
      deg <- igraph::degree(gr, v=path_nodes[2:length(path_nodes)],mode = "out")
      key_node_pos <- c(1,which(deg >= 2)+1,length(deg)+1) # index of root -> branching nodes ... ->leaf node
    }else{
      deg <- igraph::degree(gr, v=path_nodes[2:length(path_nodes)],mode = "all")
      key_node_pos <- c(1,which(deg > 2)+1,length(deg)+1) # index of root -> branching nodes ... ->leaf node
    }

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

    # Compute metrics using cells assigned to path
    for(k in 1:length(path_cells)){
      if(metric == "entropy_reduction"){
        df[[i]][[k]] <- get_ds(exp_mat = Matrix::t(dataset$counts),
                               gene_ids = gene_ids,
                               cell_ids = path_cells[[k]])
      }else if(metric == "mean"){
        df[[i]][[k]] <- Matrix::rowMeans(Matrix::t(dataset$expression)[gene_ids, path_cells[[k]]])
      }else if(metric == "median"){
        df[[i]][[k]] <- apply(Matrix::t(dataset$expression)[gene_ids, path_cells[[k]]], 1, median)
      }else if(metric == "change_rate"){
        df[[i]][[k]] <- get_cr(ptime = traj$pseudotime,
                               norm_exp_mat = Matrix::t(dataset$expression),
                               gene_ids = gene_ids,
                               cell_ids = path_cells[[k]])
      }
    }
    names(df[[i]]) <- c(paste0("segment",1:length(df[[i]])))

    # Keep track of progress
    shiny::incProgress(1/length(all_paths), detail = sprintf("%g%% done",round(100*i/length(all_paths),2)))
  }

  #names(cells) <- paste0("path",1:length(cells))
  names(df) <- path_names
  return(df)
}

#' scSTEM GUI function
#'
#' This function runs scSTEM GUI
#'
#' To use scSTEM in R, Run `run_scstem_GUI()` in R to call a shiny app based GUI. Please follow the steps described in user manual to perform scSTEM analysis.
#' @return NULL
#' @export
run_scstem_GUI <- function(){

  ui <- shiny::fluidPage(
    shiny::fluidRow(
      # The left side panel
      shiny::column(6,

                    ######## Step 1. Load expression matrix and cluster cells
                    shiny::wellPanel(
                      shiny::fluidRow(
                        shiny::column(12,
                                      align = 'center',
                                      shiny::h4("Step 1: Load input files", style="text-align:left;font-weight:bold")
                                      ),
                        shiny::column(12,
                                      align = "left",

                                      # Select count matrix file
                                      shiny::fixedRow(
                                        shiny::column(3,
                                                      align = "left",
                                                      shinyFiles::shinyFilesButton(id = "exp_file",
                                                                                   style = "margin-top: 10px;",
                                                                                   label = "Count matrix",
                                                                                   title = "Please select a *.mtx file for count matrix",
                                                                                   multiple = F,
                                                                                   icon = shiny::icon("search",lib = "glyphicon"))),
                                        shiny::column(6,
                                                      align = "left",
                                                      shiny::h5(shiny::textOutput(outputId = "exp_file_status"), style = "text-align: left; padding: 20px, 0"))
                                                    ),

                                      # Select cell metadata file
                                      shiny::fixedRow(
                                        shiny::column(3,
                                                      align = "left",
                                                      shinyFiles::shinyFilesButton(id = "cell_file",
                                                                                   style = "margin-top: 10px;",
                                                                                   label = "Cell meta",
                                                                                   title = "Please select a *.csv file for cell meta data",
                                                                                   multiple = F,
                                                                                   icon = shiny::icon("search",lib = "glyphicon"))),
                                        shiny::column(6,
                                                      align = "left",
                                                      shiny::h5(shiny::textOutput(outputId = "cell_file_status"), style = "text-align: left; padding: 20px, 0"))
                                      ),

                                      # Select gene metadata file
                                      shiny::fixedRow(
                                        shiny::column(3,
                                                      align = "left",
                                                      shinyFiles::shinyFilesButton(id = "gene_file",
                                                                                   style = "margin-top: 10px;",
                                                                                   label = "Gene meta",
                                                                                   title = "Please select a *.csv file for gene meta data",
                                                                                   multiple = F,
                                                                                   icon = shiny::icon("search",lib = "glyphicon"))),
                                        shiny::column(6,
                                                      align = "left",
                                                      shiny::h5(shiny::textOutput(outputId = "gene_file_status"), style = "text-align: left; padding: 20px, 0"))
                                      ),

                                      # Load all input files here
                                      shiny::fixedRow(
                                        style = "margin-top: 10px;",
                                        shiny::column(3,
                                                      align = 'left',
                                                      shiny::selectInput(inputId = "id_type",
                                                                         choices = c("gene symbol","ensembl ID"),
                                                                         label = "Gene ID type"),
                                                                         icon = shiny::icon("folder-open",lib = "glyphicon")),
                                        shiny::column(3,
                                                      align = 'left',
                                                      shiny::selectInput(inputId = "species",
                                                                         choices = c("homo sapiens","mus musculus"),
                                                                         label = "Species"),
                                                      icon = shiny::icon("folder-open",lib = "glyphicon")),

                                        shiny::column(3,
                                                      align = "left",
                                                      shiny::actionButton(inputId = "load",
                                                                          style = "margin-top: 24px;",
                                                                          label = "Load files",
                                                                          icon = shiny::icon("folder-open",lib = "glyphicon"))
                                                      ),
                                        shiny::column(3,
                                                      align = "left",
                                                      shiny::actionButton(inputId = "load_sample",
                                                                          style = "margin-top: 24px;",
                                                                          label = "Load sample data",
                                                                          icon = shiny::icon("folder-open",lib = "glyphicon"))
                                        )
                                      ),

                                      shiny::fixedRow(
                                        shiny::column(12,
                                                      align = "left",
                                                      shiny::h5(shiny::textOutput(outputId = "loaded_status"), style = "text-align: left; padding: 20px, 0"))
                                      )
                                      ),

                        )
                        ),

                    ####### Step 2. Visualization and cell clustering
                    shiny::wellPanel(
                      shiny::fluidRow(
                        shiny::column(12,
                                      align = 'center',
                                      shiny::h4("Step 2: Visualization and cell clustering (Optional)", style="text-align:left;font-weight:bold")
                                      ),
                        shiny::column(3,
                                      align = "left",
                                      shiny::actionButton(inputId = "umap",
                                                          label = "Run UMAP",
                                                          )),
                        shiny::column(3,
                                      align = "left",
                                      shiny::actionButton(inputId = "cluster",
                                                          label = "Run Clustering",
                                                         )),
                        shiny::column(3,
                                      align = "left",
                                      shiny::actionButton(inputId = "vis_partition",
                                                          label = "Visualize Results",icon = shiny::icon("sunglasses",lib = "glyphicon")))
                        )
                    ),

                    ######## Step 3. Infer trajectory for the given partition using user selected method
                    shiny::wellPanel(
                      shiny::fluidRow(
                        shiny::column(12,
                                      align = 'center',
                                      shiny::h4("Step 3: Infer trajecotries", style="text-align:left;font-weight:bold")
                                      ),
                        shiny::column(3,
                                      align = "left",
                                      shinyWidgets::pickerInput(inputId = "partition_select",
                                                                label = "Partition",
                                                                choices = c("all"),
                                                                multiple = T)
                                      ),
                        shiny::column(3,
                                      align = "left",
                                      shiny::selectInput(
                                        inputId = "method",
                                        label = "Method",
                                        choices = c("slingshot",
                                                    "paga_tree",
                                                    "paga",
                                                    "slice",
                                                    "pcreode",
                                                    "celltree_maptpx",
                                                    "scuba",
                                                    "celltree_vem",
                                                    "sincell",
                                                    "raceid_stemid",
                                                    "elpigraph",
                                                    "celltrails",
                                                    "urd",
                                                    "cellrouter",
                                                    "celltree_gibbs",
                                                    "slicer",
                                                    "calista",
                                                    "monocle_ddrtree",
                                                    "monocle3"
                                                    )
                                      )
                        ),
                        shiny::column(3,
                                      align = "left",
                                      shiny::selectInput(
                                        inputId = 'use_partition',
                                        label = "Use partitions?",
                                        choices = c("Yes","No")
                                      )
                        ),
                        shiny::column(3,
                                      align = "left",
                                      style = "margin-top: 25px;",
                                      shiny::actionButton(inputId = "infer",
                                                          label = "Infer trajectory",
                                                          icon = shiny::icon("folder-open",lib = "glyphicon"))
                        ),
                      )
                    ),

                    ######## Step 4. Visualize paths
                    shiny::wellPanel(
                      shiny::fluidRow(

                          shiny::column(12,
                                        align = 'center',
                                        shiny::h4("Step 4: Visualize paths (optional)", style="text-align:left;font-weight:bold")
                                        ),

                          shiny::column(3,
                                        align = 'left',
                                        shiny::selectInput(
                                          inputId = "vis_path_select",
                                          label = "Select path",
                                          choices = c()
                                        )),
                          shiny::column(3,
                                        align = 'left',
                                        style = "margin-top: 25px;",
                                        shiny::actionButton(inputId = "vis_path_umap",
                                                            label = "View by UMAP",
                                                            ))

                        )
                      ),
                    ######## Step 5. Run STEM
                    shiny::wellPanel(
                      shiny::fluidRow(
                        shiny::column(12,
                                      align = 'center',
                                      shiny::h4("Step 5: Run STEM analysis", style = 'text-align:left;font-weight:bold')
                                      ),

                        shiny::column(3,
                                      align = "left",
                                      shinyWidgets::pickerInput(inputId = "run_path_select",
                                                         label = "Select path(s)",
                                                         choices = c(),
                                                         multiple = T
                                      )
                        ),

                        shiny::column(3,
                                      align = 'left',
                                      shiny::selectInput(
                                        inputId = "metric",
                                        label = "Metric",
                                        choices = c("mean","entropy_reduction","change_rate","meidan")
                                        )
                                      )

                      ),
                      shiny::fluidRow(
                        shiny::column(3,
                                      align = "left",
                                      style = "margin-top: 25px;",
                                      shinyFiles::shinyDirButton(
                                                                 id = "tmp_folder",
                                                                 label = "Output folder",
                                                                 title = "Select a folder for saving temporary files",
                                                                 multiple = F,
                                                                 icon = shiny::icon("search",lib = "glyphicon")
                                                                )
                        ),
                        shiny::column(6,
                                      align = "left",
                                      style = "margin-top: 25px;",
                                      shiny::textOutput(outputId = "outdir_name")
                                      )
                      ),
                      shiny::fluidRow(
                        shiny::column(3,
                                      align = "left",
                                      style = "margin-top: 25px",
                                      shiny::actionButton(
                                        inputId = "run_stem",
                                        label = "Run STEM",
                                        icon = shiny::icon("play",lib = "glyphicon")
                                      )
                        )
                      )
                      ),

                      ######## Step 6. Run cluster comparison (optional)
                      shiny::wellPanel(
                        shiny::fluidRow(
                          shiny::column(12,
                                        align = 'center',
                                        shiny::h4("Step 6: Run comparison of clusters (optional)", style = 'text-align:left;font-weight:bold')
                          ),

                          shiny::column(3,
                                        align = "left",
                                        shiny::selectInput(
                                          inputId = "compare1_name",
                                          label = "Path name 1",
                                          choices = c()
                                        )
                          ),

                          shiny::column(3,
                                        align = "left",
                                        shiny::selectInput(
                                          inputId = "compare2_name",
                                          label = "Path name 2",
                                          choices = c()
                                        )
                          ),
                        ),
                        shiny::fluidRow(
                          shiny::column(3,
                                        align = "left",
                                        style = "margin-top: 25px",
                                        shiny::actionButton(
                                          inputId = "run_comparison",
                                          label = "Run comparison",
                                          icon = shiny::icon("play",lib = "glyphicon")
                                        )
                                      )
                        )
                      )
                    ),
      # The right plot area
      shiny::column(6,
                    shiny::fluidRow(shiny::plotOutput("partition_plot")),
                    shiny::fluidRow(shiny::plotOutput("path_plot"))
      )
    )
  )

  server <- function(input, output){

    rv <- shiny::reactiveValues(counts = NULL,
                         cds = NULL,
                         dataset = NULL,
                         traj = NULL,
                         traj_simp = NULL,
                         traj_proj = NULL,
                         ms_plot_data = NULL,
                         species = NULL,
                         all_paths = NULL,
                         all_path_cells = NULL,
                         closest_ms = NULL,
                         metric = NULL,
                         cell_meta = NULL,
                         gene_meta = NULL,
                         outdir_name = NULL,
                         stem_folder = system.file('STEM',package = 'scSTEM'),
                         current_dir = "",
                         current_root = "")

    output$exp_file_status <- shiny::renderPrint({
      cat(sprintf("No file selected"))
    })
    output$cell_file_status <- shiny::renderPrint({
      cat(sprintf("No file selected"))
    })
    output$gene_file_status <- shiny::renderPrint({
      cat(sprintf("No file selected"))
    })
    output$outdir_name <- shiny::renderPrint({
      cat(sprintf("No output folder selected"))
    })

    ############################################################
    # Step 1. Load all input files
    ############################################################
    # Here user specifies expression file name to be loaded
    shiny::observeEvent(input$exp_file, {
      in_file <- shinyFiles::parseFilePaths(roots = shinyFiles::getVolumes(), input$exp_file)
      if(rv$current_dir == ""){
        shinyFiles::shinyFileChoose(input = input, id = "exp_file", roots = shinyFiles::getVolumes())
      }else{
        shinyFiles::shinyFileChoose(input = input, id = "exp_file", roots = shinyFiles::getVolumes(), defaultPath = rv$current_dir)
      }

      if(length(in_file$name)>0){

        # rv$current_dir stores the relative path relative to the current root.
        rv$current_dir <- gsub(paste0(in_file$name,"$"),"",in_file$datapath) %>%
                            R.utils::getRelativePath(shinyFiles::getVolumes()()[input$exp_file$root])
        output$exp_file_status <- shiny::renderPrint({
          cat(sprintf("Selected: %s",in_file$name))
        })
      }
    })

    # Here user specifies cell meta file name to be loaded
    shiny::observeEvent(input$cell_file, {
      in_file <- shinyFiles::parseFilePaths(roots = shinyFiles::getVolumes(), input$cell_file)
      if(rv$current_dir == ""){
        shinyFiles::shinyFileChoose(input = input, id = "cell_file", roots = shinyFiles::getVolumes())
      }else{
        shinyFiles::shinyFileChoose(input = input, id = "cell_file", roots = shinyFiles::getVolumes(), defaultPath = rv$current_dir)
      }

      if(length(in_file$name)>0){

        # rv$current_dir stores the relative path relative to the current root.
        rv$current_dir <- gsub(paste0(in_file$name,"$"),"",in_file$datapath) %>%
          R.utils::getRelativePath(shinyFiles::getVolumes()()[input$cell_file$root])
        output$cell_file_status <- shiny::renderPrint({
          cat(sprintf("Selected: %s",in_file$name))
        })
      }
    })

    # Here user specifies gene meta file name to be loaded
    shiny::observeEvent(input$gene_file, {
      in_file <- shinyFiles::parseFilePaths(roots = shinyFiles::getVolumes(), input$gene_file)
      if(rv$current_dir == ""){
        shinyFiles::shinyFileChoose(input = input, id = "gene_file", roots = shinyFiles::getVolumes())
      }else{
        shinyFiles::shinyFileChoose(input = input, id = "gene_file", roots = shinyFiles::getVolumes(), defaultPath = rv$current_dir)
      }

      if(length(in_file$name)>0){

        # rv$current_dir stores the relative path relative to the current root.
        rv$current_dir <- gsub(paste0(in_file$name,"$"),"",in_file$datapath) %>%
          R.utils::getRelativePath(shinyFiles::getVolumes()()[input$gene_file$root])
        output$gene_file_status <- shiny::renderPrint({
          cat(sprintf("Selected: %s",in_file$name))
        })
      }
    })

    shiny::observeEvent(input$load, {

      # Get path and file names
      exp_fname <- shinyFiles::parseFilePaths(roots = shinyFiles::getVolumes(), input$exp_file)
      cell_fname <- shinyFiles::parseFilePaths(roots = shinyFiles::getVolumes(), input$cell_file)
      gene_fname <- shinyFiles::parseFilePaths(roots = shinyFiles::getVolumes(), input$gene_file)

      # Check if files are already specified
      if(length(exp_fname$datapath) == 0 | length(cell_fname$datapath) == 0 | length(gene_fname$datapath) == 0){
        shiny::showModal(shiny::modalDialog(title = "All input files should be specified before loading input files.",
                                            footer = modalButton("OK"),
                                            easyClose = F))
      }else{

        output$loaded_status <- shiny::renderPrint({
          cat(sprintf("reading input files ..."))
        })

        # Tell the users files are being loaded.
        shiny::showModal(shiny::modalDialog(title = "Loading input files",
                                            footer = NULL,
                                            easyClose = F))

        # Read input files
        rv$counts <- Matrix::readMM(as.character(exp_fname$datapath))
        rv$cell_meta <- read.csv(as.character(cell_fname$datapath))
        rv$gene_meta <- read.csv(as.character(gene_fname$datapath))

        # Check if meta data meets the requirements
        if(!("cell_id" %in% colnames(rv$cell_meta)) | !("time_point" %in% colnames(rv$cell_meta))){
          shiny::removeModal()
          shiny::showModal(shiny::modalDialog(title = "Error: Cell meta data table should contain at least two columns named 'cell_id' and 'time_point'",
                                              footer = modalButton("OK"),
                                              easyClose = F))
          output$loaded_status <- renderPrint({cat(sprintf("Failed to load input files"))})
        }else if(!("gene_id" %in% colnames(rv$gene_meta))){
          shiny::removeModal()
          shiny::showModal(shiny::modalDialog(title = "Error: Gene meta data table should contain at least one column named 'gene_id'",
                                              footer = modalButton("OK"),
                                              easyClose = F))
          output$loaded_status <- renderPrint({cat(sprintf("Failed to load input files"))})
        }else if(length(rv$cell_meta$cell_id) != ncol(rv$counts)){
          shiny::removeModal()
          shiny::showModal(shiny::modalDialog(title = "Error: Number of rows in cell meta data should be equal to number of columns in count matrix",
                                              footer = modalButton("OK"),
                                              easyClose = F))
          output$loaded_status <- renderPrint({cat(sprintf("Failed to load input files"))})
        }else if(length(rv$gene_meta$gene_id) != nrow(rv$counts)){
          shiny::removeModal()
          shiny::showModal(shiny::modalDialog(title = "Error: Number of rows in gene meta data should be equal to number of rows in count matrix",
                                              footer = modalButton("OK"),
                                              easyClose = F))
          output$loaded_status <- renderPrint({cat(sprintf("Failed to load input files"))})
        }else if(!is.numeric(rv$cell_meta$time_point)){
          shiny::removeModal()
          shiny::showModal(shiny::modalDialog(title = "Error: column 'time_point' in cell meta data table should only contain numeric data",
                                              footer = modalButton("OK"),
                                              easyClose = F))
          output$loaded_status <- renderPrint({cat(sprintf("Failed to load input files"))})
        }else{

          # Files are loaded. Remove file loaidng pop-up window.
          shiny::removeModal()

          # Tell the users when ID conversion is going on.
          shiny::showModal(shiny::modalDialog(title = "Converting gene IDs...",
                                              footer = NULL,
                                              easyClose = F))

          # Convert gene id to gene names (for STEM to run GO analysis)
          map_info <- list(
            "homo sapiens" = list(dataset="hsapiens_gene_ensembl",
                                  stem_species="Human (EBI)",
                                  attr="hgnc_symbol"),
            "mus musculus" = list(dataset="mmusculus_gene_ensembl",
                                  stem_species="Mouse (EBI)",
                                  attr="mgi_symbol")
          )

          if(input$id_type == "ensembl ID"){
            rv$gene_meta$gene_id <- gsub("[.].*$","",rv$gene_meta$gene_id)
            mart <- biomaRt::useDataset(map_info[[input$species]]$dataset, biomaRt::useMart("ensembl"))
            converted <- biomaRt::getBM(filters= "ensembl_gene_id",
                                        attributes= c("ensembl_gene_id",map_info[[input$species]]$attr),
                                        values = rv$gene_meta$gene_id, # The version number along with '.' will be stripped off
                                        mart= mart)

            # Remove genes with duplicated gene symbols
            converted <- converted[!(duplicated(converted[,2]) | duplicated(converted[,2], fromLast = T)),]

            res <- dplyr::left_join(rv$gene_meta,
                                    converted,
                                    by = c("gene_id" = "ensembl_gene_id"))
            replace_ensembl <- !(is.na(res[,2]) | res[,2] == "")
            res$gene_id[replace_ensembl] <- res[replace_ensembl,2]
            rv$gene_meta$gene_id <- res$gene_id
          }
          rownames(rv$cell_meta) <- rv$cell_meta$cell_id
          rownames(rv$gene_meta) <- rv$gene_meta$gene_id
          rownames(rv$counts) <- rv$gene_meta$gene_id
          colnames(rv$counts) <- rv$cell_meta$cell_id

          # Gene ID conversion done. Remove pop-up window.
          shiny::removeModal()

          # Generate cds object using input files
          rv$cds <- monocle3::new_cell_data_set(
              expression_data = rv$counts,
              cell_metadata = rv$cell_meta,
              gene_metadata = rv$gene_meta
            )

          # Tell STEM what species is used for GO annotations.
          rv$species <- map_info[[input$species]]$stem_species

          output$loaded_status <- renderPrint({
            cat(sprintf("All files successfully loaded"))
          })
        }
      }
    })

    # If user chooses to load sample data set.
    shiny::observeEvent(input$load_sample, {

      # Get sample data folder
      sample_folder <- system.file('sample_data',package = 'scSTEM')

      output$loaded_status <- shiny::renderPrint({
        cat(sprintf("reading input files ..."))
      })

      # Tell the users files are being loaded.
      shiny::showModal(shiny::modalDialog(title = "Loading input files",
                                          footer = NULL,
                                          easyClose = F))

      # Read input files
      rv$counts <- Matrix::readMM(file.path(sample_folder,"counts_3kgenes_final.mtx"))
      rv$cell_meta <- read.csv(file.path(sample_folder,"cell_meta_3kgenes_final.csv"))
      rv$gene_meta <- read.csv(file.path(sample_folder,"gene_meta_3kgenes_final.csv"))

      # Files are loaded. Remove file loaidng pop-up window.
      shiny::removeModal()

      # Tell the users when ID conversion is going on.
      shiny::showModal(shiny::modalDialog(title = "Converting gene IDs...",
                                          footer = NULL,
                                          easyClose = F))

      # Convert gene id to gene names (for STEM to run GO analysis)
      map_info <- list(
        "homo sapiens" = list(dataset="hsapiens_gene_ensembl",
                              stem_species="Human (EBI)",
                              attr="hgnc_symbol"),
        "mus musculus" = list(dataset="mmusculus_gene_ensembl",
                              stem_species="Mouse (EBI)",
                              attr="mgi_symbol")
      )

      if(input$id_type == "ensembl ID"){
        rv$gene_meta$gene_id <- gsub("[.].*$","",rv$gene_meta$gene_id)
        mart <- biomaRt::useDataset(map_info[[input$species]]$dataset, biomaRt::useMart("ensembl"))
        converted <- biomaRt::getBM(filters= "ensembl_gene_id",
                                    attributes= c("ensembl_gene_id",map_info[[input$species]]$attr),
                                    values = rv$gene_meta$gene_id, # The version number along with '.' will be stripped off
                                    mart= mart)

        # Remove genes with duplicated gene symbols
        converted <- converted[!(duplicated(converted[,2]) | duplicated(converted[,2], fromLast = T)),]

        res <- dplyr::left_join(rv$gene_meta,
                                converted,
                                by = c("gene_id" = "ensembl_gene_id"))
        replace_ensembl <- !(is.na(res[,2]) | res[,2] == "")
        res$gene_id[replace_ensembl] <- res[replace_ensembl,2]
        rv$gene_meta$gene_id <- res$gene_id
      }
      rownames(rv$cell_meta) <- rv$cell_meta$cell_id
      rownames(rv$gene_meta) <- rv$gene_meta$gene_id
      rownames(rv$counts) <- rv$gene_meta$gene_id
      colnames(rv$counts) <- rv$cell_meta$cell_id

      # Gene ID conversion done. Remove pop-up window.
      shiny::removeModal()

      # Generate cds object using input files
      rv$cds <- monocle3::new_cell_data_set(
        expression_data = rv$counts,
        cell_metadata = rv$cell_meta,
        gene_metadata = rv$gene_meta
      )

      # Tell STEM what species is used for GO annotations.
      rv$species <- map_info[[input$species]]$stem_species

      output$loaded_status <- renderPrint({
        cat(sprintf("All files successfully loaded"))
      })

      # Update the file selection status
      output$exp_file_status <- shiny::renderPrint({
        cat(sprintf("Selected: %s","counts_3kgenes_final.mtx"))
      })

      output$cell_file_status <- shiny::renderPrint({
        cat(sprintf("Selected: %s","cell_meta_3kgenes_final.csv"))
      })

      output$gene_file_status <- shiny::renderPrint({
        cat(sprintf("Selected: %s","counts_3kgenes_final.mtx"))
      })

    })

    ############################################################
    # Step 2. Visualize partitions
    ############################################################
    shiny::observeEvent(input$umap, {

      shiny::showModal(shiny::modalDialog(title = "Running UMAP, please wait...",
                                          footer = NULL,
                                          easyClose = F))
      rv$cds <- monocle3::preprocess_cds(rv$cds, num_dim = 100)
      rv$cds <- monocle3::reduce_dimension(rv$cds, reduction_method = 'UMAP')
      shiny::removeModal()
      shiny::showModal(shiny::modalDialog(title = "UMAP is done.",
                                          footer = modalButton("OK"),
                                          easyClose = F))
    })

    shiny::observeEvent(input$cluster, {
      shiny::showModal(shiny::modalDialog(title = "Clustering cells, please wait...",
                                          footer = NULL,
                                          easyClose = F))
      rv$cds <- monocle3::cluster_cells(rv$cds)
      shiny::removeModal()
      shiny::showModal(shiny::modalDialog(title = "Clustering is done.",
                                          footer = modalButton("OK"),
                                          easyClose = F))

      # Update partition selector when clustering is done.
      shinyWidgets::updatePickerInput(session = shiny::getDefaultReactiveDomain(),
                                      inputId = "partition_select",
                                      choices = c(unique(monocle3::partitions(rv$cds)),"all"))
    })

    shiny::observeEvent(input$vis_partition, {

      # Check if UMAP has been done
        if(is.null(rv$cds)){
          shiny::showModal(shiny::modalDialog(title = "To visualize cells, please first run UMAP and clustering",
                                              footer = modalButton("OK"),
                                              easyClose = F))
        }else{
          shiny::showModal(shiny::modalDialog(title = "Generating plot...",
                                              footer = NULL,
                                              easyClose = F))
          shiny::removeModal()
          output$partition_plot <- shiny::renderPlot({
            monocle3::plot_cells(
                                   rv$cds,
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

          })
        }

    })

    ############################################################
    # Step 3. Infer trajectories
    ############################################################
    shiny::observeEvent(input$infer, {
        pars <- monocle3::partitions(rv$cds)

        # If users did not run UMAP and proceed directly to run inference.
        if("all" %in% input$partition_select){
          selected_cells <- names(pars)
        }else{
          selected_cells <- names(pars[pars %in% input$partition_select])
        }

        tp_table <- rv$cell_meta[selected_cells,]
        rv$dataset <- dynwrap::wrap_expression(
          counts = Matrix::t(rv$counts[, selected_cells]),
          expression = Matrix::t(monocle3::normalized_counts(rv$cds)[, selected_cells])
        )
        start_cells <- tp_table[
                                  tp_table$time_point == min(tp_table$time_point),
                               ]$cell_id
        rv$dataset <- rv$dataset %>%
                      dynwrap::add_prior_information(start_id = start_cells)
        shiny::showModal(shiny::modalDialog(title = "Inferring trajectories, please wait...",
                                            footer = NULL,
                                            easyClose = F))
        if(input$method == "monocle3"){
          rv$cds <- monocle3::learn_graph(rv$cds,
                                          use_partition = ifelse(input$use_partition == "Yes",T,F),
                                          close_loop = F)
          gr <- monocle3::principal_graph(rv$cds)[["UMAP"]]

          closest_nodes <- rv$cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex[selected_cells,] %>%
                           as.character() %>%
                           paste0("Y_",.)

          # Use the connected components that have nodes mapped to the selected cells
          selected_comp <- igraph::components(gr)$membership %>%
          .[unique(closest_nodes)] %>%
          unique()

          selected_nodes <- igraph::components(gr)$membership %in% selected_comp %>%
                            igraph::components(gr)$membership[.] %>%
                            names()

          # Get UMAP coordinates for cells and milestone nodes
          # Use only the nodes and cells for current selected partitions
          dimred <- SingleCellExperiment::reducedDim(rv$cds, "UMAP") %>%
            magrittr::set_colnames(c("comp_1","comp_2"))
          dimred <- dimred[rownames(dimred) %in% selected_cells,]

          dimred_milestones <- t(rv$cds@principal_graph_aux$UMAP$dp_mst) %>%
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
            #tibble::as_tibble()

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

          # Wrap up trajectory as dynverse trajectory object
          rv$traj <- dynwrap::wrap_data(cell_ids = selected_cells) %>%
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

          # add pseudotime to the model
          root_id <- get_root_node(rv$traj, tp_table)
          rv$cds <- monocle3::order_cells(rv$cds,
                                          root_pr_nodes = root_id)
          rv$traj <- dynwrap::add_pseudotime(
            rv$traj,
            monocle3::pseudotime(rv$cds)[selected_cells]
            ) %>%
            dynwrap::add_root(
              root_milestone_id = root_id
            )
        }else{
          ######################################################################
          # saveRDS(rv,"C:/Users/sqsq3/Documents/stem_analysis/tmp/mouse_data_neural_crest2/slingshot/rv.rds")
          # rr <- readRDS("C:/Users/sqsq3/Documents/stem_analysis/tmp/human_data_immune/slingshot/rv.rds")
          # rv$all_paths <- rr$all_paths
          # rv$dataset <- rr$dataset
          # rv$cds <- rr$cds
          # rv$traj <- rr$traj
          # rv$all_path_cells <- rr$all_path_cells
          # rv$dir_name <- '/home/alex/tmp/'
          # rv$counts <- rr$counts
          # rv$species <- rr$species
          # rv$closest_ms <- rr$closest_ms
          # mp <- rv$traj$milestone_percentages
          # mp <- mp[with(mp, order(cell_id,percentage,decreasing = T)),]
          # rv$closest_ms <- mp[match(unique(mp$cell_id), mp$cell_id),]
          # shiny::updateSelectInput(inputId = 'partition_select',
          #                         choices = unique(monocle3::partitions(rv$cds)))
          #####################################################################
          #saveRDS(rv,"C:/Users/sqsq3/Documents/stem_analysis/tmp/mouse_data_neural_crest2/slingshot/rv.rds")

          # Use dimred results from UMAP step for visualization later
          dimred <- SingleCellExperiment::reducedDim(rv$cds, "UMAP") %>%
            magrittr::set_colnames(c("comp_1","comp_2"))
          dimred <- dimred[rownames(dimred) %in% selected_cells,]

          # Infer trajectory and add calculate pseudotime
          rv$traj <- dynwrap::infer_trajectory(
            rv$dataset,
            input$method,
            give_priors = c("start_id"),
            verbose = T
            ) %>%
          dynwrap::add_root(
            root_milestone_id = get_root_node(., tp_table)
            ) %>%
          dynwrap::add_pseudotime(
            pseudotime = dynwrap::calculate_pseudotime(.)
            ) %>%
          dynwrap::add_dimred(
            dimred = dimred
          )
        }

        # Use this root node to rewire the trajectory network
        if(rv$traj$directed){
          rv$traj <- rewire_by_root(rv$traj, rv$traj$root_milestone_id)
        }

        # Get paths
        rv$all_paths <- get_all_paths(traj = rv$traj, root_id = rv$traj$root_milestone_id)

        # Assign cells to closest milestone node (highest percentage)
        mp <- rv$traj$milestone_percentages
        mp <- mp[with(mp, order(cell_id,percentage,decreasing = T)),]
        rv$closest_ms <- mp[match(unique(mp$cell_id), mp$cell_id),]

        # Get cells mapped to each path
        rv$all_path_cells <- list()
        for(path_node_ids in rv$all_paths){
           rv$all_path_cells[[length(rv$all_path_cells)+1]] <- rv$closest_ms$cell_id[rv$closest_ms$milestone_id %in% path_node_ids]
        }
        names(rv$all_path_cells) <- names(rv$all_paths)

        # Update path info when inference is done
        shiny::updateSelectInput(inputId = 'vis_path_select', choices = names(rv$all_paths))

        shinyWidgets::updatePickerInput(session = shiny::getDefaultReactiveDomain(),
                                        inputId = 'run_path_select',
                                        choices = names(rv$all_paths))

        shiny::updateSelectInput(inputId = 'compare1_name', choices = names(rv$all_paths))
        shiny::updateSelectInput(inputId = 'compare2_name', choices = names(rv$all_paths))

        # Get a simplified trajectory tree
        rv$traj_simp <- dynwrap::simplify_trajectory(rv$traj)
        rv$traj_simp$directed <- rv$traj$directed

        # Get coordinates of trajectory projections
        rv$traj_proj <- dynwrap::project_trajectory(
          trajectory = rv$traj_simp,
          dimred = rv$traj_simp$dimred
        )

        # Get coordinates of milestones
        rv$ms_plot_data <- as.data.frame(rv$traj_proj$dimred_milestones)

        shiny::removeModal()
        shiny::showModal(shiny::modalDialog(title = "Inference is done.",
                                            footer = modalButton("OK"),
                                            easyClose = F))
    })

    ################################################################
    # Step 4. Visualize the selected path and the selected partition
    ################################################################
    shiny::observeEvent(input$vis_path_umap, {

      # Check whether trajectory has been generated ?
      if(is.null(rv$traj) | is.null(rv$traj_proj) | is.null(rv$traj_simp)){
        shiny::showModal(shiny::modalDialog(title = "Trajectory has not been constructed. Please infer trajectory first.",
                                            footer = modalButton("OK"),
                                            easyClose = F))
      }else{

      # Tell the user that plot is being generated. Please wait.
      shiny::showModal(shiny::modalDialog(title = "Ploting path...",
                                          footer = NULL,
                                          easyClose = F))

      # Remove cells not in the selected partition
      pars <- monocle3::partitions(rv$cds)
      cell_coord <- SingleCellExperiment::reducedDim(rv$cds,'UMAP')
      cell_coord <- cell_coord[pars %in% input$partition_select,]
      path_cells <- rv$all_path_cells[[input$vis_path_select]]
      path_nodes <- rv$all_paths[[input$vis_path_select]]

      # Get a simplified trajectory tree
      # traj_simp <- dynwrap::simplify_trajectory(rv$traj)
      # traj_simp$directed <- rv$traj$directed

      # Get coordinates of trajectory projections
      # traj_proj <- dynwrap::project_trajectory(
      #  trajectory = traj_simp,
      #  dimred = traj_simp$dimred
      # )

      # Get coordinates of milestones
      # ms_plot_data <- as.data.frame(traj_proj$dimred_milestones)

      # Get groupings of the trajectory edges
      idx <- grep("BEGIN|END",rownames(rv$traj_proj$dimred_segment_points))
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
            rownames(rv$traj_proj$dimred_segment_points)[idx[i]]
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
        rv$traj_proj$dimred_segment_points,
        edge_group = edge_group,
        path_group = path_group
      ) %>%
        as_tibble()

      # Get cell coordinates and cells mapped to currently selected path
      cell_plot_data <- as_tibble(rv$traj_simp$dimred, rownames = NA)
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
          data = rv$ms_plot_data,
          size = 3,
          colour = "#9C9C9C"
        ) +
        ggplot2::geom_point( # Plot milestone borders
          ggplot2::aes(
            x = comp_1,
            y = comp_2,
          ),
          data = rv$ms_plot_data,
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

        # Plot in plot area
        output$path_plot <- shiny::renderPlot({
          plt
        })

        # Remove the plotting reminder window
        shiny::removeModal()
      }
    })

    # Here user specifies the temporary output folder
    shiny::observeEvent(input$tmp_folder,{
      rv$outdir_name <- shinyFiles::parseDirPath(roots = shinyFiles::getVolumes(), input$tmp_folder)
      if(length(rv$outdir_name) == 0){
        shinyFiles::shinyDirChoose(input = input, id = 'tmp_folder', roots = shinyFiles::getVolumes()())
      }else{
        output$outdir_name <- shiny::renderPrint({
          cat(sprintf("Output folder selected: %s.", rv$outdir_name))
        })
      }
    })

    ############################################################
    # Step 5. RUN STEM program
    ############################################################
    shiny::observeEvent(input$run_stem, {
      if(length(input$run_path_select) == 0){
        shiny::showModal(shiny::modalDialog(title = "Please select at least one path for STEM",
                                            footer = modalButton("OK"),
                                            easyClose = F))
      }else if(is.null(rv$outdir_name)){
        shiny::showModal(shiny::modalDialog(title = "Please specify an output folder for STEM",
                                            footer = modalButton("OK"),
                                            easyClose = F))
      }else{

        # compute metrics as input files to stem
        shiny::withProgress(message = 'Generating input for STEM', value = 0, {
          df <- get_stem_input(
            dataset = rv$dataset,
            traj = rv$traj,
            all_paths = rv$all_paths[input$run_path_select],
            path_names = input$run_path_select,
            closest_ms = rv$closest_ms,
            metric = input$metric
          )
        })

        # Keep track of the metric method used for current STEM run.
        rv$metric = input$metric

        # save metrics as tsv files
        file_names <- paste0(names(df),"_",input$metric,".tsv")
        full_file_names <- file.path(rv$outdir_name, paste0(names(df),"_",input$metric,".tsv"))

        for(i in 1:length(df)){
          data <- do.call(cbind, df[[i]])
          data <- data.frame(gene_id = rownames(data), data)
          write.table(data,
                      full_file_names[i],
                      row.names = F,
                      col.names = T,
                      quote=F,
                      sep = "\t",
                      na = "")
        }

        # Run STEM
        # res_table <- list()
        shiny::withProgress(message = "Running STEM", value = 0, {
          shiny::incProgress(0, detail = sprintf("%g%% done",0))
          for(i in 1:length(file_names)){

            # Keep track of progress
            shiny::incProgress(1/length(file_names),
                               detail = sprintf("%g%% done",round(100*i/length(file_names),2)))
            # Run STEM and get output table
            stem_analysis(
                            input_path = paste0(rv$outdir_name,"/",file_names[i]),
                            tmp_folder = rv$outdir_name,
                            stem_path = file.path(rv$stem_folder,"stem-jt.jar"),
                            setting_path = file.path(rv$outdir_name,gsub(".tsv$","",file_names[i])),
                            setting_template_path = file.path(rv$stem_folder, "stem_setting_template"),
                            species = rv$species
                          )

            # After STEM clustering is done, remove the input file, setting file.
            unlink(file.path(rv$outdir_name,gsub(".tsv$","",file_names[i])))
            unlink(file.path(rv$outdir_name, file_names[i]))
          }

        })

        # Remove GO files
        unlink(file.path(rv$outdir_name,"go-basic.obo"))
        go_file = grep("goa.*gaf\\.gz",list.files(rv$outdir_name),perl = T,value=T)
        unlink(file.path(rv$outdir_name,go_file))

      }
    })

    ############################################################
    # Step 6. Run cluster comparison
    ############################################################
    shiny::observeEvent(input$run_comparison, {
      if(is.null(input$compare1_name) | is.null(input$compare2_name)){
        shiny::showModal(shiny::modalDialog(title = "To run comparison, both paths should be specified.",
                                            footer = modalButton("OK"),
                                            easyClose = F))
      }else if(is.null(rv$outdir_name)){
        shiny::showModal(shiny::modalDialog(title = "Please specify an output folder for STEM",
                                            footer = modalButton("OK"),
                                            easyClose = F))
      }else{

        shiny::withProgress(message = 'Generating input for STEM', value = 0, {
          df <- get_stem_input(
            dataset = rv$dataset,
            traj = rv$traj,
            all_paths = rv$all_paths[c(input$compare1_name,input$compare2_name)],
            path_names = c(input$compare1_name,input$compare2_name),
            closest_ms = rv$closest_ms,
            metric = input$metric
          )
        })

        # Keep track of the metric method used for current STEM run.
        rv$metric = input$metric

        # save metrics as tsv files
        file_names <- paste0(names(df),"_",input$metric,".tsv")
        full_file_names <- file.path(rv$outdir_name, paste0(names(df),"_",input$metric,".tsv"))

        for(i in 1:length(df)){
          data <- do.call(cbind, df[[i]])
          data <- data.frame(gene_id = rownames(data), data)
          write.table(data,
                      full_file_names[i],
                      row.names = F,
                      col.names = T,
                      quote=F,
                      sep = "\t",
                      na = "")
        }

        compare1_path = paste0(rv$outdir_name,"/",input$compare1_name,"_",rv$metric,".tsv")
        compare2_path = paste0(rv$outdir_name,"/",input$compare2_name,"_",rv$metric,".tsv")

        stem_analysis(
          input_path = NULL,
          tmp_folder = rv$outdir_name,
          stem_path = file.path(rv$stem_folder,"stem-jt.jar"),
          setting_path = file.path(rv$outdir_name,"setting"),
          setting_template_path = file.path(rv$stem_folder, "stem_setting_template"),
          species = rv$species,
          compare1_path = compare1_path,
          compare2_path = compare2_path
        )

        # After STEM clustering is done, removed the input file, setting file.
        unlink(compare1_path)
        unlink(compare2_path)
        unlink(file.path(rv$outdir_name,"setting"))

        # Remove GO files
        unlink(file.path(rv$outdir_name,"go-basic.obo"))
        go_file = grep("goa.*gaf\\.gz",list.files(rv$outdir_name),perl = T,value=T)
        unlink(file.path(rv$outdir_name,go_file))
      }
    })

    # Remove files when app is closed
    # onStop(function(){
    #  file_names = grep("path.*.tsv",list.files(rv$outdir_name),perl = T,value=T)
    #  if(length(file_names) > 0){
    #    for(i in 1:length(file_names)){
    #      unlink(file.path(rv$outdir_name, file_names[i]))
    #    }
    #  }
    #})
  }
  shiny::shinyApp(ui = ui, server = server)
}

stem_analysis<-function(
  input_path, # stem input file path (when STEM is used to run a single file)
  tmp_folder, # temporary output folder
  stem_path, # stem program path
  setting_path, # stem current setting file path
  setting_template_path, # stem setting template file path
  species, # species for GO Annotations
  compare1_path = NULL, # file path for the first set of stem clusters
  compare2_path = NULL # file path for the second set of stem clusters
){

  # Change working dir to the tmp folder
  wd = getwd()
  setwd(tmp_folder)

  # Read setting file template
  settings <- readLines(setting_template_path)

  # Specify species for GO annotations
  settings[3] <- paste0("Gene_Annotation_Source\t",species)

  # Run stem in command line mode
  # Decide what operating system current R is running on.
  osname = Sys.info()['sysname']

  # Whether to perform regular clustering, or clustering + cluster comparison?
  if(is.null(compare1_path) | is.null(compare2_path)){

    # Specify input file path
    settings[2] <- paste0("Data_File\t",input_path)

    # Write current setting file
    writeLines(settings, setting_path)

    # Perform regular clustering
    # Enclose the file paths so that the spaces would not affect.
    if(osname == "Windows"){
      cmd = paste("java", "-mx1024M", "-jar",
                  paste0('"',stem_path,'"'),
                  "-d", paste0('"',setting_path,'"'),
                  "-a")
    }else{

      # on linux/MAC we need to escape the spaces.
      stem_path = gsub(" ","\\\ ", stem_path, fixed = TRUE)
      setting_path = gsub(" ","\\\ ",setting_path, fixed = TRUE)
      cmd = paste("java", "-mx1024M", "-jar",
                  stem_path,
                  "-d", setting_path,
                  "-a")
    }
  }else{

    # Specify input file path as empty (as we are running comparison now)
    settings[2] <- ""

    # Write current setting file
    writeLines(settings, setting_path)

    # Perform clustering + cluster comparison
    # Enclose the file paths so that the spaces would not affect.
    if(osname == "Windows"){
      cmd = paste("java -mx1024M -jar",
                    paste0('"',stem_path,'"'),
                    "-d", paste0('"',setting_path,'"'),
                    "-c", paste0('"',compare1_path,'"'),
                    paste0('"',compare2_path,'"'),
                    "-a","-C")
    }else{

      # on linux/MAC we need to escape the spaces.
      stem_path = gsub(" ","\\\ ", stem_path, fixed = TRUE)
      setting_path = gsub(" ","\\\ ",setting_path, fixed = TRUE)
      compare_path1 = gsub(" ","\\\ ",compare_path1, fixed = TRUE)
      compare_path2 = gsub(" ","\\\ ",compare_path2, fixed = TRUE)
      cmd = paste("java -mx1024M -jar",
                  stem_path,
                  "-d", setting_path,
                  "-c", compare1_path,
                  compare2_path,
                  "-a","-C")
    }
  }
  system(cmd)

  # remove current setting file
  unlink(setting_path)

  # Change back the working dir
  setwd(wd)
}
