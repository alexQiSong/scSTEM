# Folder for the files
library(ggplot2)
library(Matrix)

# Specify root folder for analysis
root_folder <- "scSTEM/analysis/"

traj_types = c("linear","bifurcating","multifurcating")
de_rates = c(0.1,0.3,0.5)
bar_plot_objs <- list()
box_plot_objs <- list()

for(traj_type in traj_types){
  bar_plot_data <- list()
  box_plot_data <- list()
  for(de_rate in de_rates){
    cat("=========================================\n")
    cat(sprintf("Getting results for %s %s\n",traj_type, de_rate))
    cat("=========================================\n")
    
    # Read clustering results and simulated DE genes.
    clusters <- read.csv(file.path(root_folder, "data", "input","simulation", paste0("cluster_",traj_type,"_DE",de_rate,".csv")))
    
    # Read simulated DE genes
    DE <- read.csv(file.path(root_folder,
                                        "data",
                                        "input",
                                        "simulation",
                                        paste0("geneMeta_",traj_type,"_DE",de_rate,"_.csv")
                                        )
                              )
    
    # Read expression counts
    counts <- readMM(file.path(root_folder,
                            "data",
                            "input",
                            "simulation", 
                            paste0("counts_",traj_type,"_","DE",de_rate,"_.mtx")
                            )
                  )
    
    # Normalize the counts
    norm_expressions <- log10(counts/colSums(counts)*1e6 + 1)
    
    # Read cell info (including psedutotime infomration)
    cell_info <- read.csv(file.path(root_folder, "data", "input", "simulation", paste0("cellMeta_",traj_type,"_","DE",de_rate,"_.csv")))
    
    # Find common genes between DEGs and cluster genes
    cluster_genes <- unique(clusters$gene_id)
    DE_genes <- unique(DE$gene_id[DE$DEG])
    cluster_uniq <- setdiff(cluster_genes, DE_genes)
    DE_uniq <- setdiff(DE_genes, cluster_genes)
    comm <- intersect(cluster_genes, DE_genes)
    
    # Make plot data for bar plot
    bar_plot_data[[length(bar_plot_data)+1]] <- data.frame(Category=c("scSTEM_specific","DEG_specific","Common"),
                            Number=c(length(cluster_uniq), length(DE_uniq), length(comm)),
                            DEG_percentage = scales::label_percent()(de_rate))
    
    sub_expressions1 <- as.matrix(norm_expressions[,DE$gene_id %in% comm])
    sub_expressions2 <- as.matrix(norm_expressions[,!(DE$gene_id %in% DE_uniq)])
    
    # Compute correlations between the common gene expression and pseudotime
    cor1 <- apply(sub_expressions1,2, function(x){
      select <- x!=0
      return(cor(x[select], cell_info$pseudotime[select]))
    }) %>% abs()
      
    # Compute correlations between the non-cluster gene expressions and pseudotime
    cor2 <- apply(sub_expressions2,2, function(x){
      select <- x!=0
      return(cor(x[select], cell_info$pseudotime[select]))
    }) %>% abs()
    
    # Make plot data for boxplot
    box_plot_data[[length(box_plot_data)+1]] <- data.frame(Correlation = c(cor1,cor2),
                                                           Category = c(rep("Common",length(cor1)),rep("DEG_specific",length(cor2))),
                                                           DEG_percentage = scales::label_percent()(de_rate))
      
  }
  bar_plot_data <- do.call(rbind, bar_plot_data)
  box_plot_data <- do.call(rbind, box_plot_data)
  
  # Make bar plot
  bar_plot_objs[[length(bar_plot_objs)+1]] <- ggplot(bar_plot_data, aes(fill=Category, y=Number, x=DEG_percentage)) + 
    scale_fill_manual(values = c("Common"="coral1","scSTEM_specific"="deepskyblue1","DEG_specific"="green3")) +
    geom_bar(position="stack", stat="identity") +
    geom_text(aes(label = Number),
              size = 6,
              position = position_stack(vjust = 0.5)) +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          panel.background = element_rect(fill='transparent'),
          panel.border= element_blank(),
          axis.line.x = element_line(color="black", size = 2),
          axis.line.y = element_line(color="black", size = 2),
          axis.ticks = element_line(color = "black", size = 2),
          axis.ticks.length=unit(.5, "cm")
          )

  # Make box plot
  box_plot_objs[[length(box_plot_objs)+1]] <- ggplot(box_plot_data, aes(fill=Category, y=Correlation, x=DEG_percentage)) + 
    geom_boxplot(outlier.size = 2, lwd = 1.5) +
    scale_fill_manual(values = c("Common"="coral1","scSTEM_specific"="deepskyblue1","DEG_specific"="green3")) +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          panel.background = element_rect(fill='transparent'),
          panel.border= element_blank(),
          axis.line.x = element_line(color="black", size = 2),
          axis.line.y = element_line(color="black", size = 2),
          axis.ticks = element_line(color = "black", size = 2),
          axis.ticks.length=unit(.5, "cm")
    )
  
  # Hide legends
  bar_plot_objs[[length(bar_plot_objs)]] = bar_plot_objs[[length(bar_plot_objs)]] + theme(legend.position = "none")
  box_plot_objs[[length(box_plot_objs)]] = box_plot_objs[[length(box_plot_objs)]] + theme(legend.position = "none")
}

p <- do.call(gridExtra::grid.arrange,
              list(grobs = append(bar_plot_objs,box_plot_objs),
                ncol=3,   
                nrow=2,
                top=grid::textGrob("DEGs VS scSTEM genes",gp=grid::gpar(fontsize=27)
                )
            )
        )
ggsave(file.path(root_folder,"results","sim_results.svg"),
       p,
       width = 20,
       height = 10,
       dpi = 300)
