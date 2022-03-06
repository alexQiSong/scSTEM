#################################
# Generate all GO plots
#################################
library("magrittr")
library("ggplot2")
library("dplyr")
datasets <- c("human_data_immune",
              "mouse_data_hema",
              "mouse_data_neural_crest2"
              )

trajs <- c("monocle3",
          "slingshot",
          "paga_tree"
          )

metrics <- c("mean",
            "DS",
            "change_rate"
            )

paths <- paste0("path",1:12)
clusters <- paste0("c",0:10)

for(d in datasets){
  for(t in trajs){
    for(m in metrics){
      for(path in paths){
        in_file <- paste0("../data/cluster_go/",d,"-",t,"-",m,"-",path,".csv")
        if(file.exists(in_file)){
          for(c in clusters){
            out_file <- paste0("../results/",d,"-",t,"-",m,"-",path,"_",c,".svg")
            
            # See if the cluster exists
            go <- tryCatch({
                      read.csv(in_file, row.names = 1, check.names = F)
                    },
                    error = function(e) e)
            
            # If the cluster exists, plot the GO.
            if(is.data.frame(go)){
              if(nrow(go) > 0){
                # sort by fold change and get the top 10
                go <- go %>% arrange(desc(Fold))
                if(nrow(go) >= 10){
                  go <- go[1:10,]
                }
                
                # Format GO name texts
                go$`Category Name` <- stringr::str_wrap(go$`Category Name`, width = 35)
                
                # Replace the "<0.001" p-values with 0.001 for easier visualization 
                go$`Corrected p-value` <- as.numeric(sub("<0.001",0.001, go$`Corrected p-value`))
                
                # Generate Go barlot
                go$`Category Name` <- factor(go$`Category Name`, levels = go$`Category Name`)
                p <- ggplot(data = go, aes(x = `Category Name`, y = `#Genes Enriched`, fill = `Corrected p-value`)) +
                  geom_bar(stat = "identity") + 
                  scale_fill_continuous(name = "P-value", limits = c(0,0.05), breaks = c(0, 0.05)) +
                  theme_minimal(base_size = 22) + 
                  theme(legend.position="bottom") +
                  coord_flip() +
                  xlab(NULL) +
                  ylab(NULL)
                 ggsave(out_file, width = 7, height = 7, dpi = 300)
              }
            }
          }
        }
      }
    }
  }
}
