library(Seurat)
library(magrittr)
library(readxl)

# These two files were downloaded from https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/
exp_counts <- readRDS("../data/input/human_immune/gene_count_blood.RDS")
cell_meta <- readRDS("../data/input/human_immune/df_cell_blood.RDS")

# Marker gene file, downloaded from supplementary file of https://www.science.org/doi/10.1126/science.aba7721
diff_info <- read_xlsx("C:/Users/sqsq3/Documents/stem_analysis/tmp/human_data_immune/aba7721_TablesS1-S16_character_converted.xlsx",sheet = "Table_S7",skip = 1)

# Iterate over different number of top marker genes and produce input files for scSTEM
gene_num <- c(1000,2000,3000,4000,5000,6000)
for(n in gene_num){
  if(n < 5000){
    gene_names <- diff_info[diff_info$fold.change > 2 & diff_info$qval < 0.05, ][1:n,]$gene_id
  }else{
    gene_names <- diff_info[diff_info$fold.change > 2 & diff_info$qval < 0.05, ][1:4000,]$gene_id
    gene_names <- diff_info[diff_info$fold.change > 1.5 & diff_info$fold.change < 2 & diff_info$qval < 0.05, ][1:(n-4000),]$gene_id %>%
                  c(gene_names,.)
  }
    gene_meta <- data.frame(gene_id = gene_names)
  
  colnames(cell_meta)[colnames(cell_meta) == "Development_day"] = "time_point"
  colnames(cell_meta)[colnames(cell_meta) == "sample"] = "cell_id"
  
  Matrix::writeMM(exp_counts[gene_meta$gene_id,],
                  paste0("../data/input/human_immune/counts_",n,"genes_final.mtx"))
  
  write.csv(gene_meta,
            paste0("../data/input/human_immune/gene_meta_",n,"genes_final.csv"),
            quote = F,
            row.names = F)
  
  write.csv(cell_meta,
            paste0("../data/input/human_immune/cell_meta_",n,"genes_final.csv"),
            quote = F,
            row.names = F)
}