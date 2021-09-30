# scSTEM 0.1.0
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Getting Started
## 1. Introduction
Single cell STEM (scSTEM) is an shiny app based R package for visualizing and clustering genes in pseudotime ordered single cell RNA-seq data. scSTEM is a GUI based tool and thus does not require any coding experience.  
## 2. Installation
### 2.1 Dependencies
scSTEM works with R version 4.1.0 or higher (https://cran.r-project.org/). We would recommend installing Rstudio to interact with R in an easy-to-use GUI (https://www.rstudio.com/).scSTEM also relies on several key R packages that need to be manually installed:
- devtools
- Monocle3
- Dynverse
- ROGUE

To install **devtools**, simply type the following command in R:
```R
install.pacakges('devtools')
```
To install **Monocle3** for R, follow the installation steps described in **Monocle 3** github repos-itory page:  https://cole-trapnell-lab.github.io/monocle3/docs/installation/

To install **dynverse** for R, follow the instructions provided here: https://dynverse.org/users/1-installation.

**ROGUE** can be installed from GitHub. Simply type the following command in R:
```R
devtools::install_github("PaulingLiu/ROGUE")
```
See ROGUE GitHub repository for more details:https://github.com/PaulingLiu/ROGUE.
### 2.2 Install scSTEM
After required dependencies have been successfully installed.  scSTEM can be easily installed through github by executing the following code in R:
```R
devtools::install_github("alexQiSong/scSTEM")
```
## 3. Sample data sets
We have provided a sample data set for running the analysis.  You may download the sampledata set from:  https://github.com/alexQiSong/scSTEMsampledata.
## 4. Input and output files
scSTEM uses three files as input (See also the sample data set in **3. Sample data sets**):
- a **expression count matrix in \*.mtx format** with rows representing genes andcolumns representing cells.
- a **cell meta data table file in \*.csv format**, in which each row represents a cell. The  file  should  contain  at  least  two  columns  named  ’cell_id’  and  ’time_point’,where column ’timepoint’ contains time point information for each cell and should only include numeric data. Cells will be sorted by column ’time_point’ by ascending order and time point at the top is considered as the earliest time point. Cells and milestone nodes mapped to that time point are considered as starting cells and root node.
- a **gene meta data table file in \*.csv format**, in which each row represents a gene.The file should contain at least a column named ’geneid’.

scSTEM output files:
- a table file (tentative)
- clustering plots (tentative)
## 5. Analysis steps
