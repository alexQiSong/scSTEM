# scSTEM 0.1.0
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Getting Started
## 1. Introduction
Single cell STEM (scSTEM) is an shiny app based R package for visualizing and clustering genes in pseudotime ordered single cell RNA-seq data. scSTEM is a GUI based tool and thus does not require any coding experience.  
## 2. Installation
### 2.1 Dependencies
scSTEM works with R version 4.1.0 or higher (https://cran.r-project.org/). We would recommend installing Rstudio to interact with R in an easy-to-use GUI (https://www.rstudio.com/). scSTEM also relies on several key R packages that need to be manually installed:
- devtools
- Monocle3
- Dynverse
- ROGUE

To install **devtools**, simply type the following code in R:
```R
install.pacakges('devtools')
```
To install **Monocle3** for R, follow the installation steps described in **Monocle 3** github repos-itory page:  https://cole-trapnell-lab.github.io/monocle3/docs/installation/

To install **dynverse** for R, follow the instructions provided here: https://dynverse.org/users/1-installation.

**ROGUE** can be installed from GitHub. Simply type the following code in R:
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
scSTEM analysis involves running dimensionality reduction, cell clustering, trajectory inference, and gene clustering. While the functionality of dimensionality reduction and cell clustering is provided here to facilitate users to choose cells of interest, they are optional for the following steps. The steps of running analysis with GUI is described below. **Note that internet connection is required as scSTEM may convert gene IDs using biomaRt which converts the gene IDs using ensembl database.**

To launch the scSTEM GUI, execute the following code in R:

```R
library(scSTEM)
run_scSTEM_GUI()
```

1. **Step 1:  Load input files.** Click each of the three buttons on the top of the GUI to select corresponding input files. Click `Load files` to load all selected files. If files were successfully loaded, there will be a line of text showing `All files successfully loaded` in the top panel. We have provided a sample dataset for downloading and reproducing the clustering results (See **3 Sample data set**)
2. **Step 2: Visualizing and clustering cells.** In the second panel, scSTEM can perform dimensionality reduction and cell clustering. However, this step is optional. This step is  only intended to assist users to  select cells of interest by 2D UMAP visualization. If this step is skipped, scSTEM will take all cells from the input expression count matrix. To visualize dimensionality reduction results, first  click  `Run  UMAP`. After UMAP  is  done,  click `Run Clustering` and then `Visualize Results`.
3. **Step 3: Infer trajectories.** In the third panel, scSTEM can perform trajectory inference. To infer trajectories for the input data, simply select cell partitions of interest from the `partition` drop-down list (generated from the second step) and then select trajectory inference method from the `Method` drop-down list. Once partition and `Method` were selected, users can click `infer trajectory` to infer trajectories.
4. **Step 4: Visualize paths.** In the fourth panel, users may visualize the paths inferred by the previous step. Selected path is highlighted by red color and cells mapped to the path is marked by yellow color.
5. **Step 5: Run STEM analysis.** In the fifth panel, users can perform scSTEM analysis. Currently there are two parameters `metric` and `path` to be specified for scSTEM. `metric` specifies the method to summarize expression data from trajectories. Users can select `mean`, `entropy_reduction`, or `change_rate` from the `metric` drop-down list. In the `select a path` drop-down list, users can specify a path of interest to run scSTEM clustering (which can be visualized in **step 4**). Users will need to provide an output folder where output files and temporary files can be saved. Click the `output folder` button to specify an output folder.
6. **Step 6: Run comparative analysis.** In the sixth panel, users can perform comparison of gene clustering results from different trajectory paths. To do this, first select the two paths to be compared in the two drop-down lists `Path name 1` and `Path name 2`, then click the `Run comparison` button. scSTEM will call STEM java program to run gene clustering and comparison of clustering results for the selected paths.

# Contact
Contact us if you have any questions:
Qi (Alex) Song: qisong@andrew.cmu.edu; sqsq3178@gmail.com
Ziv Bar-Joseph: zivbj@andrew.cmu.edu

# Copyright
©2021 Qi Song, Ziv Bar-Joseph. [Systems Biology Group at Carnegie Mellon University](http://www.sb.cs.cmu.edu/). All rights reserved.
