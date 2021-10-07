# scSTEM 0.1.0
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Getting Started
## 1. Introduction
Single cell STEM (scSTEM) is a shiny app based R package for visualizing and clustering genes in pseudotime ordered single cell RNA-seq data. scSTEM is a GUI based tool and thus does not require any coding experience.  
## 2. Installation
### 2.1 (Probably) the easy way...
1. To install scSTEM in a painless way, we recommend the users to use `renv` package to automatically install the dependencies before installing scSTEM. Additionally we will also need to install `devtools`. In R, execute the following code to install `renv` and `devtools`:
```R
install.packages("renv")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```
2. Download [this **lockfile**](https://github.com/alexQiSong/scSTEM/blob/15c6b7abf609d9f21262308fe5acfc78e426ce2a/renv.lock) which contains the dependency information.
3. Once lock file is downloaded. execute the following code to install all dependencies:
```R
renv::restore("path/to/lockfile")
# "path/to/lockfile" is the path of the downloaded lock file, e.g. /home/alex/renv.lockfile.  
```
4. After all dependencies have been installed, install singularity or docker as suggested here https://dynverse.org/users/1-installation.
5. Install `scSTEM` by executing the following code:
```R
devtools::install_github("alexQiSong/scSTEM")
```
### 2.2 (Probably) the hard way...
scSTEM works with R version 4.1.0 or higher (https://cran.r-project.org/). We would recommend installing Rstudio to interact with R in an easy-to-use GUI (https://www.rstudio.com/). scSTEM also relies on several key R packages that need to be manually installed:
- devtools
- Monocle3
- Dynverse
- ROGUE
#### 2.2.1 Install devtools
To install **devtools**, simply type the following code in R:
```R
install.packages('devtools')
```
#### 2.2.2 Install Monocle3
To install **Monocle3**, execute the following in R to install **Monocle3**. For more details, please refer to the installation steps described in **Monocle 3** github repository page:  https://cole-trapnell-lab.github.io/monocle3/docs/installation/
```R
# Install Bioconductor (for R > 4.1.0)
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(version = "3.13")

# Install Bioconductor dependencies
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))

# Now, install monocle3 through the cole-trapnell-lab GitHub, execute:
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')

# To ensure that Monocle 3 was installed correctly, start a new R session and run:
library(monocle3)
```
**Note that MAC users may see the following error messages:**
```shell
make: gfortran: No such file or directory
make: *** [cigraph/src/AMD/Source/amd.o] Error 1
ERROR: compilation failed for package 'leidenbase'
```
The solutions are described below, quoted from https://cole-trapnell-lab.github.io/monocle3/docs/installation/
> * The above error indicates that you need to install gfortran on your computer (for Mac users only). In order to do so,
> * Make sure that you have Xcode command line tools installed on your computer.
> * Remove other gfortran installations if they exist. For this, you can launch a terminal window and type "which gfortran". If you see a path returned (e.g. /usr/local/bin/gfortran) you have a previous installation of gfortran that needs to be removed.
> * Download new gfortran binaries for your operating system from here and decompress the folder (eg: gunzip gfortran-8.3-bin.tar.gz).
> * Then run, sudo tar -xvf gfortran-8.3-bin.tar -C / which will install everything in /usr/local/bin/gfortran.

#### 2.2.3 Install Dynverse
To install **Dynverse** for R, follow the instructions provided here: https://dynverse.org/users/1-installation.
#### 2.2.4 Install ROGUE
**ROGUE** can be installed from GitHub. Simply type the following code in R:
```R
devtools::install_github("PaulingLiu/ROGUE")
```
See ROGUE GitHub repository for more details:https://github.com/PaulingLiu/ROGUE.
#### 2.2.5 Install scSTEM
After required dependencies have been successfully installed.  scSTEM can be easily installed through github by executing the following code in R:
```R
devtools::install_github("alexQiSong/scSTEM")
```
## 3. Sample data sets
We have provided a sample data set for running the analysis.  You may download the sampledata set from:  https://github.com/alexQiSong/scSTEM_sample_data.
## 4. Input and output files
scSTEM uses three files as input (See also the sample data set in **3. Sample data sets**):
- a **expression count matrix in \*.mtx format** with rows representing genes andcolumns representing cells.
- a **cell meta data table file in \*.csv format**, in which each row represents a cell. The  file  should  contain  at  least  two  columns  named  ’cell_id’  and  ’time_point’,where column ’timepoint’ contains time point information for each cell and should only include numeric data. Cells will be sorted by column ’time_point’ by ascending order and time point at the top is considered as the earliest time point. Cells and milestone nodes mapped to that time point are considered as starting cells and root node.
- a **gene meta data table file in \*.csv format**, in which each row represents a gene.The file should contain at least a column named ’geneid’.

scSTEM output files:
- a table file (tentative)
- clustering plots (tentative)
## 5. Analysis steps
scSTEM analysis involves running dimensionality reduction, cell clustering, trajectory inference, and gene clustering. While the functionality of dimensionality reduction and cell clustering is provided here to facilitate users to choose cells of interest, they are optional for the following steps. The steps of running analysis with GUI is described below. We will use the sample data set at https://github.com/alexQiSong/scSTEM_sample_data as an example showing the analysis functionality. This data set contains human fetal immune blood cells and was published with [Cao et al., 2020](https://www.science.org/doi/10.1126/science.aba7721). **Note that internet connection is required as scSTEM may convert gene IDs using biomaRt which converts the gene IDs using ensembl database.**

To launch the scSTEM GUI, execute the following code in R:

```R
library(scSTEM)
run_scSTEM_GUI()
```

1. **Step 1:  Load input files.** Click each of the three buttons on the top of the GUI to select corresponding input files. Click `Load files` to load all selected files. If files were successfully loaded, there will be a line of text showing `All files successfully loaded` in the top panel. Once the sample data set is downloaded, select the path for each input file and select the gene ID type and species. For ID type, gene symbol and ensmebl ID are supported and for species, human and mouse are supported. 
![Alt text](img/step1.png?raw=true "Load input files")  
2. **Step 2: Visualizing and clustering cells.** In the second panel, scSTEM can perform dimensionality reduction and cell clustering. However, this step is optional. This step is only intended to assist users to  select cells of interest by 2D UMAP visualization. If this step is skipped, scSTEM will take all cells from the input expression count matrix. To visualize dimensionality reduction results, first  click  `Run  UMAP`. After UMAP  is  done,  click `Run Clustering` and then `Visualize Results`. For the sample data set, we can see the cell partition number is displayed along with the cells. Let's focus on partition 4 highlited by green color here.  
![Alt text](img/UMAP1.png?raw=true "umap1")
4. **Step 3: Infer trajectories.** In the third panel, scSTEM can perform trajectory inference. To infer trajectories for the input data, simply select cell partitions of interest from the `partition` drop-down list (generated from the second step) and then select trajectory inference method from the `Method` drop-down list. Once partition and `Method` were selected, users can click `infer trajectory` to infer trajectories. For this sample data set, let's select partition 4 and `monocle3` from the `Method` list to infer trajectoryies. The `use partition` option is only valid for `monocle3` as inference method. Once selected, monocle3 will infer disjoint trajectory graph separately for each cell partition.
![Alt text](img/step3.png?raw=true "step3")
6. **Step 4: Visualize paths.** In the fourth panel, users may visualize the paths inferred by the previous step. Selected path is highlighted by red color and cells mapped to the path is marked by yellow color. Let's select `path1` to visualize. Selected path will be shown on the right part of the GUI.  
![Alt text](img/UMAP2.png?raw=true "umap2")  
8. **Step 5: Run STEM analysis.** In the fifth panel, users can perform scSTEM analysis. Currently there are two parameters `metric` and `path` to be specified for scSTEM. `metric` specifies the method to summarize expression data from trajectories. Users can select `mean`, `entropy_reduction`, or `change_rate` from the `metric` drop-down list. In the `select a path` drop-down list, users can specify a path of interest to run scSTEM clustering (which can be visualized in **step 4**). Users will need to provide an output folder where output files and temporary files can be saved. Click the `output folder` button to specify an output folder. In this sample data set, let's select `path2` and `path3` in the `select a path` drop-down list, and `mean` in the `metric` drop-down list to run the STEM analysis.
![Alt text](img/step5.png?raw=true "umap2")  
 scSTEM will iterate through the selected paths and generate clusters for each of the path. This is done by calling the STEM java program. The STEM java GUI will show the profile plots after clustering is finished. In the screenshot shown below, **profiles with the same color together represent a gene cluster.**
![Alt text](img/stem_plot1.png?raw=true "stem_plot1")  
Users may further click the profile to have a detailed view of gene expression pattern in the pop-up window, in which users may also view and save the cluster gene and GO annotations by clicking `Cluster Gene Table` and `Cluster GO Table`, respectively.
![Alt text](img/stem_plot2.png?raw=true "stem_plot2")   
10. **Step 6: Run comparative analysis.** In the sixth panel, users can perform comparison of gene clustering results from different trajectory paths. To do this, first select the two paths to be compared in the two drop-down lists `Path name 1` and `Path name 2`, then click the `Run comparison` button. scSTEM will call STEM java program to run gene clustering and comparison of clustering results for the selected paths. For this sample data set, let's select `path2` and `path3` from the drop-down lists to compare. After clustering is done, STEM java program will show a cluster comparison plot with each row representing the matched clusters from the two paths.  
![Alt text](img/comparison.png?raw=true "comparison")  
# Contact
Contact us if you have any questions:  
Qi (Alex) Song: qisong@andrew.cmu.edu; sqsq3178@gmail.com  
Ziv Bar-Joseph: zivbj@andrew.cmu.edu  

# Copyright
©2021 Qi Song, Ziv Bar-Joseph. [Systems Biology Group at Carnegie Mellon University](http://www.sb.cs.cmu.edu/). All rights reserved.
