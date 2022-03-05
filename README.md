# scSTEM 0.1.0
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
- [Getting Started](#getting-started)
  * [1. Introduction](#1-introduction)
  * [2. Installation](#2-installation)
    + [2.1 OS specific reqiurements](#21-os-specific-reqiurements)
      - [Windows](#windows)
      - [Linux](#linux)
      - [MacOS](#macos)
      - [Next steps](#next-steps)
    + [2.2 Use renv to install scSTEM.](#22-use-renv-to-install-scstem)
    + [2.3 (Skip this if 2.1 is successful) Manually install all dependencies.](#23--skip-this-if-21-is-successful--manually-install-all-dependencies)
      - [2.3.1 Install R dependencies.](#231-install-r-dependencies)
    + [2.4 (Skip this if installation is successful) Trouble shootings](#24--skip-this-if-installation-is-successful--trouble-shootings)
  * [3. Sample data sets](#3-sample-data-sets)
  * [4. Input and output files](#4-input-and-output-files)
  * [5. Analysis steps](#5-analysis-steps)
  * [6. Upgrade to new version of scSTEM](#6-upgrade-to-new-version-of-scstem)
- [STEM Java program](#stem-java-program)
- [Contact](#contact)
- [Copyright](#copyright)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>
# Getting Started
## 1. Introduction
Single cell STEM (scSTEM) is a shiny app based R package for visualizing and clustering genes in pseudotime ordered single cell RNA-seq data. scSTEM is a GUI based tool and thus does not require any coding experience.  
## 2. Installation
### 2.1 OS specific reqiurements
scSTEM may have system dependencies that differ for different OS. Different requirments for different OS are described here.
#### Windows
- Install docker from https://docs.docker.com/get-docker/. You may install scSTEM without installing docker. However, many trajectory inference methods might not be running correctly except for Monocle3.
- Instal Rtools from https://cran.r-project.org/bin/windows/Rtools/  
#### Linux
Run the following commands to install Docker and other system depedencies (tested with Ubuntu 18.04):
```bash
# Install docker
sudo apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl gnupg lsb-release
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update && sudo apt-get install -y --no-install-recommends docker-ce docker-ce-cli containerd.io

# Install system dependencies
sudo apt-get install -y --no-install-recommends \
      	imagemagick libgdal-dev libhdf5-dev libudunits2-dev \
        libgeos-dev libproj-dev libgmp3-dev

```
You may install scSTEM without installing docker. However, many trajectory inference methods might not be running correctly except for Monocle3.
#### MacOS
Install docker from https://docs.docker.com/get-docker/. You may install scSTEM without installing docker. However, many trajectory inference methods might not be running correctly except for Monocle3.

#### Next steps
Once system dependencies have been installed, you may choose to install scSTEM by method **2.2** or **2.3** or **2.4**. If you encounter any issues during installation, you may refer to **2.4 Trouble shootings** or create a new `issue` under this repository.
### 2.2 Use renv to install scSTEM.
For easy installation and reproducibility, you may use `renv` to install all dependencies and `scSTEM`. Make sure that you have 1) **R version >= 4.1.0, [download R here](https://cran.r-project.org/)**, 2) **Java, [download Java here](https://java.com/en/download/help/download_options.html)** installed and 3) **internet access**. We would recommend installing Rstudio to interact with R in an easy-to-use GUI. If R asks "Do you want to install from sources the package which needs compilation?", it is recommended to select No (or N). Below are installation steps using `renv`:
1.In R, execute the following code to install `renv`:
```R
install.packages("renv")
```
2. Create a new directory to host dependency files for scSTEM (e.g. `/home/alex/scstem`). `renv` will later install dependencies into this personal folder. Download the `renv.lock` file in this repository. `renv.lock` is the lock file which contains the dependency information.
3. Once lock file is downloaded, execute the following code to activate the environment in this personal folder.
```R
install_folder = "/home/alex/scstem/"
renv::activate(install_folder)
# "/home/alex/scstem/" is the path to the directory we just created (The direcotry where the renv.lock file should be saved.)
# You will need to replace them with your own directory and path.
```
4. Once the environment is activated. Install all necessary files into this personal folder simply by using the following code.
```
renv::restore(prompt = F)
```
5. **How to run and exit scSTEM.**
To activate the project environment and run scSTEM, simply executing the following code. R will then load all dependencies from the folder we just created.
```R
renv::activate(install_folder)
library("scSTEM")
run_scstem_GUI()
```
Where `install_folder` is the project folder we just created to install scSTEM.
After analysis is done, you may deactivate the project environment for `scSTEM`, simiply by executing the following code:
```R
renv::deactivate()
```
### 2.3 (Skip this if 2.1 is successful) Manually install all dependencies.
#### 2.3.1 Install R dependencies.
Make sure that you have 1) **R version >= 4.1.0, [download R here](https://cran.r-project.org/)**, 2) **Java, [download Java here](https://java.com/en/download/help/download_options.html)** installed and 3) **internet access**. We would recommend installing Rstudio to interact with R in an easy-to-use GUI. Execute the following code to install all R dependencies. If R asks "Do you want to install from sources the package which needs compilation?", it is recommended to select No (or N). If windows users are seeing "Rtools is required to build R packages but is not currently installed" during installation, you may instal Rtools from https://cran.r-project.org/bin/windows/Rtools/
```R
# Install devtools
install.packages('devtools')

# Install Bioconductor (for R >= 4.1.0)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13", ask = FALSE)

# Install Bioconductor dependencies
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'), ask = FALSE)

# Now, install monocle3 through the cole-trapnell-lab GitHub, execute:
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase', upgrade = "always")
devtools::install_github('cole-trapnell-lab/monocle3', upgrade = "always")

# Install biomaRt
BiocManager::install("biomaRt", ask = FALSE)

# Install ROGUE
devtools::install_github("PaulingLiu/ROGUE", upgrade = "always")

# Install Dynverse
devtools::install_github("dynverse/dyno", upgrade = "always")

# Install scSTEM
devtools::install_github("alexQiSong/scSTEM", upgrade = "always")
```
### 2.4 (Skip this if installation is successful) Trouble shootings
1. `make: gfortran: No such file or directory`. MAC users may see the following error message:
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
2. `API rate limit exceeded`. This is because the dependent packages for `scSTEM` use `GitHub` API, which has been limited 60 requests by default. To increase the limit, you can execute the following code in R:
```R
usethis::create_github_token()
```
This will take you to your `GitHub` account page (you may register a new account if you haven't done so), where you can generate a new token. Once this new token has been created. You may run the following code in R:
```R
usethis::edit_r_environ()
```
This will take you to the R environment variable editor, where you can specify your token by `GITHUB_PAT = 'your_token'`. Replace `your_token` by the token you just created then close restart R and resume the installation steps. 
3. If you encounter error message related to Xcode, open terminal and type xcode-select --install and then follow the prompts. This is required by Monocle 3
## 3. Sample data sets
We have provided a sample data set for running the analysis. You may click `Load sample data` in the GUI to automatically load the sample data set. Alternatively, you may also download the sample data set (https://github.com/alexQiSong/scSTEM_sample_data) and load it manually.
## 4. Input and output files
scSTEM uses three files as input (See also the sample data set in **3. Sample data sets**):
- a **expression count matrix in \*.mtx format** with rows representing genes andcolumns representing cells.
- a **cell meta data table file in \*.csv format**, in which each row represents a cell. The  file  should  contain  at  least  two  columns  named  ’cell_id’  and  ’time_point’,where column ’timepoint’ contains time point information for each cell and should only include numeric data. Cells will be sorted by column ’time_point’ by ascending order and time point at the top is considered as the earliest time point. Cells and milestone nodes mapped to that time point are considered as starting cells and root node.
- a **gene meta data table file in \*.csv format**, in which each row represents a gene.The file should contain at least a column named ’geneid’.

scSTEM output files:
- cluster gene table (table can be saved inside STEM java GUI)
- Enriched GO term table (table can be saved inside STEM java GUI)
## 5. Analysis steps
scSTEM analysis involves running dimensionality reduction, cell clustering, trajectory inference, and gene clustering. While the functionality of dimensionality reduction and cell clustering is provided here to facilitate users to choose cells of interest, they are optional for the following steps. The steps of running analysis with GUI is described below. We will use the sample data set at https://github.com/alexQiSong/scSTEM_sample_data as an example showing the analysis functionality. This data set contains human fetal immune blood cells and was published with [Cao et al., 2020](https://www.science.org/doi/10.1126/science.aba7721). **Note that internet connection is required as scSTEM may convert gene IDs using biomaRt which converts the gene IDs using ensembl database.**

To launch the scSTEM GUI, execute the following code in R:

```R
library(scSTEM)
run_scstem_GUI()
```

1. **Step 1:  Load input files.** Click each of the three buttons on the top of the GUI to select corresponding input files. Click `Load files` to load all selected files. If files were successfully loaded, there will be a line of text showing `All files successfully loaded` in the top panel. Once the sample data set is downloaded, select the path for each input file and select the gene ID type and species. For ID type, gene symbol and ensmebl ID are supported and for species, human and mouse are supported.
 
![Alt text](img/step1.png?raw=true "Load input files")  

2. **Step 2: Visualizing and clustering cells.** In the second panel, scSTEM can perform dimensionality reduction and cell clustering. However, this step is optional. This step is only intended to assist users to  select cells of interest by 2D UMAP visualization. If this step is skipped, scSTEM will take all cells from the input expression count matrix. To visualize dimensionality reduction results, first  click  `Run  UMAP`. After UMAP  is  done,  click `Run Clustering` and then `Visualize Results`. For the sample data set, we can see the cell partition number is displayed along with the cells. Let's focus on partition 4 highlited by green color here. For the sample data set, partitions may bi numbered differently on different OS. The UMAP result will be identical across different platforms. You may find the partition shown in the figure below numbered by a different number.  

![Alt text](img/UMAP1.png?raw=true "umap1")

3. **Step 3: Infer trajectories.** In the third panel, scSTEM can perform trajectory inference. To infer trajectories for the input data, simply select cell partitions of interest from the `partition` drop-down list (generated from the second step) and then select trajectory inference method from the `Method` drop-down list. Once partition and `Method` were selected, users can click `infer trajectory` to infer trajectories. For this sample data set, let's select partition 4 and `monocle3` from the `Method` list to infer trajectoryies. The `use partition` option is only valid for `monocle3` as inference method. Once selected, monocle3 will infer disjoint trajectory graph separately for each cell partition.

![Alt text](img/step3.png?raw=true "step3")

4. **Step 4: Visualize paths.** In the fourth panel, users may visualize the paths inferred by the previous step. Selected path is highlighted by red color and cells mapped to the path is marked by yellow color. Let's select `path1` to visualize. Selected path will be shown on the right part of the GUI.  

![Alt text](img/UMAP2.png?raw=true "umap2")  

5. **Step 5: Run STEM analysis.** In the fifth panel, users can perform scSTEM analysis. Currently there are two parameters `metric` and `path` to be specified for scSTEM. `metric` specifies the method to summarize expression data from trajectories. Users can select `mean`, `entropy_reduction`, or `change_rate` from the `metric` drop-down list. In the `select a path` drop-down list, users can specify a path of interest to run scSTEM clustering (which can be visualized in **step 4**). Users will need to provide an output folder where output files and temporary files can be saved. Click the `output folder` button to specify an output folder. In this sample data set, let's select `path2` and `path3` in the `select a path` drop-down list, and `mean` in the `metric` drop-down list to run the STEM analysis.
![Alt text](img/step5.png?raw=true "umap2")  
 scSTEM will iterate through the selected paths and generate clusters for each of the path. This is done by calling the STEM java program. The STEM java GUI will show the profile plots after clustering is finished. In the screenshot shown below, **profiles with the same color together represent a gene cluster.**
![Alt text](img/stem_plot1.png?raw=true "stem_plot1")  
Users may further click the profile to have a detailed view of gene expression pattern in the pop-up window, in which users may also view and save the cluster gene and GO annotations by clicking `Cluster Gene Table` and `Cluster GO Table`, respectively.

![Alt text](img/stem_plot2.png?raw=true "stem_plot2")   

6. **Step 6: Run comparative analysis.** In the sixth panel, users can perform comparison of gene clustering results from different trajectory paths. To do this, first select the two paths to be compared in the two drop-down lists `Path name 1` and `Path name 2`, then click the `Run comparison` button. scSTEM will call STEM java program to run gene clustering and comparison of clustering results for the selected paths. For this sample data set, let's select `path2` and `path3` from the drop-down lists to compare. After clustering is done, STEM java program will show a cluster comparison plot with each row representing the matched clusters from the two paths.  

![Alt text](img/comparison.png?raw=true "comparison")  

To zoom-in the plot, you may hold and drag the right mouse. To move the plot in the window, you may hold and drag the left mouse.

## 6. Upgrade to new version of scSTEM
Upgrade existing version of scSTEM is quite easy with a single line of code in R:
```R
devtools::install_github("alexQiSong/scSTEM", upgrade = "always")
```
Note that you will need to activate the corresponding environment before upgrading scSTEM if you installed it by `renv` (**See 2.1, 5**):
```R
renv::activate(install_folder)
```

# STEM Java program
We have made changes to original STEM Java prgram as in https://github.com/jernst98/STEM_DREM. The modified version of STEM, along with its source code, can be downloaded from: https://github.com/alexQiSong/scSTEM-STEM-java 

# Contact
Contact us if you have any questions:  
Qi (Alex) Song: qisong@andrew.cmu.edu; sqsq3178@gmail.com  
Ziv Bar-Joseph: zivbj@andrew.cmu.edu  

# Copyright
©2021 Qi Song, Ziv Bar-Joseph. [Systems Biology Group at Carnegie Mellon University](http://www.sb.cs.cmu.edu/). All rights reserved.
