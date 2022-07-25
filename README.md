# IgIDivA

## Overview

`IgIDivA` (Immunoglobulin Intraclonal Diversification Analysis) is a purpose-built tool for the analysis of the intraclonal diversification process using high-throughput sequencing data. It is written in [shiny](https://shiny.rstudio.com/). Every step of the analysis can be performed interactively, thus not requiring any programming skills. It takes as input the output files "clonotypes_computation" and "grouped_alignment_nt" from the [`tripr`](https://bio.tools/TRIP_-_T-cell_Receptor_Immunoglobulin_Profiler) package.  
Functions for an `R` command-line use are also available.  
The UserGuide of IgIDivA can be found [here](https://rpubs.com/laura_zara/911248) and an example dataset to be used as Input of IgIDivA can be found [here](https://doi.org/10.5281/zenodo.6616046).



## Installation
The `IgIDivA` scripts can be freely downloaded [here](https://github.com/laurazara/IgIDivA).
It requires `R` (version "4.1"), which can be installed on any operating system from [CRAN](https://cran.r-project.org/). `IgIDivA` has been successfully tested in Windows 10, macOS 10.15 and Linux WSL.
All the packages that need to be installed in your `R` session are the following:

```r
install.packages("shiny")
install.packages("shinyFiles")
install.packages("fs")
install.packages("pdftools")
install.packages("purrr")
install.packages("DT")
install.packages("bslib")
install.packages("shinyhelper")
install.packages("data.table")
install.packages("stringr")
install.packages("RGenetics")
install.packages("dplyr")
install.packages("ggsci")
install.packages("tidygraph")
install.packages("ggraph")
install.packages("igraph")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("rstatix")
install.packages("shinyvalidate")

```

All the scripts need to be downloaded in the same folder, and the input should be saved in a folder called `Input`  

To run the app, open the script `app.R` in your `R` session and press the button `Run App`.   


Alternatively, `IgIDivA` can be installed using a `conda` environment. We recommend to use [Miniconda](https://docs.conda.io/en/latest/miniconda.html) to install all the dependencies. The dependencies can be found in .yml format [in the IgIDivA GitHub repository](https://github.com/laurazara/IgIDivA). The yml file and the `IgIDivA` scripts need to be stored in the working folder. After downloading all the files, a terminal should be opened and the following commands should be written:  
```r
conda env create -f IgIDivA.yml 
conda activate IgIDivA
R
install.packages(c("shinyvalidate", "RGenetics","rstatix"))
q()
Rscript app.R
```  
This will produce a url that can be copied in a web browser and will direct the user to the IgIDivA app.  




