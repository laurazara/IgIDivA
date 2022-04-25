# IgIDivA

## Overview

`IgIDivA` (Immunoglobulin Intraclonal Diversification Analysis) is a purpose-built tool for the analysis of the intraclonal diversification process using high-throughput sequencing data. It is written in [shiny](https://shiny.rstudio.com/). Every step of the analysis can be performed interactively, thus not requiring any programming skills. It takes as input the output files "clonotypes_computation" and "grouped_alignment_nt" from the `tripr` package.
Functions for an `R` command-line use are also available.



## Installation
The `IgIDivA` scripts can be freely downloaded [here](https://github.com/laurazara/IgIDivA).
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
install.packages("shinyhelper")
install.packages("RGenetics")
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

```

All the scripts need to be downloaded in the same folder, and the input should be saved in a folder called `Input`

## Launching the app
To run the app, open the script `app.R` in your `R` session and press the button `Run App`. 
