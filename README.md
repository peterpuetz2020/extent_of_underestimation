## Overview
This repository allows to replicate the results obtained in "Loenenbach et al. (2024) - Time course and extent of underestimation of notifiable COVID-19 cases using data from participatory, virological, and wastewater surveillance systems in Germany, 2020-2024". The analyses were conducted using R 4.3.0 (64 bit, Windows). You can recreate the project environment by using the `renv` package (https://rstudio.github.io/renv/articles/renv).

## Structure of this repository
* `In`: This folder contains the dataset needed for the analyses. See Table 1 in the corresponding paper for a description of the single variables.
* `Scripts`: This folders contains the R script needed to replicate the results.
* `Out`: This folder contains the figures of the manuscript which can be reproduced by running the R script.
* `renv`: This folder contains files to recreate the project environment and need not be called directly, see https://rstudio.github.io/renv/articles/renv.

## How to replicate the results
If you use RStudio, open the `R` project 
`extent_of_underestimation.Rproj` in the main directory first, then all R scripts loaded into this project should run as they are. If you do not use RStudio, you have to set your working directory at the beginning of each R script (`setwd(…)`) to the directory where the folders “In” and “Out” are located.

You can run the script
* `reproduce_results.R` to obtain the figures of the manuscript and the results of the breakpoint regression.
