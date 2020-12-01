# Libs

library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)
library(readxl)
library(purrr)
library(tidyr)
library(RColorBrewer)
library(kableExtra)
library(ips)
library(ape)
library(msa)
library(phangorn)

# Needed for windows
library(pegas)
library(adegenet)


options(shiny.maxRequestSize=10*1024^2) 

# Run this before deploying: 

#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.9/bioc"))
getOption("repos")

# Other ways to point to Bioconductor:

#options(repos = c(CRAN = "http://cran.rstudio.com"))
#options(repos = c(getOption("repos"), BiocInstaller::biocinstallRepos()))

# Not needed
# library(phytools)
# library(hierfstat)
# library(apex)

#~~ Tricky installations: 

# Installing older version worked:
# devtools::install_version("spdep", version = "0.8-1", repos = "http://cran.rstudio.com")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install(c("msa"))

# Source functions

source("scripts/1.0_LoadAlignment.R")
source("scripts/1.1_seqsimilarityID.R")
source("scripts/1.2_phyloID.R")
source("scripts/2.1_individualisation.R")
source("scripts/3.1_parentage.R")
source("scripts/4.1_haplotype_network.R")
#source("scripts/5.1_marker_selection.R")
source("scripts/theme_emily.R")

