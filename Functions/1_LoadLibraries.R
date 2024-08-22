# Load libraries for single cell RNA seq analysis
# AT Updated 7/18/24

library(RColorBrewer)
library(SeuratDisk)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(sctransform)
suppressPackageStartupMessages(library(clustree))
library(reticulate)
library(biomaRt)
library(remotes)
library(future)
library(tidyverse)
library(R.utils)
library(openxlsx) 
library(ggprism)
library(scCustomize)
library(scales)
library(scater)
library(fgsea)
library("msigdbr")
library(lisi) 
library(easyGgplot2)
library(devtools)
library(ggbreak)
library(reshape2)
library(hrbrthemes)
library(tidyr)
library(readr)
library(trqwe)
library(DESeq2)
library(pheatmap)
library(apeglm)

# Increase memory usage
options(future.globals.maxSize=50000*1024^2)

# Set ggplot theme
theme_set(theme_classic())

