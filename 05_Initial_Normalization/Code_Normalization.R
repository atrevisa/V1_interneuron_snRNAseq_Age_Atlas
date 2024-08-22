# 2_Negative_Selection_code on soupX corrected data
# AT Updated 6/11/24

#### Setup ####




# Set up default directory
setwd("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/02_NegativeSelection/soupx/")

# Install necessary packages if not already installed
#source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/0_IstallPackages.R")

# load packages if necessary
source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/1_LoadLibraries.R")

# Import custom functions
sourceDirectory("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/", modifiedOnly=FALSE)




#### 2_Normalization ####




# Load data, if necessary
# Data_QC <- readRDS("1_QC/Data_QC_merged.rds")
Data_QC
DefaultAssay(Data_QC) <- "RNA"

# Normalize, no regression
Data_norm <- Norm(Data_QC)

# Save and clean up
# dir.create("2_Normalization")
saveRDS(Data_norm, file = "2_Normalization/Data_norm_merged.rds")
remove(Data_QC)




