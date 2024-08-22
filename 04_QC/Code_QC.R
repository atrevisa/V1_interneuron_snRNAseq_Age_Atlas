# 2_Negative_Selection_code on soupX corrected data: QC
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




#### QC ####




# Load data, if necessary
Data_raw <- readRDS("0_Prefiltering/Data_Raw_Merged.rds")
head(Data_raw@meta.data)

# InitialQC, I set to 8000 cells for this round
Data_QC <- ApplyQC(Data_raw, features_cutoff = 8000)
Data_QC 

# Visuals stats after QC filtering - all cells grouped together
Graphs_QC_overall = qc_visuals(Data_QC, spliton = "dummy_variable", filename = "1_QC/Graphs_QC_Overall.pdf")
dev.off()
Graphs_QC_overall[[1]]
Graphs_QC_overall[[2]]
Graphs_QC_overall[[3]]
Graphs_QC_overall[[4]]
Graphs_QC_overall[[5]]

# Visuals stats after QC filtering - cells grouped by orig.ident
Graphs_QC_sample = qc_visuals(Data_QC, spliton = "orig.ident", filename = "1_QC/Graphs_QC_Sample.pdf")
dev.off()
Graphs_QC_sample[[1]]
Graphs_QC_sample[[2]]
Graphs_QC_sample[[3]]
Graphs_QC_sample[[4]]
Graphs_QC_sample[[5]]

# Save data and cleanup workspace
saveRDS(Data_QC, file = "1_QC/Data_QC_merged.rds")
remove(Data_raw, Graphs_QC_overall, Graphs_QC_sample)






