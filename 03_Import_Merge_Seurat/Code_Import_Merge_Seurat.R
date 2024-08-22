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




#### 0_All_Bikoff_Cells (pre-filtering) ####




# Load the data, there are cutoffs for min_cells and min_features
# This may take a while if you have a lot of cellranger files
# strained counts is the soupx data
Data_raw = Load_and_merge("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/01_Cellranger/cellranger_forcecells/", min_cells = 10, min_features = 1000, data_type = "strainedCounts")
Data_raw

# Add Bikoff metadata. This includes, age, region, replicate, and a dummy variable that is the same for all cells that can be used when you want to pool all of them together
head(Data_raw@meta.data)
Data_raw = Add_metadata(Data_raw)
head(Data_raw@meta.data)

# View QC on all cells in total - graphs will be exported into PDF
# Note this function creates a list of different QC graphs
Graphs_raw_overall <- qc_visuals(Data_raw, spliton = "dummy_variable", filename = '0_Prefiltering/Graphs_Raw_Overall.pdf')

Graphs_raw_overall[[1]] # nFeature_RNA, nCount_RNA, percent_mt
Graphs_raw_overall[[2]] # Total cells
Graphs_raw_overall[[3]] # log(Features) vs log(counts)
Graphs_raw_overall[[4]] # Counts vs features
Graphs_raw_overall[[5]] # Novelty score histogram

# View QC on all cells broken down by orig.ident, graphs will be exported to PDF
# Note this function creates a list of different QC graphs
Graphs_raw_sample <- qc_visuals(Data_raw, spliton = "orig.ident", filename = "0_Prefiltering/Graphs_Raw_Sample.pdf")

Graphs_raw_sample[[1]] # nFeature_RNA, nCount_RNA, percent_mt
Graphs_raw_sample[[2]] # Cells per sample
Graphs_raw_sample[[3]] # Features vs counts
Graphs_raw_sample[[4]] # Counts vs features
Graphs_raw_sample[[5]] # Novelty score

# Export the pre-filtered data and clean up workspace
saveRDS(Data_raw, file = "0_Prefiltering/Data_Raw_Merged.rds")
remove(Graphs_raw_sample, Graphs_raw_overall)



