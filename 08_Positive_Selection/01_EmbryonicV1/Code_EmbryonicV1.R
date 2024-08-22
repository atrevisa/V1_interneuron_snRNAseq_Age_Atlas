# 3_PositiveSelectionDelile
# 2_EmbryonicBikoffOnly
# AT Updated 1/2/24

#### Setup ####




# Set up default directory
setwd("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/03_PositiveSelectionDelile/1_Bikoff_embryonic_alone/")

# Install necessary packages if not already installed
#source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/0_IstallPackages.R")

# load packages if necessary
source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/1_LoadLibraries.R")

# Import custom functions
sourceDirectory("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/", modifiedOnly=FALSE)





#### Process Bikoff embryonic data ####




# Load in the Bikoff embryonic data
Data_raw = Load_and_merge("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/03_PositiveSelectionDelile/1_Bikoff_embryonic_alone/", min_cells = 20, min_features = 1000)
Data_raw 

# Add metadata
colnames(Data_raw@meta.data)
Data_raw <- Add_metadata(Data_raw)
colnames(Data_raw@meta.data)

# View QC on all cells in total - graphs will be exported into PDF
# Note this function creates a list of different QC graphs
Graphs_raw <- qc_visuals(Data_raw, spliton = "dummy_variable")

Graphs_raw[[1]] # nFeature_RNA, nCount_RNA, percent_mt
Graphs_raw[[2]] # Total cells
Graphs_raw[[3]] # Features vs counts
Graphs_raw[[4]] # Counts vs features

# View QC on all cells broken down by orig.ident, graphs will be exported to PDF
Graphs_raw_sample <- qc_visuals(Data_raw, spliton = "orig.ident")
Graphs_raw_sample[[1]] # nFeature_RNA, nCount_RNA, percent_mt
Graphs_raw_sample[[2]] # Cells per sample
Graphs_raw_sample[[3]] # Features vs counts
Graphs_raw_sample[[4]] # Counts vs features

# Export the pre-filtered data and clean up workspace
saveRDS(Data_raw, file = "Data_Raw.rds")
remove(Data_raw, Graphs_raw_sample, Graphs_raw)




#### Remove bad quality cells #####





#InitialQC, I set to 10,000 cells for this round
Data_QC <- ApplyQC(Data_raw, features_cutoff = 10000, max_mt = 2.5, low = 2.5, high = 2.5)
Data_QC 

# Visuals stats after QC filtering - all cells grouped together
Graphs_QC = qc_visuals(Data_QC, spliton = "dummy_variable")
# dev.off()
Graphs_QC[[1]]
Graphs_QC[[2]]
Graphs_QC[[3]]
Graphs_QC[[4]]

# Visuals stats after QC filtering - cells grouped by orig.ident
Graphs_QC_sample = qc_visuals(Data_QC, spliton = "orig.ident")
# dev.off()
Graphs_QC_sample[[1]]
Graphs_QC_sample[[2]]
Graphs_QC_sample[[3]]
Graphs_QC_sample[[4]]

# Save data and cleanup workspace
saveRDS(Data_QC, file = "Data_QC.rds")
remove(Data_raw, Graphs_QC, Graphs_QC_sample, Data_QC)




#### Normalize and scale ####





# Normalize
Data_norm <- Norm(Data_QC, regression = NULL) # no regression, much faster

# Save and clean up
saveRDS(Data_norm, file = "Data_norm.rds")
remove(Data_QC, Data_norm)




#### Dimensionality reduction without integration ####




# Load data if necessary
Data_norm <- readRDS("Data_norm.rds")

# FIND VARIABLE FEATURES AND PLOT THEM
# Find variable features
Data_NoInt <- FindVariableFeatures(Data_norm)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Data_NoInt), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Data_NoInt)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# LINEAR DIMENIONALITY REDUCTION
# Run PCA
Data_NoInt <- RunPCA(Data_NoInt, features = VariableFeatures(object = Data_NoInt))
# View the genes contributing to the top PCAs
print(Data_NoInt[["pca"]], dims = 1:5, nfeatures = 5)
# View the top contributing genes for the first two PCAs
VizDimLoadings(Data_NoInt, dims = 1:2, reduction = "pca")
# View PCA plot
DimPlot(Data_NoInt, reduction = "pca")

# Choose PCs to use for downstream analysis
pcs <- find_PC(Data_NoInt)

# NON LINEAR DIMENSIONALITY REDUCTION
# Run UMAP
Data_NoInt <- RunUMAP(Data_NoInt, dims = 1:pcs)
# visualize the UMAP
DimPlot(Data_NoInt, reduction = "umap") & NoAxes()

# Save the data and clean up
# Save and clean up
saveRDS(Data_NoInt, file = "Data_NoIntegration.rds")
remove(Data_NoInt, Data_norm, plot1, plot2, pcs, top10)




#### Visualze GOI ####




# Load data, if necessary
Data_NoInt <- readRDS("Data_NoIntegration.rds")

# QC numbers
FeaturePlot(Data_NoInt, reduction = "umap", features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), pt.size = 1, label = FALSE, raster = FALSE) & NoAxes()

# En1
FeaturePlot(Data_NoInt, features = c("En1"), label = FALSE, slot = "data") & NoAxes()

# Delile V1 markers
FeaturePlot(Data_NoInt, slot = "data", features = c("En1", "Lhx1", "Lhx5", "Otp")) & NoAxes()

# Clade markers
FeaturePlot(Data_NoInt, features = c("Foxp2", "Mafa", "Pou6f2", "Sp8"), label = FALSE, slot = "data") & NoAxes()

# Renshaw markers
FeaturePlot(Data_NoInt, features = c("Calb1", "Mafa", "Chrna2", "Onecut2"), label = FALSE, slot = "data") & NoAxes()







