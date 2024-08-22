# Delile analysis from scratch
# AT Updated 10/3/23

#### Setup ####





# Set up default directory
setwd("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/3_PositiveSelectionDelile/0_Delile_alone/")

# Install necessary packages if not already installed
#source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/0_IstallPackages.R")

# load packages if necessary
source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/1_LoadLibraries.R")

# Import custom functions
sourceDirectory("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/", modifiedOnly=FALSE)




#### Import raw Delile Data ####




# Load the data, there are cutoffs for min_cells and min_features. percent_mt auto calculated and added into metadata
Data_raw = Load_and_merge("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/3_PositiveSelectionDelile/0_Delile_alone/", min_cells = 10, min_features = 200)
Data_raw

# Add a dummy variable
metadata <- Data_raw@meta.data
metadata$dummy_variable = "dummy_variable"
Data_raw <- AddMetaData(object = Data_raw, metadata = metadata)

# View QC on all cells in total - graphs will be exported into PDF
Graphs_raw <- qc_visuals(Data_raw, spliton = "dummy_variable")
dev.off()

Graphs_raw[[1]] # nFeature_RNA, nCount_RNA, percent_mt
Graphs_raw[[2]] # Total cells
Graphs_raw[[3]] # Features vs counts
Graphs_raw[[4]] # Counts vs features

# View QC on all cells broken down by orig.ident, graphs will be exported to PDF
Graphs_raw_sample <- qc_visuals(Data_raw, spliton = "orig.ident")
dev.off()

Graphs_raw_sample[[1]] # nFeature_RNA, nCount_RNA, percent_mt
Graphs_raw_sample[[2]] # Cells per sample
Graphs_raw_sample[[3]] # Features vs counts
Graphs_raw_sample[[4]] # Counts vs features

# Export the pre-filtered data and clean up workspace
saveRDS(Data_raw, file = "Data_Raw.rds")
remove(Data_raw, Graphs_raw_sample, Graphs_raw, metadata)




#### Remove bad quality cells #####




# Load data, if necessary
Data_raw <- readRDS("Data_Raw.rds")

#InitialQC, I set to 10,000 cells for this round
Data_QC <- ApplyQC(Data_raw, features_cutoff = 10000, max_mt = 10, low = 1, high = 2.5)
Data_QC 

# Visuals stats after QC filtering - all cells grouped together
Graphs_QC = qc_visuals(Data_QC, spliton = "dummy_variable")
dev.off()
Graphs_QC[[1]]
Graphs_QC[[2]]
Graphs_QC[[3]]
Graphs_QC[[4]]

# Visuals stats after QC filtering - cells grouped by orig.ident
Graphs_QC_sample = qc_visuals(Data_QC, spliton = "orig.ident")
dev.off()
Graphs_QC_sample[[1]]
Graphs_QC_sample[[2]]
Graphs_QC_sample[[3]]
Graphs_QC_sample[[4]]

# Save data and cleanup workspace
saveRDS(Data_QC, file = "Data_QC.rds")
remove(Data_raw, Graphs_QC, Graphs_QC_sample, Data_QC)




#### Normalize and scale ####



# Load data, if necessary
Data_QC <- readRDS("Data_QC.rds")
DefaultAssay(Data_QC) <- "RNA"

# Log normalize and scale the data (note that you are scaling all genes, not just the most variable ones)
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
top10 <- head(VariableFeatures(Data_NoInt), 10) # globin genes
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
DimPlot(Data_NoInt, reduction = "umap")

# CONCLUSION  
# Merging the data without integration does not remove batch effect

# Save the data and clean up
# Save and clean up
saveRDS(Data_NoInt, file = "Data_NoIntegration.rds")
remove(Data_NoInt, Data_norm, plot1, plot2, pcs, top10)




#### Integration without SCT ####




# Load data if necessary
Data_norm <- readRDS("Data_norm.rds")

# split the dataset into a list of seurat objects based on orig.ident
Data.list <- SplitObject(Data_norm, split.by = "orig.ident")

# identify variable features for each dataset independently
Data.list <- lapply(X = Data.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Data.list)

# INTEGRATION
# Identify anchors using the FindIntegrationAnchors() function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
anchors <- FindIntegrationAnchors(object.list = Data.list, anchor.features = features)
Data_integrated <- IntegrateData(anchorset = anchors)

# Set default assay
DefaultAssay(Data_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
Data_integrated <- ScaleData(Data_integrated, verbose = FALSE)
Data_integrated <- RunPCA(Data_integrated, npcs = 75, verbose = FALSE)
Data_integrated <- RunUMAP(Data_integrated, reduction = "pca", dims = 1:49, n.neighbors = 30L) # 50 is better than 16
Data_integrated <- FindNeighbors(Data_integrated, reduction = "pca", dims = 1:75, k.param = 25) # 50 is better than 16
Data_integrated <- FindClusters(Data_integrated, algorithm = 3, resolution = 0.9)

# Visualization
DimPlot(Data_integrated, reduction = "umap", group.by = "orig.ident", shuffle = TRUE)
DimPlot(Data_integrated, reduction = "umap", label = TRUE, repel = FALSE)

# Save and clean up
saveRDS(Data_integrated, file = "Data_IntegrationNoSCT.rds")
remove(anchors, Data_integrated, Data_norm, Data.list, features)




#### Broad cell type identification ####




# Load data if necessary
Data_integrated <- readRDS("Data_IntegrationNoSCT.rds")
  
# Highlight a particular cluster
Idents(Data_integrated) <- "seurat_clusters"
levels(Data_integrated)
highlights = WhichCells(Data_integrated, idents = 25)
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = highlights) & NoAxes()

# First set assay and Idents
DefaultAssay(Data_integrated) <- "RNA"
Idents(Data_integrated) <- "seurat_clusters"
levels(Data_integrated)

# IDENTIFY BROAD CELL TYPES FIRST

# NULL, an identity used in the paper, referring to low quality cells
FeaturePlot(Data_integrated, reduction = "umap", features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), pt.size = 1, label = TRUE, raster = FALSE) & NoAxes()
FeaturePlot(Data_integrated, reduction = "umap", features = c("nFeature_RNA"), pt.size = 1, label = TRUE, raster = FALSE) & NoAxes()
VlnPlot(Data_integrated, features = c("nFeature_RNA"), pt.size = 0)
Null = WhichCells(Data_integrated, idents = c(0, 2, 3, 16, 20, 26))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = Null) & NoAxes()

# Progenitors, general = "Sox9", "Fabp2", "Sox2"
FeaturePlot(Data_integrated, features = c("Sox9", "Fabp2", "Sox2"), label = TRUE, slot = "data") & NoAxes()
Progenitors = WhichCells(Data_integrated, idents = c(5, 9, 12, 13, 21, 22, 25, 36))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = Progenitors) & NoAxes()

# Spinal neurons general (NOT DRG neurons, see below for those) = Tubb3, Elavl3
FeaturePlot(Data_integrated, features = c("Elavl3", "Tubb3"), label = TRUE, slot = "data") & NoAxes()
Neurons = WhichCells(Data_integrated, idents = c(7, 8, 10, 11, 17, 18, 19, 23, 24, 27, 29, 30, 31, 32, 34, 38, 40))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = Neurons) & NoAxes()

# Blood-related = Sox17, Fermt3, Klf1, Hemgn, Car2
FeaturePlot(Data_integrated, features = c("Sox17", "Fermt3", "Klf1", "Hemgn", "Car2"), label = FALSE, slot = "data") & NoAxes()
FeaturePlot(Data_integrated, features = c("Fermt3"), label = FALSE, slot = "data") & NoAxes()
FeaturePlot(Data_integrated, features = c("Car2"), label = TRUE, slot = "data") & NoAxes()
Blood = WhichCells(Data_integrated, idents = c(28, 35, 37))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = Blood) & NoAxes()

# Mesoderm = Foxc1, Foxc2, Twist1, Twist2, Meox1, Meox2, Myog
FeaturePlot(Data_integrated, features = c("Foxc1", "Foxc2", "Twist1", "Twist2", "Meox1", "Meox2", "Myog"), label = TRUE, slot = "data") & NoAxes()
Mesoderm = WhichCells(Data_integrated, idents = c(1, 6, 14, 33, 39))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = Mesoderm) & NoAxes()

# Neural crest = Sox10, Sox2
FeaturePlot(Data_integrated, features = c("Sox10", "Sox2"), label = TRUE, slot = "data") & NoAxes()
NC = WhichCells(Data_integrated, idents = c(4))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = NC) & NoAxes()

# DRG = Tubb3, Elavl3, Sox10, Tlx2, Six1
FeaturePlot(Data_integrated, features = c("Tubb3", "Elavl3", "Sox10", "Tlx2", "Six1"), label = TRUE, slot = "data") & NoAxes()
DRG = WhichCells(Data_integrated, idents = c(15))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = DRG) & NoAxes()

# Skin = Krt8
FeaturePlot(Data_integrated, features = c("Krt8"), label = FALSE, slot = "data") & NoAxes()

# Add these identities to the metadata table
metadata <- Data_integrated@meta.data
head(metadata)
metadata$Broad_Identities = NaN
metadata$Broad_Identities[metadata$seurat_clusters %in% c(0, 2, 3, 16, 20, 26)] <- "Null"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(5, 9, 12, 13, 21, 22, 36)] <- "Progenitors"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(7, 8, 10, 11, 17, 18, 19, 23, 24, 25, 27, 29, 30, 31, 32, 34, 38, 40)] <- "SC_Neurons"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(28, 35, 37)] <- "Blood"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(1, 6, 14, 33, 39)] <- "Mesoderm"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(4)] <- "NC"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(15)] <- "DRG"
unique(metadata$Broad_Identities)
Data_integrated<- AddMetaData(object = Data_integrated, metadata = metadata)

# Plot the new identityes
Idents(Data_integrated) <- "Broad_Identities"
levels(Data_integrated)
DimPlot(Data_integrated, reduction = "umap", label = TRUE, repel = TRUE)

# Save and clean up
saveRDS(Data_integrated, file = "Data_IntegrationNoSCT_Labeled.rds")
remove(Data_integrated, metadata, Blood, DRG, highlights, Mesoderm, NC, Neurons, Null, Progenitors)




#### Fine cell type identification ####




# Load data if necessary
Data_integrated <- readRDS("Data_IntegrationNoSCT_Labeled.rds")

# First set identity and assay
DefaultAssay(Data_integrated) <- "RNA"
Idents(Data_integrated) <- "Broad_Identities"
levels(Data_integrated)
Idents(Data_integrated) <- "seurat_clusters"
levels(Data_integrated)

# Highlight a particular cluster
Idents(Data_integrated) <- "seurat_clusters"
levels(Data_integrated)
highlights = WhichCells(Data_integrated, idents = "21")
DimPlot(Data_integrated,reduction = "umap", label = TRUE, repel = FALSE, label.size = 6, raster = FALSE, cells.highlight = highlights) & NoAxes()

# dI1 = Pou4f1, Lhx2, Lhx9, Barhl1, Barhl2, Atoh1
FeaturePlot(Data_integrated, slot = "data", features = c("Pou4f1", "Lhx2", "Lhx9", "Barhl1", "Barhl2", "Atoh1"), label = TRUE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Atoh1"), slot = "data", label = TRUE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Lhx2"), slot = "data", label = TRUE) & NoAxes()
FeaturePlot(delile_mariano_integrated_final, features = c("Pou4f1", "Lhx2", "Lhx9", "Barhl1", "Barhl2", "Atoh1"), label = TRUE) & NoAxes()
dI1 = WhichCells(Data_integrated, idents = c(19))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = dI1) & NoAxes()

# dI2 = Pou4f1, Foxd3, Lhx1, Lhx5 
FeaturePlot(Data_integrated, features = c("Pou4f1", "Foxd3", "Lhx1", "Lhx5"), label = TRUE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Pou4f1"), label = TRUE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Foxd3"), label = TRUE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Lhx1", "Lhx5"), label = TRUE) & NoAxes()

# dI3 = Pou4f1, Isl1, Tlx3, Prrxl1, Otp
FeaturePlot(Data_integrated, features = c("Pou4f1", "Isl1", "Tlx3", "Prrxl1", "Otp"), label = TRUE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Otp"), label = TRUE) & NoAxes()
dI3 = WhichCells(Data_integrated, idents = c(32))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = dI3) & NoAxes()

# dI4 = Pax3, Pax7, Gsx2, Gsx1, Lhx1, Lhx5, Pax8, Lbx1, Pax2, Gbx1, Bhlhe22, Ptf1a
FeaturePlot(Data_integrated, features = c("Pax3", "Pax7", "Gsx2", "Gsx1", "Lhx1", "Lhx5", "Pax8", "Lbx1", "Pax2", "Gbx1", "Bhlhe22", "Ptf1a"), label = TRUE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Lhx1", "Lhx5", "Pax8", "Lbx1", "Pax2", "Gbx1"), label = TRUE) & NoAxes()

# dI5 = Pax3, Pax7, Gsx2, Gsx1, Lmx1b, Pou4f1, Tlx3, Prrxl1, Lbx1
FeaturePlot(Data_integrated, features = c("Pax3", "Pax7", "Gsx2", "Gsx1", "Lmx1b", "Pou4f1", "Tlx3", "Prrxl1", "Lbx1"), label = TRUE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Lmx1b", "Pou4f1", "Tlx3", "Prrxl1", "Lbx1"), label = TRUE) & NoAxes()
dI5 = WhichCells(Data_integrated, idents = c(7, 11, 17, 40))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = FALSE, label.size = 6, raster = FALSE, cells.highlight = dI5) & NoAxes()

# dI6 = Lhx2, Lhx5, Lbx1, Bhlhe22, Dmrt3, Wt1
FeaturePlot(Data_integrated, features = c("Lhx2", "Lhx5", "Lbx1", "Bhlhe22", "Dmrt3", "Wt1"), label = TRUE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Lhx2", "Lhx5", "Lbx1", "Bhlhe22", "Dmrt3", "Wt1"), label = FALSE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Dmrt3", "Wt1"), label = TRUE) & NoAxes()

# V0 = Lhx1, Lhx5, Evx1, Evx2, Pitx2
FeaturePlot(Data_integrated, features = c("Lhx1", "Lhx5", "Evx1", "Evx2", "Pitx2"), label = FALSE) & NoAxes()
FeaturePlot(Data_integrated, features = c("Evx1", "Evx2"), label = TRUE, repel = FALSE) & NoAxes()
V0 = WhichCells(Data_integrated, idents = c(38))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = V0) & NoAxes()

# V1 = En1, Lhx1, Lhx5, Otp
FeaturePlot(Data_integrated, slot = "data", features = c("En1", "Lhx1", "Lhx5", "Otp")) & NoAxes()
FeaturePlot(Data_integrated, slot = "data", features = c("En1"), label = TRUE) & NoAxes()
V1 = WhichCells(Data_integrated, idents = c(23))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = V1) & NoAxes()

# V2 = Msx1, Ascl1, Foxn4, Lhx1, Lhx5, Bhlhe22, Lhx3, Vsx2, Sox14 ,Sox21, Vsx1, Tal1, Gata2, Gata3
FeaturePlot(Data_integrated, slot = "data", features = c("Msx1", "Ascl1", "Foxn4", "Lhx1", "Lhx5", "Bhlhe22", "Lhx3", "Vsx2", "Sox14","Sox21", "Vsx1", "Tal1", "Gata2", "Gata3")) & NoAxes()
FeaturePlot(Data_integrated, slot = "data", features = c("Msx1", "Lhx1", "Lhx5", "Bhlhe22", "Tal1", "Gata2", "Gata3"), label = FALSE) & NoAxes()
FeaturePlot(Data_integrated, slot = "data", features = c("Tal1", "Gata2", "Gata3"), label = TRUE) & NoAxes()
V2 = WhichCells(Data_integrated, idents = c(31))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = V2) & NoAxes()

# MN = Olig2, Isl1, Lhx3, Isl2, Mnx1, Slc10a4, Slc18a3, Aldh1a2, Arhgap36
FeaturePlot(Data_integrated, slot = "data", features = c("Olig2", "Isl1", "Isl2", "Lhx3", "Mnx1", "Slc10a4", "Slc18a3", "Aldh1a2", "Arhgap36"), label = FALSE) & NoAxes()
FeaturePlot(Data_integrated, slot = "data", features = c("Isl1", "Isl2", "Lhx3", "Mnx1", "Slc10a4", "Slc18a3", "Arhgap36"), label = FALSE) & NoAxes()
FeaturePlot(Data_integrated, slot = "data", features = c("Mnx1"), label =TRUE)  & NoAxes()
MN = WhichCells(Data_integrated, idents = c(18, 29))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = MN) & NoAxes()

# V3 = Sim1, Nkx2-2
FeaturePlot(Data_integrated, slot = "data", features = c("Sim1", "Nkx2-2"), label = TRUE) & NoAxes()
V3 = WhichCells(Data_integrated, idents = c(30))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = V3) & NoAxes()

# Add these identities to the metadata table
metadata <- Data_integrated@meta.data
head(metadata)
metadata$Specific_Identities = metadata$Broad_Identities
tail(metadata)
metadata$Specific_Identities[metadata$seurat_clusters %in% c(19)] <- "dI1"
metadata$Specific_Identities[metadata$seurat_clusters %in% c(27)] <- "dI2"
metadata$Specific_Identities[metadata$seurat_clusters %in% c(32)] <- "dI3"
metadata$Specific_Identities[metadata$seurat_clusters %in% c(8, 10, 25, 34)] <- "dI4"
metadata$Specific_Identities[metadata$seurat_clusters %in% c(7, 11, 17, 40)] <- "dI5"
metadata$Specific_Identities[metadata$seurat_clusters %in% c(24)] <- "dI6"
metadata$Specific_Identities[metadata$seurat_clusters %in% c(38)] <- "V0"
metadata$Specific_Identities[metadata$seurat_clusters %in% c(23)] <- "V1"
metadata$Specific_Identities[metadata$seurat_clusters %in% c(31)] <- "V2"
metadata$Specific_Identities[metadata$seurat_clusters %in% c(18, 29)] <- "MN"
metadata$Specific_Identities[metadata$seurat_clusters %in% c(30)] <- "V3"
Data_integrated<- AddMetaData(object = Data_integrated, metadata = metadata)

# Plot the new identityes
Idents(Data_integrated) <- "Specific_Identities"
levels(Data_integrated)
DimPlot(Data_integrated, reduction = "umap", label = TRUE, repel = TRUE)

# save and clean up
saveRDS(Data_integrated, file = "Data_IntegrationNoSCT_Labeled.rds")
remove(highlights, metadata, dI1, dI2, dI3, dI4, dI5, dI6, MN, V0, V1, V2, V3)



#### Remove NULL cells from the seurat object ####




# Load labeled data, if necessary
Data_integrated <- readRDS("Data_IntegrationNoSCT_Labeled.rds")

# set default ident
Idents(Data_integrated) <- "Specific_Identities"
levels(Data_integrated)

# Select all cells that are NOT Null
Data_integrated #59671
Cells_NotNull <- WhichCells(Data_integrated, idents = "Null", invert = TRUE)
length(Cells_NotNull) 

# Load original data again and subset from there
Data_raw <- readRDS("Data_Raw.rds")
Data_integrated <- subset(Data_raw, cells = Cells_NotNull)
Data_integrated 



#### Re-analyze again ####



# Normalize
Data_integrated <- Norm(Data_integrated, regression = NULL) # no regression, much faster

# split the dataset into a list of seurat objects based on orig.ident
Data.list <- SplitObject(Data_integrated, split.by = "orig.ident")

# identify variable features for each dataset independently
Data.list <- lapply(X = Data.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Data.list)

# INTEGRATION
# identify anchors using the FindIntegrationAnchors() function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
anchors <- FindIntegrationAnchors(object.list = Data.list, anchor.features = features)
# this command creates an 'integrated' data assay
Data_integrated <- IntegrateData(anchorset = anchors)

# Set default assay
DefaultAssay(Data_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
Data_integrated <- ScaleData(Data_integrated, verbose = FALSE)
Data_integrated <- RunPCA(Data_integrated, npcs = 75, verbose = FALSE)
Data_integrated <- RunUMAP(Data_integrated, reduction = "pca", dims = 1:49, n.neighbors = 30L)
Data_integrated <- FindNeighbors(Data_integrated, reduction = "pca", dims = 1:75, k.param = 25) 
Data_integrated <- FindClusters(Data_integrated, algorithm = 3, resolution = 0.)

# Visualization
DimPlot(Data_integrated, reduction = "umap", group.by = "orig.ident", shuffle = TRUE)
DimPlot(Data_integrated, reduction = "umap", label = TRUE, repel = FALSE)


#### Broad cell type identification ####




# Load data if necessary
# Data_integrated <- readRDS("Data_IntegrationNoSCT.rds")

# Highlight a particular cluster
Idents(Data_integrated) <- "seurat_clusters"
levels(Data_integrated)
highlights = WhichCells(Data_integrated, idents = 25)
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = highlights) & NoAxes()

# First set assay and Idents
DefaultAssay(Data_integrated) <- "RNA"
Idents(Data_integrated) <- "seurat_clusters"
levels(Data_integrated)

# IDENTIFY BROAD CELL TYPES FIRST

# NULL, an identity used in the paper, referring to low quality cells
FeaturePlot(Data_integrated, reduction = "umap", features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), pt.size = 1, label = TRUE, raster = FALSE) & NoAxes()
FeaturePlot(Data_integrated, reduction = "umap", features = c("nFeature_RNA"), pt.size = 1, label = TRUE, raster = FALSE) & NoAxes()
VlnPlot(Data_integrated, features = c("nFeature_RNA"), pt.size = 0)
Null = WhichCells(Data_integrated, idents = c(28, 36))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = Null) & NoAxes()

# Progenitors, general = Sox2
FeaturePlot(Data_integrated, features = c("Sox2"), label = TRUE, slot = "data") & NoAxes()
Progenitors = WhichCells(Data_integrated, idents = c(5, 9, 12, 13, 21, 22, 25, 36))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = Progenitors) & NoAxes()

# Spinal neurons general (NOT DRG neurons, see below for those) = Tubb3, Elavl3
FeaturePlot(Data_integrated, features = c("Elavl3", "Tubb3"), label = TRUE, slot = "data") & NoAxes()
Neurons = WhichCells(Data_integrated, idents = c(7, 8, 10, 11, 17, 18, 19, 23, 24, 27, 29, 30, 31, 32, 34, 38, 40))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = Neurons) & NoAxes()

# Blood related= Sox17, Fermt3, Klf1, Hemgn, Car2
FeaturePlot(Data_integrated, features = c("Sox17", "Fermt3", "Klf1", "Hemgn", "Car2"), label = FALSE, slot = "data") & NoAxes()
FeaturePlot(Data_integrated, features = c("Fermt3"), label = FALSE, slot = "data") & NoAxes()
FeaturePlot(Data_integrated, features = c("Car2"), label = TRUE, slot = "data") & NoAxes()
Blood = WhichCells(Data_integrated, idents = c(28, 35, 37))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = Blood) & NoAxes()

# Mesoderm = Foxc1, Foxc2, Twist1, Twist2, Meox1, Meox2, Myog
FeaturePlot(Data_integrated, features = c("Foxc1", "Foxc2", "Twist1", "Twist2", "Meox1", "Meox2", "Myog"), label = TRUE, slot = "data") & NoAxes()
Mesoderm = WhichCells(Data_integrated, idents = c(1, 6, 14, 33, 39))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = Mesoderm) & NoAxes()

# Neural crest = Sox10, Sox2
FeaturePlot(Data_integrated, features = c("Sox10", "Sox2"), label = TRUE, slot = "data") & NoAxes()
NC = WhichCells(Data_integrated, idents = c(4))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = NC) & NoAxes()

# DRG = Tubb3, Elavl3, Sox10, Tlx2, Six1
FeaturePlot(Data_integrated, features = c("Tubb3", "Elavl3", "Sox10", "Tlx2", "Six1"), label = TRUE, slot = "data") & NoAxes()
DRG = WhichCells(Data_integrated, idents = c(15))
DimPlot(Data_integrated,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = DRG) & NoAxes()

# Skin = Krt8
FeaturePlot(Data_integrated, features = c("Krt8"), label = FALSE, slot = "data") & NoAxes()

# Add these identities to the metadata table
metadata <- Data_integrated@meta.data
head(metadata)
metadata$Broad_Identities = NaN
metadata$Broad_Identities[metadata$seurat_clusters %in% c(0, 2, 3, 16, 20, 26)] <- "Null"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(5, 9, 12, 13, 21, 22, 36)] <- "Progenitors"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(7, 8, 10, 11, 17, 18, 19, 23, 24, 25, 27, 29, 30, 31, 32, 34, 38, 40)] <- "SC_Neurons"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(28, 35, 37)] <- "Blood"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(1, 6, 14, 33, 39)] <- "Mesoderm"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(4)] <- "NC"
metadata$Broad_Identities[metadata$seurat_clusters %in% c(15)] <- "DRG"
unique(metadata$Broad_Identities)
Data_integrated<- AddMetaData(object = Data_integrated, metadata = metadata)

# Plot the new identityes
Idents(Data_integrated) <- "Broad_Identities"
levels(Data_integrated)
DimPlot(Data_integrated, reduction = "umap", label = TRUE, repel = TRUE)

# Save and clean up
saveRDS(Data_integrated, file = "Data_IntegrationNoSCT_Labeled.rds")
remove(Data_integrated, metadata, Blood, DRG, highlights, Mesoderm, NC, Neurons, Null, Progenitors)




