# Integrate final dataset
# AT Updated 8/8/24

#### Setup #####




# Set up default directory
setwd("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/03_PositiveSelectionDelile/3_PostnatalV1Selection/soupx/LabelTransfer/")

# Install necessary packages if not already installed
# source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/0_IstallPackages.R")

# load packages if necessary
source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/1_LoadLibraries.R")

# Import custom functions
sourceDirectory("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/", modifiedOnly=FALSE)




#### Subset V1 barcodes from original data ####




# Import raw Data
Data_raw <- readRDS("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/02_NegativeSelection/soupx/0_Prefiltering/Data_Raw_Merged.rds")
Data_raw

# Import barcodes from positive selection
Cells_V1_pos <- read.csv("Cells_pos.csv")
Cells_V1_pos <- Cells_V1_pos[["x"]]
head(Cells_V1_pos)

#subset
length(Cells_V1_pos)
Data_V1 <- subset(Data_raw, cells = Cells_V1_pos)
Data_V1 




#### Log normalize ####





# Normalize the 'traditional' way, no regression
Data_V1 <- Norm(Data_V1)




#### Integrate with SCT ####




# Note that with SCT you cannot integrate without references. would have more than 2^31-1 nonzero entries
Data_V1 <- Parallel_integrate(Data_V1,  sreferences = c(1,3,5))

# Save integrated data
saveRDS(Data_V1, file = "Data_V1_SCT_Ref.rds")

# Set default assay
DefaultAssay(Data_V1) <- "integrated"

# First run PCA
Data_V1 <- RunPCA(Data_V1, verbose = TRUE, assay = "integrated")
pcs = find_PC(Data_V1)
pcs #18

# Run UMAP n.neighbors
Data_V1 <- RunUMAP(Data_V1, dims = 1:10, reduction = "pca", assay = "integrated")

# Determine the K-nearest neighbor graph
Data_V1 <- FindNeighbors(object = Data_V1,dims = 1:pcs, verbose = TRUE)

# Determine the clusters for various resolutions                                
Data_V1 <- FindClusters(object = Data_V1, resolution = 0.24, algorithm = 1)

# Visualize UMAP
DimPlot(Data_V1, group.by = "orig.ident", reduction = "umap",  label = FALSE, repel = TRUE, shuffle = TRUE, raster = FALSE) + NoAxes() + theme(aspect.ratio=1) & NoLegend()
DimPlot(Data_V1, group.by = "age", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "replicate", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "region", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "seurat_clusters", reduction = "umap", raster = TRUE, label = TRUE) & NoAxes() + theme(aspect.ratio=1)

# Highlight cells
levels(Data_V1)
V1_selection = WhichCells(Data_V1, ident = "10")
V1_selection
V2_selection = WhichCells(Data_V1, ident = "5")
# Highlight V1's in the data
DimPlot(object = Data_V1, cells.highlight = list(V1_selection, V2_selection), cols.highlight = c("red", "blue"), cols = "gray", order = TRUE)

# Quickly view cell types of interest
DefaultAssay(Data_V1) <- "RNA"
FeaturePlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
FeaturePlot(Data_V1, features = c("Foxp2", "St18", "Calb1", "Pou6f2", "Sp8", "Piezo2", "Nr5a2", "Nos1", "Reln"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Foxp2", "St18", "Calb1", "Pou6f2", "Sp8", "Piezo2", "Nr5a2", "Rnf220"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Rnf220"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend() + theme(aspect.ratio=1)
FeaturePlot(Data_V1, features = c("Sp8"), label = FALSE, raster = TRUE, slot = "data", max.cutoff = 0.3, order = TRUE) & NoAxes()
FeaturePlot(Data_V1, features = c("Pou6f2"), label = FALSE, raster = TRUE, slot = "data", max.cutoff = NA, order = FALSE) & NoAxes()

# Integration removes scaling in the RNA assay so scale the data again
DefaultAssay(Data_V1) <- "RNA"
tail(Data_V1[["RNA"]]@counts['Foxp2',]) # look at raw counts, will be integers
tail(Data_V1[["RNA"]]@data['Foxp2',]) # notice that the log normalized data looks exactly the same as the counts right now, because we haven't done the normalization yet
tail(Data_V1[["RNA"]]@scale.data['Foxp2',]) # this will throw an error because we haven't scaled the data yet either

# Normalize the 'traditional' way, no regression
Data_V1 <- Norm(Data_V1)

# Check the data after normalizing for comparison
tail(Data_V1[["RNA"]]@counts['Foxp2',]) # look at raw counts, will be integers
tail(Data_V1[["RNA"]]@data['Foxp2',]) # this should be a decimal
tail(Data_V1[["RNA"]]@scale.data['Foxp2',]) # this will no longer throw an error and will contain negative numbers

# Save integrated data
saveRDS(Data_V1, file = "Data_V1_SCT_Ref.rds")
