# Negative selection on integrated soupX data
# AT Updated 6/11/24

#### Setup ####




# Set up default directory
setwd("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/02_NegativeSelection/soupx/4_Integration")

# Install necessary packages if not already installed
#source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/0_IstallPackages.R")

# load packages if necessary
source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/1_LoadLibraries.R")

# Import custom functions
sourceDirectory("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/", modifiedOnly=FALSE)




#### Load data ####




# Load the raw soupx data
Data_raw = Load_and_merge("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/01_Cellranger/cellranger_forcecells/", min_cells = 20, min_features = 1000, data_type = "strainedCounts")
Data_raw

# Add Bikoff metadata. This includes, age, region, replicate, and a dummy variable that is the same for all cells that can be used when you want to pool all of them together
head(Data_raw@meta.data)
Data_raw = Add_metadata(Data_raw)
head(Data_raw@meta.data)

# First load the list of barcodes of interest (V1s)
Cells_V1 <- read_csv("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/02_NegativeSelection/soupx/3_NegSelection_NoInt/Cells_PutativeV1_NegSelect.csv")
length(Cells_V1$x)

# subset the V1's identified through negative selection
Data_V1 <- subset(Data_raw, cells = Cells_V1$x)
Data_V1
remove(Data_raw, Cells_V1)




#### Log Normalize ####




# Normalize the 'traditional' way, no regression
Data_V1 <- Norm(Data_V1)




#### Integrate no SCT ####




# Split the dataset into a list of seurat objects 
Data_list <- SplitObject(Data_V1, split.by = "orig.ident")
Data_list

# Find individual variable features for each dataset independently
Data_list <- lapply(X = Data_list, FUN = function(x) {
  # be sure to set assay to RNA
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, assay ="RNA")
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Data_list)

# FindIntegrationAnchors
anchors <- FindIntegrationAnchors(object.list = Data_list, anchor.features = features)

# This command creates an 'integrated' data assay
Data_V1 <- IntegrateData(anchorset = anchors)

# set default assay
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Data_V1) <- "integrated"

# Run the standard workflow for visualization and clustering
Data_V1 <- ScaleData(Data_V1, verbose = TRUE)
Data_V1 <- RunPCA(Data_V1, npcs = 50, verbose = TRUE)

# Visualize PCA results
VizDimLoadings(Data_V1, dims = 1:2, reduction = "pca")
DimPlot(Data_V1, group.by = "age", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "replicate", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "region", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "orig.ident", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimHeatmap(Data_V1, dims = 1, cells = 500, balanced = TRUE, assay = "integrated")
DimHeatmap(Data_V1, dims = 1:9, cells = 500, balanced = TRUE)

# Pick no PC's to use
ElbowPlot(Data_V1)
pcs = find_PC(Data_V1)
pcs

# Non-linear dimensionality reduction
Data_V1 <- RunUMAP(Data_V1, reduction = "pca", dims = 1:50)

# Clustering
Data_V1 <- FindNeighbors(Data_V1, reduction = "pca", dims = 1:50) # 
Data_V1 <- FindClusters(Data_V1, algorithm = 3, resolution = 0.6)

# Visualize UMAP
DimPlot(Data_V1, group.by = "orig.ident", reduction = "umap",  label = FALSE, repel = TRUE, shuffle = TRUE, raster = TRUE) + NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "age", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "replicate", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "region", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "seurat_clusters", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)

# Save integrated data
saveRDS(Data_V1, file = "Data_V1_negSelection1.rds")
remove(anchors, Data_list, features, pcs)




#### Contamination ID ####





# Set default assay to RNA
DefaultAssay(Data_V1) <- "RNA"

# Set current ident to clusters
Idents(Data_V1) <- "seurat_clusters"
levels(Data_V1)

# Vusualize clusters again
DimPlot(Data_V1, group.by = "seurat_clusters", reduction = "umap", raster = TRUE, label = TRUE) & NoAxes() + theme(aspect.ratio=1)

# pull out metadata dataframe
metadata <- Data_V1@meta.data
head(metadata)
# make a new column called ID and fill all values with "unknown"
metadata$ID <- "Putative-V1"

# Empty GEMS
FeaturePlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
FeaturePlot(Data_V1, features = c("nFeature_RNA"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
VlnPlot(Data_V1, features = c("lognFeature_RNA", "lognCount_RNA", "percent_mt", "percent_ribo"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident") +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()
VlnPlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident", add.noise = FALSE) +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(19)] <- "Null"
dplyr::count(metadata, ID, sort = TRUE)

# Glia, general
FeaturePlot(Data_V1, features = c("Sox6"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes()
VlnPlot(Data_V1, features = c("Sox6"), pt.size = 0, slot = "data", stack = FALSE, fill.by = "ident") & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(28, 24)] <- "Glia_general"
dplyr::count(metadata, ID, sort = TRUE)

# Blood-related (endothelial cells, pericytes, etc.)
FeaturePlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', 'Fli1', 'Sox17', 'Fermt3', 'Klf1', 'Hemgn', 'Car2', "Pecam1", "Tek", "My19", "Pdgfra", "Flt4"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', 'Fli1', "Pecam1", "Tek"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Flt1"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', "Pecam1", "Tek"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Microglia
FeaturePlot(Data_V1, features = c("Aif1", "Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1", "Cd68", "Adgre1", "Mertk", "Fcer1g", "Fcrls", "Hexb", "Sall1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Sall1"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Aif1", "Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1", "Cd68", "Adgre1", "Mertk", "Fcer1g", "Fcrls", "Hexb", "Sall1"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Astrocytes
FeaturePlot(Data_V1, features = c("Sox6", "Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Glul", "Gja1", "Slc1a3", "Slc4a4", "Sox2", "Slc1a2", "S100b", "Ndrg2", "Nkx6-1", "Sox9"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Slc1a3"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Oligodendrocyte lineage cells
# overall
FeaturePlot(Data_V1, features = c("Olig2", "Pdgfra", "Mbp", "Mog", "Mag", "Bcas1", "Mobp", "Plp1", "Mpz", "Pmp22", "Prx", "Cspg4", "Gpr17", "Sox10", "Cnp", "Olig1", "Nkx2-2", "Cd9", "Zfp488", "Zfp536","Nkx6-2", "Cd82", "Mal", "Bmp4","Aspa"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
# opc
FeaturePlot(Data_V1, features = c("Pdgfra", "Cspg4", "Gpr17", "Nkx2-2", "Cd9", "Bmp4"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
# mature
FeaturePlot(Data_V1, features = c("Mbp", "Mog", "Mag", "Cnp", "Mbp", "Mobp", "Plp1", "Zfp536", "Nkx6-2", "Cd82","Mal","Aspa"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
# both
FeaturePlot(Data_V1, features = c("Sox10", "Olig2", "Bcas1", "Olig1", "Zfp488", "St18"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Mag", "Mog",  "Bcas1","Mobp", "Plp1", "Mbp", "Cspg4", "Pdgfra"), pt.size = 0, slot = "data", fill.by = "ident", stack = TRUE, flip = TRUE) & NoLegend()
VlnPlot(Data_V1, features = c("Mag",  "Mobp", "Mog","Bcas1", "Cspg4", "Pdgfra"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Excitatory / Inhibitory neuron markers
FeaturePlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2", "Sox5", "Slc17a7", "Grin1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Slc17a6"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip =TRUE) & NoLegend()
FeaturePlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)
metadata$ID[metadata$seurat_clusters %in% c(13, 19)] <- "Excitatory-Neurons"
dplyr::count(metadata, ID, sort = TRUE)

# Motor
FeaturePlot(Data_V1, features = c("Chat", "Slc18a3", "Slc5a7", "Isl1"), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Chat", "Slc18a3", "Slc5a7", "Isl1"), pt.size = 0, slot = "counts")
metadata$ID[metadata$seurat_clusters %in% c(30)] <- "MN"
dplyr::count(metadata, ID, sort = TRUE)

# CSF contacting Neurons
FeaturePlot(Data_V1, features = c("Pkd2l1", "Pkd1l2", "Myo3b"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pkd2l1", "Pkd1l2", "Myo3b"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)
metadata$ID[metadata$seurat_clusters %in% c(31)] <- "CSF-cN"
dplyr::count(metadata, ID, sort = TRUE)

# Mesoderm/ myoblast
FeaturePlot(Data_V1, features = c('Foxc1', 'Foxc2', 'Twist1', 'Twist2', 'Meox1', 'Meox2', 'Myog'), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c('Foxc1', 'Foxc2', 'Twist1', 'Twist2', 'Meox1', 'Meox2', 'Myog'), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Ependymal cells
FeaturePlot(Data_V1, features = c("Dnah12", "Spef2", "Ccdc114", "Ddo", "Cfap65", "Ak9", "Fam216b", "Zfp474", "Wdr63", "Ccdc180"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Dnah12", "Spef2", "Ccdc114", "Ddo", "Cfap65", "Ak9", "Fam216b", "Zfp474", "Wdr63", "Ccdc180"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI1
FeaturePlot(Data_V1, features = c("Pou4f1", "Lhx2", "Lhx9", "Barhl1", "Barhl2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pou4f1", "Lhx2", "Lhx9", "Barhl1", "Barhl2"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI2
FeaturePlot(Data_V1, features = c("Pou4f1", "Foxd3", "Lhx1", "Lhx5"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pou4f1", "Foxd3", "Lhx1", "Lhx5"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI3 
FeaturePlot(Data_V1, features = c("Pou4f1", "Isl1", "Tlx3", "Prrxl1", "Otp"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pou4f1", "Isl1", "Tlx3", "Prrxl1", "Otp"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI4
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Pax8", "Lbx1", "Pax2", "Gbx1", "Bhlhe22"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lhx5", "Pax8", "Lbx1", "Pax2", "Gbx1", "Bhlhe22"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI5 Gsh1/2+, Lbx1+, Ptf1a−, Tlx1/3+
FeaturePlot(Data_V1, features = c("Pax2", "Pax7", "Gsx2", "Ascl1", "Gsx1", "Lmx1b", "Pou4f1", "Tlx3", "Prrxl1", "Lbx1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pax2", "Pax7", "Gsx2", "Ascl1", "Gsx1", "Lmx1b", "Pou4f1", "Tlx3", "Prrxl1", "Lbx1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI6
FeaturePlot(Data_V1, features = c("Lhx1", "Lbx1", "Pax2", "Bhlhe22", "Dmrt3", "Wt1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lbx1", "Pax2", "Bhlhe22", "Dmrt3", "Wt1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# V0
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Evx1", "Evx2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lhx5", "Evx1", "Evx2"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# V1 
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Foxd3", "En1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lhx5", "Foxd3", "En1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# V2a
FeaturePlot(Data_V1, features = c("Bhlhe22", "Lhx3", "Vsx2", "Sox14", "Sox21"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Bhlhe22", "Lhx3", "Vsx2", "Sox14", "Sox21"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# V2b
FeaturePlot(Data_V1, features = c("Bhlhe22", "Tal1", "Gata2", "Gata3"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Bhlhe22", "Tal1", "Gata2", "Gata3"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# V3
FeaturePlot(Data_V1, features = c("Nkx2-2", "Sim1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Nkx2-2", "Sim1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# RP
FeaturePlot(Data_V1, features = c("Lmx1a", "Msx1", "Msx2", "Pax3", "Wnt1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lmx1a", "Msx1", "Msx2", "Pax3", "Wnt1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# FP
FeaturePlot(Data_V1, features = c("Nkx6-2", "Foxa2", "Ferd3l", "Arx", "Shh", "Lmx1b"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Nkx6-2", "Foxa2", "Ferd3l", "Arx", "Shh", "Lmx1b"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Non MN cholinergic
FeaturePlot(Data_V1, features = c("Uncx", "Irx2", "Slc5a7", "Lpcat2","Nwd2", "2900055J20Rik", "Pid1", "Sox6"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Uncx", "Irx2", "Slc5a7", "Lpcat2","Nwd2", "2900055J20Rik", "Pid1", "Sox6"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# Skin
FeaturePlot(Data_V1, features = c('Krt8'), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c('Krt8'), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# dILA (GABAergic)= Ptf1a, Lbx1, Pax2, Lhx1/5, Gbx1, Sall3
FeaturePlot(Data_V1, features = c("Ptf1a", "Gbx1", "Lbx1", "Pax2", "Lhx1", "Lhx5", "Sall3", "Rorb", "Pdzd2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Ptf1a", "Gbx1", "Lbx1", "Pax2", "Lhx1", "Lhx5", "Sall3", "Rorb", "Pdzd2"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(15)] <- "dILA"
dplyr::count(metadata, ID, sort = TRUE)

# dILB (Excitatory = Tlx1/3, Lbx1, Lmx1b, Prrxl1, and Brn3a)
FeaturePlot(Data_V1, features = c("Tlx1", "Tlx3", "Lbx1", "Lmx1b", "Prrxl1", "Pou4f1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Tlx1", "Tlx3", "Lbx1", "Lmx1b", "Prrxl1", "Pou4f1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# Neural crest / DRG
FeaturePlot(Data_V1, features = c("Sox10", "Tlx2", 'Six1'), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Sox10", 'Tlx2', 'Six1'), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Pan neuron
FeaturePlot(Data_V1, features = c("Snap25", "Enc1", "Uchl1", "Dlg2", "Dlg4", "Rbfox3", "Mapt", "Eno2", "Dcx", 'Tubb3', 'Elavl3', "Syp", "Snhg11"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Snap25", "Enc1", "Uchl1", "Dlg2", "Dlg4", "Rbfox3", "Mapt", "Eno2", "Dcx", 'Tubb3', 'Elavl3', "Syp", "Snhg11"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Check on known V1 subsets
FeaturePlot(Data_V1, features = c("Foxp2", "Calb1", "Pou6f2", "Sp8", "Piezo2", "Nr5a2", "Nos1", "Reln"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Add metadata back to seurat object
Data_V1 <- AddMetaData(object = Data_V1, metadata = metadata)

# Plot UMAP with new identities instead of cluster number
DimPlot(Data_V1, group.by = "ID", reduction = "umap", raster = TRUE, label = FALSE, repel = TRUE, pt.size = 8, raster.dpi = c(4096, 4096)) & NoAxes() + theme(aspect.ratio=1)




#### DE ####




# Set up the seurat object for DE
DefaultAssay(Data_V1) <- "RNA"
Idents(Data_V1) <- "seurat_clusters"
levels(Data_V1)

# Find DE
markers <- FindMarkers(object = Data_V1, ident.1 = 15, only.pos = TRUE)
head(markers, 50)

# Visualize DE
FeaturePlot(Data_V1, features = c("Hoxb3"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes()




#### Selection 1 ####




# Set ident to ID
Idents(Data_V1) <- "ID"
levels(Data_V1)

# Select cell IDs that are putative V1's
Data_V1 
Cells_V1 <- WhichCells(Data_V1, idents = "Putative-V1")
length(Cells_V1)
Data_V1 <- subset(Data_V1, cells = Cells_V1)
Data_V1 

# export list of final barcodes
write.csv(Cells_V1, file ="Cells_PutativeV1_NegSelectInt.csv", row.names=FALSE)
remove(metadata, Cells_V1, markers)

# save
saveRDS(Data_V1, file = "Data_V1_negSelection.rds")








#### Integrate no SCT again ####




# Split the dataset into a list of seurat objects 
Data_list <- SplitObject(Data_V1, split.by = "orig.ident")
Data_list

# Find individual variable features for each dataset independently
Data_list <- lapply(X = Data_list, FUN = function(x) {
  # be sure to set assay to RNA
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, assay ="RNA")
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Data_list)

# FindIntegrationAnchors takes a while to run
anchors <- FindIntegrationAnchors(object.list = Data_list, anchor.features = features)

# This command creates an 'integrated' data assay. Also takes a while to run
Data_V1 <- IntegrateData(anchorset = anchors)

# Set default assay
# Specify that we will perform downstream analysis on the corrected data note that the
# Original unmodified data still resides in the 'RNA' assay
DefaultAssay(Data_V1) <- "integrated"

# Run the standard workflow for visualization and clustering
Data_V1 <- ScaleData(Data_V1, verbose = TRUE)
Data_V1 <- RunPCA(Data_V1, npcs = 50, verbose = TRUE)

# Visualize PCA results
VizDimLoadings(Data_V1, dims = 1:2, reduction = "pca")
DimPlot(Data_V1, group.by = "age", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "replicate", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "region", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "orig.ident", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimHeatmap(Data_V1, dims = 1, cells = 500, balanced = TRUE, assay = "integrated")
DimHeatmap(Data_V1, dims = 1:9, cells = 500, balanced = TRUE)

# Pick no PC's to use
ElbowPlot(Data_V1)
pcs = find_PC(Data_V1)
pcs

# Non-linear dimensionality reduction
Data_V1 <- RunUMAP(Data_V1, reduction = "pca", dims = 1:50)

# Clustering
Data_V1 <- FindNeighbors(Data_V1, reduction = "pca", dims = 1:50) # 
Data_V1 <- FindClusters(Data_V1, algorithm = 3, resolution = 1.3)

# Visualize UMAP
DimPlot(Data_V1, group.by = "orig.ident", reduction = "umap",  label = FALSE, repel = TRUE, shuffle = TRUE, raster = TRUE) + NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "age", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "replicate", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "region", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "seurat_clusters", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)

# Save integrated data
saveRDS(Data_V1, file = "Data_V1_negSelection2.rds")
remove(anchors, Data_list, features, pcs)





#### Contamination ID ####





# Set default assay to RNA
DefaultAssay(Data_V1) <- "RNA"

# Set current ident to clusters
Idents(Data_V1) <- "seurat_clusters"
levels(Data_V1)

# Vusualize clusters again
DimPlot(Data_V1, group.by = "seurat_clusters", reduction = "umap", raster = TRUE, label = TRUE) & NoAxes() + theme(aspect.ratio=1)

# pull out metadata dataframe
metadata <- Data_V1@meta.data
head(metadata)
# make a new column called ID and fill all values with "unknown"
metadata$ID <- "Putative-V1"

# Empty GEMS
FeaturePlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
FeaturePlot(Data_V1, features = c("nFeature_RNA"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
VlnPlot(Data_V1, features = c("lognFeature_RNA", "lognCount_RNA", "percent_mt", "percent_ribo"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident") +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()
VlnPlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident", add.noise = FALSE) +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(22)] <- "Null"
dplyr::count(metadata, ID, sort = TRUE)

# Glia, general
FeaturePlot(Data_V1, features = c("Sox6"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes()
VlnPlot(Data_V1, features = c("Sox6"), pt.size = 0, slot = "data", stack = FALSE, fill.by = "ident") & NoLegend()

# Blood-related (endothelial cells, pericytes, etc.)
FeaturePlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', 'Fli1', 'Sox17', 'Fermt3', 'Klf1', 'Hemgn', 'Car2', "Pecam1", "Tek", "My19", "Pdgfra", "Flt4"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', 'Fli1', "Pecam1", "Tek"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Flt1"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', "Pecam1", "Tek"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Microglia
FeaturePlot(Data_V1, features = c("Aif1", "Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1", "Cd68", "Adgre1", "Mertk", "Fcer1g", "Fcrls", "Hexb", "Sall1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Sall1"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Aif1", "Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1", "Cd68", "Adgre1", "Mertk", "Fcer1g", "Fcrls", "Hexb", "Sall1"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Astrocytes
FeaturePlot(Data_V1, features = c("Sox6", "Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Glul", "Gja1", "Slc1a3", "Slc4a4", "Sox2", "Slc1a2", "S100b", "Ndrg2", "Nkx6-1", "Sox9"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Slc1a3"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Oligodendrocyte lineage cells
# overall
FeaturePlot(Data_V1, features = c("Olig2", "Pdgfra", "Mbp", "Mog", "Mag", "Bcas1", "Mobp", "Plp1", "Mpz", "Pmp22", "Prx", "Cspg4", "Gpr17", "Sox10", "Cnp", "Olig1", "Nkx2-2", "Cd9", "Zfp488", "Zfp536","Nkx6-2", "Cd82", "Mal", "Bmp4","Aspa"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
# opc
FeaturePlot(Data_V1, features = c("Pdgfra", "Cspg4", "Gpr17", "Nkx2-2", "Cd9", "Bmp4"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
# mature
FeaturePlot(Data_V1, features = c("Mbp", "Mog", "Mag", "Cnp", "Mbp", "Mobp", "Plp1", "Zfp536", "Nkx6-2", "Cd82","Mal","Aspa"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
# both
FeaturePlot(Data_V1, features = c("Sox10", "Olig2", "Bcas1", "Olig1", "Zfp488", "St18"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Mag", "Mog",  "Bcas1","Mobp", "Plp1", "Mbp", "Cspg4", "Pdgfra"), pt.size = 0, slot = "data", fill.by = "ident", stack = TRUE, flip = TRUE) & NoLegend()
VlnPlot(Data_V1, features = c("Mag",  "Mobp", "Mog","Bcas1", "Cspg4", "Pdgfra"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Excitatory / Inhibitory neuron markers
FeaturePlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2", "Sox5", "Slc17a7", "Grin1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Slc17a6"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip =TRUE) & NoLegend()
FeaturePlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)
metadata$ID[metadata$seurat_clusters %in% c(41)] <- "Excitatory"
dplyr::count(metadata, ID, sort = TRUE)

# Motor neurons
FeaturePlot(Data_V1, features = c("Chat", "Slc18a3", "Slc5a7", "Isl1"), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Chat", "Slc18a3", "Slc5a7", "Isl1"), pt.size = 0, slot = "counts")

# CSF contacting Neurons
FeaturePlot(Data_V1, features = c("Pkd2l1", "Pkd1l2", "Myo3b"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pkd2l1", "Pkd1l2", "Myo3b"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Mesoderm/ myoblast
FeaturePlot(Data_V1, features = c('Foxc1', 'Foxc2', 'Twist1', 'Twist2', 'Meox1', 'Meox2', 'Myog'), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c('Foxc1', 'Foxc2', 'Twist1', 'Twist2', 'Meox1', 'Meox2', 'Myog'), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Ependymal cells
FeaturePlot(Data_V1, features = c("Dnah12", "Spef2", "Ccdc114", "Ddo", "Cfap65", "Ak9", "Fam216b", "Zfp474", "Wdr63", "Ccdc180"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Dnah12", "Spef2", "Ccdc114", "Ddo", "Cfap65", "Ak9", "Fam216b", "Zfp474", "Wdr63", "Ccdc180"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI1
FeaturePlot(Data_V1, features = c("Pou4f1", "Lhx2", "Lhx9", "Barhl1", "Barhl2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pou4f1", "Lhx2", "Lhx9", "Barhl1", "Barhl2"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI2
FeaturePlot(Data_V1, features = c("Pou4f1", "Foxd3", "Lhx1", "Lhx5"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pou4f1", "Foxd3", "Lhx1", "Lhx5"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI3 
FeaturePlot(Data_V1, features = c("Pou4f1", "Isl1", "Tlx3", "Prrxl1", "Otp"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pou4f1", "Isl1", "Tlx3", "Prrxl1", "Otp"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI4
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Pax8", "Lbx1", "Pax2", "Gbx1", "Bhlhe22"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lhx5", "Pax8", "Lbx1", "Pax2", "Gbx1", "Bhlhe22"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI5 Gsh1/2+, Lbx1+, Ptf1a−, Tlx1/3+
FeaturePlot(Data_V1, features = c("Pax2", "Pax7", "Gsx2", "Ascl1", "Gsx1", "Lmx1b", "Pou4f1", "Tlx3", "Prrxl1", "Lbx1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pax2", "Pax7", "Gsx2", "Ascl1", "Gsx1", "Lmx1b", "Pou4f1", "Tlx3", "Prrxl1", "Lbx1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# dI6
FeaturePlot(Data_V1, features = c("Lhx1", "Lbx1", "Pax2", "Bhlhe22", "Dmrt3", "Wt1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lbx1", "Pax2", "Bhlhe22", "Dmrt3", "Wt1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# V0
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Evx1", "Evx2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lhx5", "Evx1", "Evx2"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# V1 
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Foxd3", "En1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lhx5", "Foxd3", "En1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# V2a
FeaturePlot(Data_V1, features = c("Bhlhe22", "Lhx3", "Vsx2", "Sox14", "Sox21"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Bhlhe22", "Lhx3", "Vsx2", "Sox14", "Sox21"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# V2b
FeaturePlot(Data_V1, features = c("Bhlhe22", "Tal1", "Gata2", "Gata3"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Bhlhe22", "Tal1", "Gata2", "Gata3"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# V3
FeaturePlot(Data_V1, features = c("Nkx2-2", "Sim1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Nkx2-2", "Sim1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# RP
FeaturePlot(Data_V1, features = c("Lmx1a", "Msx1", "Msx2", "Pax3", "Wnt1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lmx1a", "Msx1", "Msx2", "Pax3", "Wnt1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# FP
FeaturePlot(Data_V1, features = c("Nkx6-2", "Foxa2", "Ferd3l", "Arx", "Shh", "Lmx1b"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Nkx6-2", "Foxa2", "Ferd3l", "Arx", "Shh", "Lmx1b"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Non MN cholinergic
FeaturePlot(Data_V1, features = c("Uncx", "Irx2", "Slc5a7", "Lpcat2","Nwd2", "2900055J20Rik", "Pid1", "Sox6"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Uncx", "Irx2", "Slc5a7", "Lpcat2","Nwd2", "2900055J20Rik", "Pid1", "Sox6"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# Skin
FeaturePlot(Data_V1, features = c('Krt8'), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c('Krt8'), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# dILA (GABAergic)= Ptf1a, Lbx1, Pax2, Lhx1/5, Gbx1, Sall3
FeaturePlot(Data_V1, features = c("Ptf1a", "Gbx1", "Lbx1", "Pax2", "Lhx1", "Lhx5", "Sall3", "Rorb", "Pdzd2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Ptf1a", "Gbx1", "Lbx1", "Pax2", "Lhx1", "Lhx5", "Sall3", "Rorb", "Pdzd2"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(15)] <- "dILA"
dplyr::count(metadata, ID, sort = TRUE)

# dILB (Excitatory = Tlx1/3, Lbx1, Lmx1b, Prrxl1, and Brn3a)
FeaturePlot(Data_V1, features = c("Tlx1", "Tlx3", "Lbx1", "Lmx1b", "Prrxl1", "Pou4f1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Tlx1", "Tlx3", "Lbx1", "Lmx1b", "Prrxl1", "Pou4f1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# Neural crest / DRG
FeaturePlot(Data_V1, features = c("Sox10", "Tlx2", 'Six1'), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Sox10", 'Tlx2', 'Six1'), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Pan neuron
FeaturePlot(Data_V1, features = c("Snap25", "Enc1", "Uchl1", "Dlg2", "Dlg4", "Rbfox3", "Mapt", "Eno2", "Dcx", 'Tubb3', 'Elavl3', "Syp", "Snhg11"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Snap25", "Enc1", "Uchl1", "Dlg2", "Dlg4", "Rbfox3", "Mapt", "Eno2", "Dcx", 'Tubb3', 'Elavl3', "Syp", "Snhg11"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Check on known V1 subsets
FeaturePlot(Data_V1, features = c("Foxp2", "Calb1", "Pou6f2", "Sp8", "Piezo2", "Nr5a2", "Nos1", "Reln"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Add metadata back to seurat object
Data_V1 <- AddMetaData(object = Data_V1, metadata = metadata)

# Plot UMAP with new identities instead of cluster number
DimPlot(Data_V1, group.by = "ID", reduction = "umap", raster = TRUE, label = FALSE, repel = TRUE, pt.size = 8, raster.dpi = c(4096, 4096)) & NoAxes() + theme(aspect.ratio=1)





#### Selection 2 from raw data ####




# Set ident to ID
Idents(Data_V1) <- "ID"
levels(Data_V1)

# Select cell IDs that are putative V1's
Data_V1
Cells_V1 <- WhichCells(Data_V1, idents = "Putative-V1")
length(Cells_V1)
# Import raw Data
Data_raw <- readRDS("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/2_NegativeSelection/soupx/0_Prefiltering/Data_Raw_Merged.rds")
#subset
Data_V1 <- subset(Data_raw, cells = Cells_V1)
Data_V1 #89545

# export list of final barcodes
write.csv(Cells_V1, file ="Cells_PutativeV1_NegSelectInt.csv", row.names=FALSE)
remove(metadata, Cells_V1, Data_raw)




#### Log normalize again ####





# Normalize the 'traditional' way, no regression
Data_V1 <- Norm(Data_V1)




#### Integrate no SCT again ####




# Split the dataset into a list of seurat objects 
Data_list <- SplitObject(Data_V1, split.by = "orig.ident")
Data_list

# Find individual variable features for each dataset independently
Data_list <- lapply(X = Data_list, FUN = function(x) {
  # be sure to set assay to RNA
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, assay ="RNA")
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Data_list)

# FindIntegrationAnchors
anchors <- FindIntegrationAnchors(object.list = Data_list, anchor.features = features)

# This command creates an 'integrated' data assay
Data_V1 <- IntegrateData(anchorset = anchors)

# Set default assay
# Specify that we will perform downstream analysis on the corrected data note that the
# Original unmodified data still resides in the 'RNA' assay
DefaultAssay(Data_V1) <- "integrated"

# Run the standard workflow for visualization and clustering
Data_V1 <- ScaleData(Data_V1, verbose = TRUE)
Data_V1 <- RunPCA(Data_V1, npcs = 50, verbose = TRUE)

# Visualize PCA results
VizDimLoadings(Data_V1, dims = 1:2, reduction = "pca")
DimPlot(Data_V1, group.by = "age", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "replicate", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "region", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "orig.ident", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimHeatmap(Data_V1, dims = 1, cells = 500, balanced = TRUE, assay = "integrated")
DimHeatmap(Data_V1, dims = 1:9, cells = 500, balanced = TRUE)

# Pick no PC's to use
ElbowPlot(Data_V1)
pcs = find_PC(Data_V1)
pcs

# Non-linear dimensionality reduction
Data_V1 <- RunUMAP(Data_V1, reduction = "pca", dims = 1:50)

# Clustering
Data_V1 <- FindNeighbors(Data_V1, reduction = "pca", dims = 1:50) # 
Data_V1 <- FindClusters(Data_V1, algorithm = 3, resolution = 1)

# Visualize UMAP
DimPlot(Data_V1, group.by = "orig.ident", reduction = "umap",  label = FALSE, repel = TRUE, shuffle = TRUE, raster = TRUE) + NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "age", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "replicate", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "region", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "seurat_clusters", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)

# Save integrated data
saveRDS(Data_V1, file = "Data_V1_negSelection3.rds")




#### Contamination ID ####





# Set default assay to RNA
DefaultAssay(Data_V1) <- "RNA"

# Set current ident to clusters
Idents(Data_V1) <- "seurat_clusters"
levels(Data_V1)

# Vusualize clusters again
DimPlot(Data_V1, group.by = "seurat_clusters", reduction = "umap", raster = TRUE, label = TRUE) & NoAxes() + theme(aspect.ratio=1)

# pull out metadata dataframe
metadata <- Data_V1@meta.data
head(metadata)
# make a new column called ID and fill all values with "unknown"
metadata$ID <- "Putative-V1"

# Empty GEMS
FeaturePlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
FeaturePlot(Data_V1, features = c("nFeature_RNA"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
VlnPlot(Data_V1, features = c("lognFeature_RNA", "lognCount_RNA", "percent_mt", "percent_ribo"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident") +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()
VlnPlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident", add.noise = FALSE) +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()

# Glia, general
FeaturePlot(Data_V1, features = c("Sox6"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes()
VlnPlot(Data_V1, features = c("Sox6"), pt.size = 0, slot = "data", stack = FALSE, fill.by = "ident") & NoLegend()

# Blood-related (endothelial cells, pericytes, etc.)
FeaturePlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', 'Fli1', 'Sox17', 'Fermt3', 'Klf1', 'Hemgn', 'Car2', "Pecam1", "Tek", "My19", "Pdgfra", "Flt4"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', 'Fli1', "Pecam1", "Tek"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Flt1"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', "Pecam1", "Tek"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Microglia
FeaturePlot(Data_V1, features = c("Aif1", "Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1", "Cd68", "Adgre1", "Mertk", "Fcer1g", "Fcrls", "Hexb", "Sall1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Sall1"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Aif1", "Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1", "Cd68", "Adgre1", "Mertk", "Fcer1g", "Fcrls", "Hexb", "Sall1"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Astrocytes
FeaturePlot(Data_V1, features = c("Sox6", "Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Glul", "Gja1", "Slc1a3", "Slc4a4", "Sox2", "Slc1a2", "S100b", "Ndrg2", "Nkx6-1", "Sox9"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Slc1a3"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Oligodendrocyte lineage cells
# overall
FeaturePlot(Data_V1, features = c("Olig2", "Pdgfra", "Mbp", "Mog", "Mag", "Bcas1", "Mobp", "Plp1", "Mpz", "Pmp22", "Prx", "Cspg4", "Gpr17", "Sox10", "Cnp", "Olig1", "Nkx2-2", "Cd9", "Zfp488", "Zfp536","Nkx6-2", "Cd82", "Mal", "Bmp4","Aspa"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
# opc
FeaturePlot(Data_V1, features = c("Pdgfra", "Cspg4", "Gpr17", "Nkx2-2", "Cd9", "Bmp4"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
# mature
FeaturePlot(Data_V1, features = c("Mbp", "Mog", "Mag", "Cnp", "Mbp", "Mobp", "Plp1", "Zfp536", "Nkx6-2", "Cd82","Mal","Aspa"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
# both
FeaturePlot(Data_V1, features = c("Sox10", "Olig2", "Bcas1", "Olig1", "Zfp488", "St18"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Mag", "Mog",  "Bcas1","Mobp", "Plp1", "Mbp", "Cspg4", "Pdgfra"), pt.size = 0, slot = "data", fill.by = "ident", stack = TRUE, flip = TRUE) & NoLegend()
VlnPlot(Data_V1, features = c("Mag",  "Mobp", "Mog","Bcas1", "Cspg4", "Pdgfra"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Excitatory / Inhibitory neuron markers
FeaturePlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2", "Sox5", "Slc17a7", "Grin1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Slc17a6"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip =TRUE) & NoLegend()
FeaturePlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Motor
FeaturePlot(Data_V1, features = c("Chat", "Slc18a3", "Slc5a7", "Isl1"), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Chat", "Slc18a3", "Slc5a7", "Isl1"), pt.size = 0, slot = "counts", stack = TRUE, flip = TRUE) & NoLegend()

# CSF contacting Neurons. Notice below when you look at inhibitory markers that cluster 19 is also gaba-ergic, which agrees with published literature
FeaturePlot(Data_V1, features = c("Pkd2l1", "Pkd1l2", "Myo3b"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pkd2l1", "Pkd1l2", "Myo3b"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# Mesoderm/ myoblast
FeaturePlot(Data_V1, features = c('Foxc1', 'Foxc2', 'Twist1', 'Twist2', 'Meox1', 'Meox2', 'Myog'), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c('Foxc1', 'Foxc2', 'Twist1', 'Twist2', 'Meox1', 'Meox2', 'Myog'), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Ependymal cells
FeaturePlot(Data_V1, features = c("Dnah12", "Spef2", "Ccdc114", "Ddo", "Cfap65", "Ak9", "Fam216b", "Zfp474", "Wdr63", "Ccdc180"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Dnah12", "Spef2", "Ccdc114", "Ddo", "Cfap65", "Ak9", "Fam216b", "Zfp474", "Wdr63", "Ccdc180"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# dI1
FeaturePlot(Data_V1, features = c("Pou4f1", "Lhx2", "Lhx9", "Barhl1", "Barhl2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pou4f1", "Lhx2", "Lhx9", "Barhl1", "Barhl2"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# dI2
FeaturePlot(Data_V1, features = c("Pou4f1", "Foxd3", "Lhx1", "Lhx5"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pou4f1", "Foxd3", "Lhx1", "Lhx5"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# dI3 
FeaturePlot(Data_V1, features = c("Pou4f1", "Isl1", "Tlx3", "Prrxl1", "Otp"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pou4f1", "Isl1", "Tlx3", "Prrxl1", "Otp"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# dI4
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Pax8", "Lbx1", "Pax2", "Gbx1", "Bhlhe22"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lhx5", "Pax8", "Lbx1", "Pax2", "Gbx1", "Bhlhe22"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# dI5 Gsh1/2+, Lbx1+, Ptf1a−, Tlx1/3+
FeaturePlot(Data_V1, features = c("Pax2", "Pax7", "Gsx2", "Ascl1", "Gsx1", "Lmx1b", "Pou4f1", "Tlx3", "Prrxl1", "Lbx1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pax2", "Pax7", "Gsx2", "Ascl1", "Gsx1", "Lmx1b", "Pou4f1", "Tlx3", "Prrxl1", "Lbx1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# dI6
FeaturePlot(Data_V1, features = c("Lhx1", "Lbx1", "Pax2", "Bhlhe22", "Dmrt3", "Wt1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lbx1", "Pax2", "Bhlhe22", "Dmrt3", "Wt1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# V0
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Evx1", "Evx2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lhx5", "Evx1", "Evx2"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# V1 
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Foxd3", "En1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lhx1", "Lhx5", "Foxd3", "En1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# V2a
FeaturePlot(Data_V1, features = c("Bhlhe22", "Lhx3", "Vsx2", "Sox14", "Sox21"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Bhlhe22", "Lhx3", "Vsx2", "Sox14", "Sox21"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# V2b
FeaturePlot(Data_V1, features = c("Bhlhe22", "Tal1", "Gata2", "Gata3"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Bhlhe22", "Tal1", "Gata2", "Gata3"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# V3
FeaturePlot(Data_V1, features = c("Nkx2-2", "Sim1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Nkx2-2", "Sim1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# RP
FeaturePlot(Data_V1, features = c("Lmx1a", "Msx1", "Msx2", "Pax3", "Wnt1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Lmx1a", "Msx1", "Msx2", "Pax3", "Wnt1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# FP
FeaturePlot(Data_V1, features = c("Nkx6-2", "Foxa2", "Ferd3l", "Arx", "Shh", "Lmx1b"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Nkx6-2", "Foxa2", "Ferd3l", "Arx", "Shh", "Lmx1b"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# Non MN cholinergic
FeaturePlot(Data_V1, features = c("Uncx", "Irx2", "Slc5a7", "Lpcat2","Nwd2", "2900055J20Rik", "Pid1", "Sox6"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Uncx", "Irx2", "Slc5a7", "Lpcat2","Nwd2", "2900055J20Rik", "Pid1", "Sox6"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# Skin
FeaturePlot(Data_V1, features = c('Krt8'), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c('Krt8'), pt.size = 0, slot = "data", stack = FALSE, flip = TRUE) & NoLegend()

# dILA (GABAergic)= Ptf1a, Lbx1, Pax2, Lhx1/5, Gbx1, Sall3
FeaturePlot(Data_V1, features = c("Ptf1a", "Gbx1", "Lbx1", "Pax2", "Lhx1", "Lhx5", "Sall3", "Rorb", "Pdzd2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Ptf1a", "Gbx1", "Lbx1", "Pax2", "Lhx1", "Lhx5", "Sall3", "Rorb", "Pdzd2"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# dILB (Excitatory = Tlx1/3, Lbx1, Lmx1b, Prrxl1, and Brn3a)
FeaturePlot(Data_V1, features = c("Tlx1", "Tlx3", "Lbx1", "Lmx1b", "Prrxl1", "Pou4f1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Tlx1", "Tlx3", "Lbx1", "Lmx1b", "Prrxl1", "Pou4f1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()

# Neural crest / DRG
FeaturePlot(Data_V1, features = c("Sox10", "Tlx2", 'Six1'), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Sox10", 'Tlx2', 'Six1'), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Pan neuron
FeaturePlot(Data_V1, features = c("Snap25", "Enc1", "Uchl1", "Dlg2", "Dlg4", "Rbfox3", "Mapt", "Eno2", "Dcx", 'Tubb3', 'Elavl3', "Syp", "Snhg11"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Snap25", "Enc1", "Uchl1", "Dlg2", "Dlg4", "Rbfox3", "Mapt", "Eno2", "Dcx", 'Tubb3', 'Elavl3', "Syp", "Snhg11"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE)

# Check on known V1 subsets
FeaturePlot(Data_V1, features = c("Foxp2", "St18", "Calb1", "Pou6f2", "Sp8", "Piezo2", "Nr5a2", "Nos1", "Reln"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Add metadata back to seurat object
Data_V1 <- AddMetaData(object = Data_V1, metadata = metadata)

# Plot UMAP with new identities instead of cluster number
DimPlot(Data_V1, group.by = "ID", reduction = "umap", raster = TRUE, label = FALSE, repel = TRUE, pt.size = 8, raster.dpi = c(4096, 4096)) & NoAxes() + theme(aspect.ratio=1)






#### Check positive selection ####



# Load negative selection data
Data_V1 <- readRDS("Data_V1_negSelection3.rds")
DimPlot(Data_V1, group.by = "seurat_clusters", reduction = "umap", raster = TRUE, label = TRUE) & NoAxes() + theme(aspect.ratio=1)

# Load csv file containing barcodes that 
Cells_pos <- read_csv("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/3_PositiveSelectionDelile/3_PostnatalV1Selection/soupx/LabelTransfer/Cells_pos.csv")
Cells_pos$x

highlights = Cells_pos$x
highlights = WhichCells(Data_V1, cells = Cells_pos$x, invert = TRUE)
DimPlot(Data_V1,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = highlights) & NoAxes()



  