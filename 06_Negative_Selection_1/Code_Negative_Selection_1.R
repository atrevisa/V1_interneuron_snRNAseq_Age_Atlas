# 2_Negative_Selection_code on soupX corrected data
# Negative selection without integration
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




#### Analysis NO integration ####




# Load the normalized data the set default assay to RNA
Data_norm <- readRDS("2_Normalization/Data_norm_merged.rds")
DefaultAssay(Data_norm) <- "RNA"

# Select most variable features 
# Note that we are not doing this on each individual dataset but rather all datasets combined
Data_V1 <- FindVariableFeatures(Data_norm, selection.method = "vst", nfeatures = 2000, assay = "RNA")

# Identify the 10 most highly variable genes, check that it worked
top10 <- head(VariableFeatures(Data_V1), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Data_V1)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# PCA
Data_V1 <- RunPCA(Data_V1, features = VariableFeatures(object = Data_V1))

# Visualize PCA results
VizDimLoadings(Data_V1, dims = 1:2, reduction = "pca")
DimPlot(Data_V1, group.by = "age", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "replicate", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "region", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "orig.ident", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimHeatmap(Data_V1, dims = 1, cells = 500, balanced = TRUE, assay = "RNA")
DimHeatmap(Data_V1, dims = 1:9, cells = 500, balanced = TRUE)

# Pick no PC's to use
ElbowPlot(Data_V1)
pcs = find_PC(Data_V1)
pcs

# cluster the cells
Data_V1 <- FindNeighbors(Data_V1, dims = 1:pcs)
Data_V1 <- FindClusters(Data_V1, resolution = 0.1)

# non linear dimension reduction
Data_V1 <- RunUMAP(Data_V1, dims = 1:pcs)

# visualize UMAP
DimPlot(Data_V1, group.by = "seurat_clusters", reduction = "umap", raster = TRUE, label = FALSE, repel = TRUE, raster.dpi = c(4096, 4096), pt.size = 8) & NoAxes() + theme(aspect.ratio=1) + theme(legend.title = element_blank())
ggsave("3_NegSelection_NoInt/Graph_NoIntegration1_UMAP_clusters.pdf", bg = "transparent", device = "pdf")
DimPlot(Data_V1, group.by = "age", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "replicate", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "region", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "orig.ident", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)

# Clean up and save
saveRDS(Data_V1, file = "3_NegSelection_NoInt/Data_noIntegration1.rds")
remove(plot1, pcs, top10, Data_norm)




#### Contaminating cell type ID 1 ####




# Set default assay to RNA
DefaultAssay(Data_V1) <- "RNA"

# Set current ident to clusters
Idents(Data_V1) <- "seurat_clusters"
levels(Data_V1)

# Plot the clusters again
DimPlot(Data_V1, reduction = "umap", raster = TRUE, label = TRUE, repel = TRUE, raster.dpi = c(4096, 4096), pt.size = 3) & NoAxes() + theme(aspect.ratio=1)

# Pull out metadata dataframe
metadata <- Data_V1@meta.data
head(metadata)
# Make a new column called ID and fill all values with "unknown"
metadata$ID <- "Putative-V1"

# Empty GEMS
FeaturePlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
VlnPlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident", add.noise = FALSE) +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()
VlnPlot(Data_V1, features = c("lognFeature_RNA", "lognCount_RNA", "percent_mt", "percent_ribo", "novelty_score"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident") +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend() + theme(aspect.ratio=0.2)
ggsave("3_NegSelection_NoInt/Graph_Vln_QC.pdf", bg = "transparent", device = "pdf")
metadata$ID[metadata$seurat_clusters %in% c(4)] <- "Dead_or_dying"
dplyr::count(metadata, ID, sort = TRUE)

# Glia, general
FeaturePlot(Data_V1, features = c("Sox6"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes()
VlnPlot(Data_V1, features = c("Sox6"), pt.size = 0, slot = "data", stack = FALSE, fill.by = "ident") & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(7, 9, 11)] <- "Glia_general"
dplyr::count(metadata, ID, sort = TRUE)

# Blood-related (endothelial cells, pericytes, etc.)
blood_genes <- list(c("Cldn5", "Rgs5", "Flt1", "Slco1c1", "Fli1", "Sox17", "Fermt3", "Klf1", "Hemgn", "Car2", "Pecam1", "Tek", "Pdgfra", "Flt4"))
FeaturePlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', 'Fli1', 'Sox17', 'Fermt3', 'Klf1', 'Hemgn', 'Car2', "Pecam1", "Tek", "My19", "Pdgfra", "Flt4"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', 'Fli1', "Pecam1", "Tek"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
Data_V1 <- AddModuleScore(object = Data_V1, features = blood_genes,ctrl = 5,name = 'Blood_Features')
FeaturePlot(Data_V1, features = c("Blood_Features1"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes()
VlnPlot(Data_V1, features = c("Blood_Features1"), pt.size = 0, slot = "data", stack = FALSE, fill.by = "ident", flip = TRUE) & NoLegend()
VlnPlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', "Pecam1", "Tek"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(12)] <- "Blood-related"
dplyr::count(metadata, ID, sort = TRUE)

# Microglia
u_genes <- list(c("Aif1", "Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1", "Cd68", "Adgre1", "Mertk", "Fcer1g", "Fcrls", "Hexb", "Sall1"))
FeaturePlot(Data_V1, features = c("Aif1", "Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1", "Cd68", "Adgre1", "Mertk", "Fcer1g", "Fcrls", "Hexb", "Sall1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()
Data_V1 <- AddModuleScore(object = Data_V1, features = u_genes,ctrl = 5,name = 'uglia_Features')
FeaturePlot(Data_V1, features = c("uglia_Features1"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes()
VlnPlot(Data_V1, features = c("uglia_Features1"), pt.size = 0, slot = "data", stack = FALSE, fill.by = "ident", flip = TRUE) & NoLegend()
metadata$ID[metadata$seurat_clusters == 14] <- "Microglia"
dplyr::count(metadata, ID, sort = TRUE)

# Astrocytes
astro_genes <- list(c("Sox6", "Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Glul", "Gja1", "Slc1a3", "Slc4a4", "Sox2", "Slc1a2","S100b", "Ndrg2", "Nkx6-1", "Sox9"))
FeaturePlot(Data_V1, features = c("Sox6", "Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Glul", "Gja1", "Slc1a3", "Slc4a4", "Sox2", "Slc1a2", "S100b", "Ndrg2", "Nkx6-1", "Sox9"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Slc1a3"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()
Data_V1 <- AddModuleScore(object = Data_V1, features = astro_genes,ctrl = 5,name = 'Astro_Features')
FeaturePlot(Data_V1, features = c("Astro_Features1"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes()
VlnPlot(Data_V1, features = c("Astro_Features1"), pt.size = 0, slot = "data", stack = FALSE, fill.by = "ident", flip = TRUE) & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(9, 11)] <- "Astrocytes"
dplyr::count(metadata, ID, seurat_clusters, sort = TRUE)

# Oligodendrocyte lineage cells
# overall
FeaturePlot(Data_V1, features = c("Olig2", "Pdgfra", "Mbp", "Mog", "Mag", "Bcas1", "Mobp", "Plp1", "Mpz", "Pmp22", "Prx", "Cspg4", "Gpr17", "Sox10", "Cnp", "Olig1", "Nkx2-2", "Cd9", "Zfp488", "Zfp536","Nkx6-2", "Cd82", "Mal", "Bmp4","Aspa"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
# opc
FeaturePlot(Data_V1, features = c("Pdgfra", "Cspg4", "Gpr17", "Nkx2-2", "Cd9", "Bmp4"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
# mature
FeaturePlot(Data_V1, features = c("Mbp", "Mog", "Mag", "Cnp", "Mbp", "Mobp", "Plp1", "Zfp536", "Nkx6-2", "Cd82","Mal","Aspa"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
# both
FeaturePlot(Data_V1, features = c("Sox10", "Olig2", "Bcas1", "Olig1", "Zfp488", "St18"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()

VlnPlot(Data_V1, features = c("Mag", "Mog",  "Bcas1","Mobp", "Plp1", "Mbp", "Cspg4", "Pdgfra"), pt.size = 0, slot = "data", fill.by = "ident", stack = TRUE, flip = TRUE) & NoLegend()
VlnPlot(Data_V1, features = c("Mag",  "Mobp", "Mog","Bcas1", "Cspg4", "Pdgfra"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()
oligo_genes <- list(c("Sox10", "Olig2", "Bcas1", "Olig1", "Zfp488", "St18"))
Data_V1 <- AddModuleScore(object = Data_V1, features = oligo_genes,ctrl = 5,name = 'Oligo_Features')
FeaturePlot(Data_V1, features = c("Oligo_Features1"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = "Oligo_Features1", pt.size = 0, slot = "data", fill.by = "ident", stack = FALSE, flip = TRUE) & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(3, 7, 6)] <- "Oligo-lineage"
dplyr::count(metadata, ID, sort = TRUE)

# Excitatory neuron and inhibitory neuron markers
excitatory_genes <- list(c("Slc17a6", "Lmx1b", "Ebf2", "Sox5", "Slc17a7", "Ebf1", "Ebf3", "Cacna2d1"))
FeaturePlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2", "Sox5", "Slc17a7", "Grin1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Cacna2d1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
Data_V1 <- AddModuleScore(object = Data_V1, features = excitatory_genes,ctrl = 5,name = 'Excite_Features')
FeaturePlot(Data_V1, features = c("Excite_Features1"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes()
VlnPlot(Data_V1, features = "Excite_Features1", pt.size = 0, slot = "data", fill.by = "ident", stack = FALSE, flip = TRUE) & NoLegend()
VlnPlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip =TRUE) & NoLegend()
FeaturePlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE) & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(5)] <- "Excitatory-Neurons"
dplyr::count(metadata, ID, sort = TRUE)

# Motor
FeaturePlot(Data_V1, features = c("Chat", "Slc18a3", "Slc5a7", "Isl1"), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Chat", "Slc18a3", "Slc5a7", "Isl1"), pt.size = 0, slot = "counts")

# CSF contacting Neurons
FeaturePlot(Data_V1, features = c("Pkd2l1", "Pkd1l2", "Myo3b"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Pkd2l1", "Pkd1l2", "Myo3b"), pt.size = 0, slot = "counts")

# Check for V1 clusters
FeaturePlot(Data_V1, features = c("Foxp2", "Pou6f2", "Chrna2", "St18"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes() & NoLegend()

# Add metadata back to seurat object
Data_V1 <- AddMetaData(object = Data_V1, metadata = metadata)

# Export UMAP with idents
DimPlot(Data_V1, group.by = "ID", reduction = "umap", raster = TRUE, label = FALSE, repel = TRUE, pt.size = 8, raster.dpi = c(4096, 4096)) & NoAxes() + theme(aspect.ratio=1)
ggsave("3_NegSelection_NoInt/Graph_NoIntegration_UMAP_ID.pdf", bg = "transparent", device = "pdf")

# Put violin plot of module scores together
VlnPlot(Data_V1, features = c("Oligo_Features1",  "Excite_Features1", "Astro_Features1", "Blood_Features1", "uglia_Features1"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) + stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend() + theme(aspect.ratio=0.2)
ggsave("3_NegSelection_NoInt/Graph_NoIntegration_Vln_SCORES.pdf", bg = "transparent", device = "pdf")




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

# clean up
remove(astro_genes, blood_genes, excitatory_genes, metadata, oligo_genes, u_genes, Cells_V1)




#### Analysis NO integration 2 ####




# Set default assay
DefaultAssay(Data_V1) <- "RNA"

# Select most variable features 
# Note that we are not doing this on each individual dataset but rather all datasets combined
Data_V1 <- FindVariableFeatures(Data_V1, selection.method = "vst", nfeatures = 2000, assay = "RNA")

# Identify the 10 most highly variable genes, check that it worked
top10 <- head(VariableFeatures(Data_V1), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Data_V1)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# PCA
Data_V1 <- RunPCA(Data_V1, features = VariableFeatures(object = Data_V1))

# Visualize PCA results
VizDimLoadings(Data_V1, dims = 1:2, reduction = "pca")
DimPlot(Data_V1, group.by = "age", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "replicate", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "region", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimPlot(Data_V1, group.by = "orig.ident", reduction = "pca", raster = TRUE, shuffle = TRUE)
DimHeatmap(Data_V1, dims = 1, cells = 500, balanced = TRUE, assay = "RNA")
DimHeatmap(Data_V1, dims = 1:9, cells = 500, balanced = TRUE)

# Pick no PC's to use, quick
ElbowPlot(Data_V1) # 10
pcs = find_PC(Data_V1)
pcs

# cluster the cells
Data_V1 <- FindNeighbors(Data_V1, dims = 1:pcs)
Data_V1 <- FindClusters(Data_V1, resolution = 0.5)

# non linear dimension reduction
Data_V1 <- RunUMAP(Data_V1, dims = 1:pcs)

# UMAP
DimPlot(Data_V1, group.by = "seurat_clusters", reduction = "umap", raster = TRUE, label = FALSE, repel = TRUE, raster.dpi = c(4096, 4096), pt.size = 8) & NoAxes() + theme(aspect.ratio=1) + theme(legend.title = element_blank())
DimPlot(Data_V1, group.by = "age", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "replicate", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "region", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, group.by = "orig.ident", reduction = "umap", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)

# visualize features on the UMAP
colnames(Data_V1@meta.data)
metrics <-  c("nFeature_RNA", "nCount_RNA", "percent_mt")
FeaturePlot(Data_V1, reduction = "umap", features = metrics, pt.size = 1, label = TRUE, raster = TRUE, repel = TRUE) & NoAxes() + theme(aspect.ratio=1)

# Clean up and save
remove(plot1, metrics, top10, pcs)
saveRDS(Data_V1, file = "3_NegSelection_NoInt/Data_noIntegration2.rds")





#### Contaminating cell type ID 2 ####




# Load data, if necessary
Data_V1 <- readRDS("3_NegSelection_NoInt/Data_noIntegration2.rds")

# Set default assay to RNA
DefaultAssay(Data_V1) <- "RNA"

# Set current ident to clusters
Idents(Data_V1) <- "seurat_clusters"
levels(Data_V1)

# plot the clusters again
DimPlot(Data_V1, reduction = "umap", raster = TRUE, label = TRUE, repel = TRUE, raster.dpi = c(4096, 4096), pt.size = 3) & NoAxes() + theme(aspect.ratio=1)

# pull out metadata dataframe
metadata <- Data_V1@meta.data
head(metadata)
# make a new column called ID and fill all values with "unknown"
metadata$ID <- "Putative-V1"

# empty GEMS
FeaturePlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "novelty_score"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
FeaturePlot(Data_V1, features = c("nFeature_RNA"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
VlnPlot(Data_V1, features = c("lognFeature_RNA", "lognCount_RNA", "percent_mt"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident") +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()
VlnPlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident", add.noise = FALSE) +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()
metadata$ID[metadata$seurat_clusters %in% c(2, 8)] <- "Null"
dplyr::count(metadata, ID, sort = TRUE)

# Glia, general
FeaturePlot(Data_V1, features = c("Sox6"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes()
VlnPlot(Data_V1, features = c("Sox6"), pt.size = 0, slot = "data", stack = FALSE, fill.by = "ident") & NoLegend()

# Blood-related (endothelial cells, pericytes, etc.)
FeaturePlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', 'Fli1', "Pecam1", "Tek"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', "Pecam1", "Tek"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Microglia
FeaturePlot(Data_V1, features = c("Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Astrocytes
FeaturePlot(Data_V1, features = c("Sox6", "Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Glul", "Gja1", "Slc1a3", "Slc4a4", "Sox2", "Slc1a2", "S100b", "Ndrg2", "Nkx6-1", "Sox9"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Sox6", "Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Slc1a3"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Oligodendrocyte lineage cells
# overall
FeaturePlot(Data_V1, features = c("Olig2", "Pdgfra", "Mbp", "Mog", "Mag", "Bcas1", "Mobp", "Plp1", "Mpz", "Pmp22", "Prx", "Cspg4", "Gpr17", "Sox10", "Cnp", "Olig1", "Nkx2-2", "Cd9", "Zfp488", "Zfp536","Nkx6-2", "Cd82", "Mal", "Bmp4","Aspa"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Mbp", "Mog", "Mag", "Bcas1", "Mobp", "Plp1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Mbp", "Mog", "Mag", "Bcas1", "Mobp", "Plp1"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip =TRUE, add.noise = TRUE) & NoLegend()
# opc
FeaturePlot(Data_V1, features = c("Pdgfra", "Cspg4", "Gpr17", "Nkx2-2", "Cd9", "Bmp4"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
# mature
FeaturePlot(Data_V1, features = c("Mbp", "Mog", "Mag", "Cnp", "Mbp", "Mobp", "Plp1", "Zfp536", "Nkx6-2", "Cd82","Mal","Aspa"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
# all types
FeaturePlot(Data_V1, features = c("Sox10", "Olig2", "Bcas1", "Olig1", "Zfp488", "St18"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()

# Excitatory / Onhibitory markers
FeaturePlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2", "Sox5", "Slc17a7", "Grin1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip =TRUE) & NoLegend()
FeaturePlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1", "percent_mt"), pt.size = 0, slot = "data", stack = TRUE)
metadata$ID[metadata$seurat_clusters %in% c(7, 13)] <- "Excitatory-Neurons"
dplyr::count(metadata, ID, sort = TRUE)

# Meninges
FeaturePlot(Data_V1, features = c("Dcn",  "Col3a1", "Igf2"), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()

# CSF contacting neurons
FeaturePlot(Data_V1, features = c("Pkd2l1", "Pkd1l2", "Myo3b"), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()

# Mesoderm/ myoblast
FeaturePlot(Data_V1, features = c('Foxc1', 'Foxc2', 'Twist1', 'Twist2', 'Meox1', 'Meox2', 'Myog'), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()

# Motor
FeaturePlot(Data_V1, features = c("Chat", "Slc18a3", "Slc5a7", "Isl1", "Ret", "Slit3", "Prph", "Lhx3", "Isl2", "Mnx1", "Slc0a4", "Slc8a3"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Chat"), pt.size = 0, slot = "data", stack = FALSE)

# Ependymal cells
FeaturePlot(Data_V1, features = c("Dnah12", "Spef2", "Ccdc114", "Ddo", "Cfap65", "Ak9", "Fam216b", "Zfp474", "Wdr63", "Ccdc180"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI1
FeaturePlot(Data_V1, features = c("Pou4f1", "Lhx2", "Lhx9", "Barhl1", "Barhl2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI2
FeaturePlot(Data_V1, features = c("Pou4f1", "Foxd3", "Lhx1", "Lhx5"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI3
FeaturePlot(Data_V1, features = c("Pou4f1", "Isl1", "Tlx3", "Prrxl1", "Otp"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI4
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Pax8", "Lbx1", "Pax2", "Gbx1", "Bhlhe22"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI5
FeaturePlot(Data_V1, features = c("Pax2", "Pax7", "Gsx2", "Ascl1", "Gsx1", "Lmx1b", "Pou4f1", "Tlx3", "Prrxl1", "Lbx1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI6
FeaturePlot(Data_V1, features = c("Lhx1", "Lbx1", "Pax2", "Bhlhe22", "Dmrt3", "Wt1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# V0
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Evx1", "Evx2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# V1
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Foxd3", "En1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# V2a
FeaturePlot(Data_V1, features = c("Bhlhe22", "Lhx3", "Vsx2", "Sox14", "Sox21"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# V2b
FeaturePlot(Data_V1, features = c("Bhlhe22", "Tal1", "Gata2", "Gata3"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# V3
FeaturePlot(Data_V1, features = c("Nkx2-2", "Sim1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# RP
FeaturePlot(Data_V1, features = c("Lmx1a", "Msx1", "Msx2", "Pax3", "Wnt1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# FP
FeaturePlot(Data_V1, features = c("Nkx6-2", "Foxa2", "Ferd3l", "Arx", "Shh", "Lmx1b"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Non MN cholinergic
FeaturePlot(Data_V1, features = c("Uncx", "Irx2", "Slc5a7", "Lpcat2","Nwd2", "2900055J20Rik", "Pid1", "Sox6"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Skin
FeaturePlot(Data_V1, features = c('Krt8'), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Dorsal
FeaturePlot(Data_V1, features = c("Gbx1", "Sall3", "Pax2", "Celf2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Other
FeaturePlot(Data_V1, features = c("Rorb", "Gbx1", "Pdzd2", "Arpp21", "Sorcs3"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Other
FeaturePlot(Data_V1, features = c("Gria2", "Celf2", "Tll2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Other
FeaturePlot(Data_V1, features = c("Dlc1", "Cald1", "Rapgef5"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Neural crest / DRG
FeaturePlot(Data_V1, features = c("Sox10", 'Tlx2', 'Six1'), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Progenitors
FeaturePlot(Data_V1, features = c("Sox9", "Fabp2", "Sox2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Pan neuron
FeaturePlot(Data_V1, features = c("Snap25", "Enc1", "Uchl1", "Dlg2", "Dlg4", "Rbfox3", "Mapt", "Eno2", "Dcx", 'Tubb3', 'Elavl3', "Syp", "Snhg11"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Check on known V1 subsets
FeaturePlot(Data_V1, features = c("Foxp2", "Calb1", "Pou6f2", "Sp8", "Piezo2", "Nr5a2", "Nos1", "Reln"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Add metadata back to seurat object
Data_V1 <- AddMetaData(object = Data_V1, metadata = metadata)

# Plot UMAP with new identities instead of cluster number
DimPlot(Data_V1, group.by = "ID", reduction = "umap", raster = TRUE, label = FALSE, repel = TRUE, pt.size = 8, raster.dpi = c(4096, 4096)) & NoAxes() + theme(aspect.ratio=1)




#### Selection 2 ####




# Set ident to ID
Idents(Data_V1) <- "ID"
levels(Data_V1)

# Select cell IDs that are putative V1's
Data_V1
Cells_V1 <- WhichCells(Data_V1, idents = "Putative-V1")
length(Cells_V1)
Data_V1 <- subset(Data_V1, cells = Cells_V1)
Data_V1




#### Analysis NO integration 3 ####




# Select most variable features 
# Note that we are not doing this on each individual dataset but rather all datasets combined
Data_V1 <- FindVariableFeatures(Data_V1, selection.method = "vst", nfeatures = 2000, assay = "RNA")

# Identify the 10 most highly variable genes, check that it worked
top10 <- head(VariableFeatures(Data_V1), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Data_V1)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# PCA
Data_V1 <- RunPCA(Data_V1, features = VariableFeatures(object = Data_V1))

# Visualize PCA results
VizDimLoadings(Data_V1, dims = 1:2, reduction = "pca")
DimPlot(Data_V1, reduction = "pca", raster = TRUE, group.by = "age", shuffle = TRUE)
DimPlot(Data_V1, reduction = "pca", raster = TRUE, group.by = "replicate", shuffle = TRUE)
DimPlot(Data_V1, reduction = "pca", raster = TRUE, group.by = "region", shuffle = TRUE)
DimPlot(Data_V1, reduction = "pca", raster = TRUE, group.by = "orig.ident", shuffle = TRUE)
DimHeatmap(Data_V1, dims = 1, cells = 500, balanced = TRUE, assay = "RNA")
DimHeatmap(Data_V1, dims = 1:9, cells = 500, balanced = TRUE)
# age seems to have the biggest influence on PCA

# Pick no PC's to use, quick
ElbowPlot(Data_V1) # 10
pcs = find_PC(Data_V1)
pcs

# cluster the cells
Data_V1 <- FindNeighbors(Data_V1, dims = 1:pcs)
Data_V1 <- FindClusters(Data_V1, resolution = 0.8)

# non linear dimension reduction
Data_V1 <- RunUMAP(Data_V1, dims = 1:pcs)

# UMAP plots
DimPlot(Data_V1, reduction = "umap", group.by = "seurat_clusters", raster = TRUE, label = FALSE, repel = TRUE, raster.dpi = c(4096, 4096), pt.size = 8) & NoAxes() + theme(aspect.ratio=1) + theme(legend.title = element_blank())
DimPlot(Data_V1, reduction = "umap", group.by = "age", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, reduction = "umap", group.by = "replicate", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, reduction = "umap", group.by = "region", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)
DimPlot(Data_V1, reduction = "umap", group.by = "orig.ident", raster = TRUE) & NoAxes() + theme(aspect.ratio=1)

# visualize QC features on the UMAP
colnames(Data_V1@meta.data)
metrics <-  c("nFeature_RNA", "nCount_RNA", "percent_mt")
FeaturePlot(Data_V1, reduction = "umap", features = metrics, pt.size = 1, label = TRUE, raster = TRUE, repel = TRUE) & NoAxes() + theme(aspect.ratio=1)

# Clean up and save
saveRDS(Data_V1, file = "3_NegSelection_NoInt/Data_noIntegration3.rds")
remove(plot1, metrics, top10, metadata, Cells_V1, pcs)




#### Contaminating cell type ID 3 ####





# Set default assay to RNA
DefaultAssay(Data_V1) <- "RNA"

# Set current ident to clusters
Idents(Data_V1) <- "seurat_clusters"
levels(Data_V1)

# plot the clusters again
DimPlot(Data_V1, reduction = "umap", raster = TRUE, label = TRUE, repel = TRUE, raster.dpi = c(4096, 4096), pt.size = 3) & NoAxes() + theme(aspect.ratio=1)

# pull out metadata dataframe
metadata <- Data_V1@meta.data
head(metadata)
# make a new column called ID and fill all values with "unknown"
# only do this once
metadata$ID <- "Putative-V1"

# empty GEMS
FeaturePlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "novelty_score"), label = FALSE, pt.size = 1, raster = TRUE, slot = "counts") & NoAxes()
VlnPlot(Data_V1, features = c("lognFeature_RNA", "lognCount_RNA", "percent_mt"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident") +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()
VlnPlot(Data_V1, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), pt.size = 0, slot = "data", stack = TRUE, flip = TRUE, fill = "ident", add.noise = FALSE) +
  stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95) & NoLegend()

# Glia, general
FeaturePlot(Data_V1, features = c("Sox6"), label = TRUE, pt.size = 1, raster = FALSE, slot = "data") & NoAxes()
VlnPlot(Data_V1, features = c("Sox6"), pt.size = 0, slot = "data", stack = FALSE, fill.by = "ident") & NoLegend()

# Blood-related (endothelial cells, pericytes, etc.)
FeaturePlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', 'Fli1', "Pecam1", "Tek"), label = FALSE, pt.size = 1, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Cldn5", "Rgs5", "Flt1", 'Slco1c1', "Pecam1", "Tek"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Microglia
FeaturePlot(Data_V1, features = c("Tmem119", "Trem2", "Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Inpp5d", "Ctss", "Itgam", "Ptprc", "Cx3cr1"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Astrocytes
FeaturePlot(Data_V1, features = c("Sox6", "Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Glul", "Gja1", "Slc1a3", "Slc4a4", "Sox2", "Slc1a2", "S100b", "Ndrg2", "Nkx6-1", "Sox9"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Sox6", "Aldh1l1", "Fgfr3", "Aqp4", "Gfap", "Slc1a3"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip = TRUE) & NoLegend()

# Oligodendrocyte lineage
# overall
FeaturePlot(Data_V1, features = c("Olig2", "Pdgfra", "Mbp", "Mog", "Mag", "Bcas1", "Mobp", "Plp1", "Mpz", "Pmp22", "Prx", "Cspg4", "Gpr17", "Sox10", "Cnp", "Olig1", "Nkx2-2", "Cd9", "Zfp488", "Zfp536","Nkx6-2", "Cd82", "Mal", "Bmp4","Aspa"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
FeaturePlot(Data_V1, features = c("Mbp", "Mog", "Mag", "Bcas1", "Mobp", "Plp1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Mbp", "Mog", "Mag", "Bcas1", "Mobp", "Plp1"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip =TRUE, add.noise = TRUE) & NoLegend()
# opc
FeaturePlot(Data_V1, features = c("Pdgfra", "Cspg4", "Gpr17", "Nkx2-2", "Cd9", "Bmp4"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
# mature
FeaturePlot(Data_V1, features = c("Mbp", "Mog", "Mag", "Cnp", "Mbp", "Mobp", "Plp1", "Zfp536", "Nkx6-2", "Cd82","Mal","Aspa"), label = FALSE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()
# all types
FeaturePlot(Data_V1, features = c("Sox10", "Olig2", "Bcas1", "Olig1", "Zfp488", "St18"), label = TRUE, raster = TRUE, slot = "data", repel = TRUE) & NoAxes() & NoLegend()

# Excitatory / Inhibitory neuron markers
FeaturePlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2", "Sox5", "Slc17a7", "Grin1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Slc17a6", "Lmx1b", "Ebf2"), pt.size = 0, slot = "data", stack = TRUE, fill.by = "ident", flip =TRUE) & NoLegend()
FeaturePlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Gad1", "Gad2", "Slc6a5", "Slc32a1"), pt.size = 0, slot = "data", stack = TRUE)

# Meninges
FeaturePlot(Data_V1, features = c("Dcn",  "Col3a1", "Igf2"), label = TRUE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()

# CSF contacting neurons
FeaturePlot(Data_V1, features = c("Pkd2l1", "Pkd1l2", "Myo3b"), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()

# Mesoderm/ myoblast 
FeaturePlot(Data_V1, features = c('Foxc1', 'Foxc2', 'Twist1', 'Twist2', 'Meox1', 'Meox2', 'Myog'), label = FALSE, raster = FALSE, slot = "data") & NoAxes() & NoLegend()

# Motor
FeaturePlot(Data_V1, features = c("Chat", "Slc18a3", "Slc5a7", "Isl1", "Ret", "Slit3", "Prph", "Lhx3", "Isl2", "Mnx1", "Slc0a4", "Slc8a3"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()
VlnPlot(Data_V1, features = c("Chat"), pt.size = 0, slot = "data", stack = FALSE)

# Ependymal cells
FeaturePlot(Data_V1, features = c("Dnah12", "Spef2", "Ccdc114", "Ddo", "Cfap65", "Ak9", "Fam216b", "Zfp474", "Wdr63", "Ccdc180"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI1
FeaturePlot(Data_V1, features = c("Pou4f1", "Lhx2", "Lhx9", "Barhl1", "Barhl2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI2
FeaturePlot(Data_V1, features = c("Pou4f1", "Foxd3", "Lhx1", "Lhx5"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI3
FeaturePlot(Data_V1, features = c("Pou4f1", "Isl1", "Tlx3", "Prrxl1", "Otp"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI4
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Pax8", "Lbx1", "Pax2", "Gbx1", "Bhlhe22"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI5
FeaturePlot(Data_V1, features = c("Pax2", "Pax7", "Gsx2", "Ascl1", "Gsx1", "Lmx1b", "Pou4f1", "Tlx3", "Prrxl1", "Lbx1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# dI6
FeaturePlot(Data_V1, features = c("Lhx1", "Lbx1", "Pax2", "Bhlhe22", "Dmrt3", "Wt1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# V0
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Evx1", "Evx2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# V1
FeaturePlot(Data_V1, features = c("Lhx1", "Lhx5", "Foxd3", "En1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# V2a
FeaturePlot(Data_V1, features = c("Bhlhe22", "Lhx3", "Vsx2", "Sox14", "Sox21"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# V2b 
FeaturePlot(Data_V1, features = c("Bhlhe22", "Tal1", "Gata2", "Gata3"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# V3
FeaturePlot(Data_V1, features = c("Nkx2-2", "Sim1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# RP
FeaturePlot(Data_V1, features = c("Lmx1a", "Msx1", "Msx2", "Pax3", "Wnt1"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# FP
FeaturePlot(Data_V1, features = c("Nkx6-2", "Foxa2", "Ferd3l", "Arx", "Shh", "Lmx1b"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Non MN cholinergic
FeaturePlot(Data_V1, features = c("Uncx", "Irx2", "Slc5a7", "Lpcat2","Nwd2", "2900055J20Rik", "Pid1", "Sox6"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Skin
FeaturePlot(Data_V1, features = c('Krt8'), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Dorsal
FeaturePlot(Data_V1, features = c("Gbx1", "Sall3", "Pax2", "Celf2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Other
FeaturePlot(Data_V1, features = c("Rorb", "Gbx1", "Pdzd2", "Arpp21", "Sorcs3"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Other
FeaturePlot(Data_V1, features = c("Gria2", "Celf2", "Tll2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Other
FeaturePlot(Data_V1, features = c("Dlc1", "Cald1", "Rapgef5"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Neural crest / DRG
FeaturePlot(Data_V1, features = c("Sox10", 'Tlx2', 'Six1'), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Progenitors
FeaturePlot(Data_V1, features = c("Sox9", "Fabp2", "Sox2"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Pan neuron
FeaturePlot(Data_V1, features = c("Snap25", "Enc1", "Uchl1", "Dlg2", "Dlg4", "Rbfox3", "Mapt", "Eno2", "Dcx", 'Tubb3', 'Elavl3', "Syp", "Snhg11"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()

# Check on known V1 subsets
FeaturePlot(Data_V1, features = c("Foxp2", "Calb1", "Chrna2", "Pou6f2", "Sp8", "Piezo2", "Nr5a2", "Nos1", "Reln"), label = FALSE, raster = TRUE, slot = "data") & NoAxes() & NoLegend()




#### Final Negative Selection ####




# Export the final list of putative V1 barcodes
Cells_V1 <- WhichCells(Data_V1)
length(Cells_V1)
write.csv(Cells_V1, file ="3_NegSelection_NoInt/Cells_PutativeV1_NegSelect.csv", row.names=FALSE)
remove(metadata, Cells_V1)
