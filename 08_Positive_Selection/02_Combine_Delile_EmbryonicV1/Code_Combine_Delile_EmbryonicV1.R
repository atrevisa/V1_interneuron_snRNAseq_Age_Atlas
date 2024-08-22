# Delile data and embryonic Bikoff data combined
# AT Updated 5/8/24

#### Setup ####





# Set up default directory
setwd("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/3_PositiveSelectionDelile/2_Delile_Bikoff_Embryonic_together/")

# Install necessary packages if not already installed
#source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/0_IstallPackages.R")

# load packages if necessary
source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/1_LoadLibraries.R")

# Import custom functions
sourceDirectory("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/", modifiedOnly=FALSE)




#### Import data ####




# Import the Delile data, it has already been analyzed and integrated
# Note that it has been log normalized
Data_reference <- readRDS("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/3_PositiveSelectionDelile/0_Delile_alone/Data_IntegrationNoSCT_Labeled.rds")

# Peek at the Delile metadata
# I will be using the "Specific_Identities" as my reference cell types
colnames(Data_reference@meta.data)

# QUICKLY LOOK AT THE REFERENCE DATA
# Specify the order of the levels
DefaultAssay(Data_reference) <- 'integrated'
# Set default ident
Idents(Data_reference) <- "Specific_Identities"
levels(Data_reference)
# Set order that you want the cell types to appear
orders = rev(c( "Null", "Progenitors", "Mesoderm", "Blood", "NC","DRG", "dI1", "dI2","dI3", "dI4", "dI5", "dI6","V0", "V2", "MN", "V3", "V1"))
# Specify the colors you want each cluster to be 
color_list = c( "Null" = "#cdcd69", "Progenitors" = "#ff37cd", 
                "Mesoderm" = "#379bcd", "Blood" = "#ff6969", 
                "NC" = "#cd69cd","DRG" = "#cd6937", 
                "dI1" = "#ff3769", "dI2" = "#9b9bff","dI3" = "#ffcd69", "dI4" = "#ff6937", 
                "dI5" = "#9b69ff", "dI6" = "#69ffff","V0" = "#37ff69",
                "V1" = "#ff9b37", "V2" = "#37ffcd", 
                "MN" = "#cd3769", "V3" = "#37379b")
# Visualize the query cells with their cell type identities
DimPlot(Data_reference, reduction = "umap", group.by = "Specific_Identities", label = TRUE, repel = FALSE, pt.size = 0.5, label.size = 4, 
        cols =  color_list, order = orders) + ggtitle("Reference annotations") + NoAxes() + theme(aspect.ratio=1)
# Visualize QC stats on the query umap
FeaturePlot(Data_reference, features = c("nFeature_RNA", "percent_mt")) & NoAxes()




# Import the embryonic Bikoff data, it has already been analyzed
# Note that it has been log normalized
Data_query <- readRDS("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/3_PositiveSelectionDelile/1_Bikoff_embryonic_alone/Data_NoIntegration.rds")

# Peek at the Bikoff metadata
colnames(Data_query@meta.data)

# QUICKLY LOOK AT THE QUIERY DATA
# set default assay (there is only one)
DefaultAssay(Data_query) <- "RNA"
# set identity
Idents(Data_query) <- "source"
levels(Data_query)
# Visualize the Bikoff data before cell type identification (black)
DimPlot(Data_query, group.by = "dummy_variable", label = FALSE, repel = TRUE, pt.size = 0.5, cols = "black") + 
  ggtitle("Query dataset") + NoAxes() + theme(aspect.ratio=1)
# Visualize QC metrics on Bikoff data
FeaturePlot(Data_query, features = c("nFeature_RNA", "percent_mt")) & NoAxes()
# Visualize En1 on Bikoff data
FeaturePlot(Data_query, features = c("En1"),order = TRUE, min.cutoff = 0, max.cutoff = 2, pt.size = 1) & NoAxes()




#### Label transfer without integration ####




# Find the query anchors. Be sure to set normalization method to SCT, since that's how both datasets were normalized
Data_query.anchors <- FindTransferAnchors(reference = Data_reference, 
                                          reference.assay = "integrated", 
                                          query = Data_query, 
                                          query.assay = "RNA", 
                                          k.anchor = 16, 
                                          reference.reduction = "pca")

# Predict cell types using the anchors from above
Data_predictions <- TransferData(anchorset = Data_query.anchors, refdata = Data_reference$Specific_Identities)

# Add the predicted cell types back to the Bikoff metadata
Data_query <- AddMetaData(Data_query, metadata = Data_predictions)

# examine the results
table(Data_query$predicted.id)

# View predict cell type
DimPlot(Data_query, group.by = "predicted.id", label = TRUE, repel = TRUE, cols =  color_list, order = orders) & NoAxes()

# save the reference data with the predicted id's back in the original folder and this folder
saveRDS(Data_query, "/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/3_PositiveSelectionDelile/1_Bikoff_embryonic_alone/Data_LabelTransfer.rds")
saveRDS(Data_query, "/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/3_PositiveSelectionDelile/2_Delile_Bikoff_Embryonic_together/Data_Bikoff_LabelTransfer.rds")




# RunUMAP
Data_reference <- RunUMAP(Data_reference, dims = 1:49, n.neighbors = 30L, reduction = "pca", return.model = TRUE)

# Project the unknown dataset onto the reference UMAP
Data_query <- MapQuery(anchorset = Data_query.anchors, reference = Data_reference, query = Data_query, refdata = list(celltype = "Specific_Identities"), reference.reduction = "pca", reduction.model = "umap")

# Visualize the query cells with their cell type identities
DimPlot(Data_reference, reduction = "umap", group.by = "Specific_Identities", label = TRUE, repel = FALSE, pt.size = 0.5, label.size = 4, cols =  color_list, order = orders) + ggtitle("Reference annotations") + NoAxes() + theme(aspect.ratio=1)

# Visualize the Bikoff data before cell type identification (black)
DimPlot(Data_query, reduction = "ref.umap", group.by = "source", label = FALSE, repel = TRUE, pt.size = 0.5, cols = "black") + ggtitle("Query dataset") + NoAxes() + theme(aspect.ratio=1)

# Visualize Bikoff data after cell type identification
DimPlot(Data_query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, repel = TRUE, pt.size = 0.5, label.size = 4,cols = color_list, order = orders, shuffle = FALSE) + ggtitle("Query transferred labels") + NoAxes() + theme(aspect.ratio=1)

# Now graph the confidence of the label transfer score
FeaturePlot(Data_query, reduction = "ref.umap", features = "prediction.score.V1", label = FALSE, repel = TRUE, pt.size = 0.5) + ggtitle("V1 prediction score") + NoAxes() + theme(aspect.ratio=1)

# Now merge the label transferred datasets back together - note these are NOT integrated
Data_label_transfer <- merge(Data_reference, y = Data_query)
Data_label_transfer 
head(Data_label_transfer@meta.data)
tail(Data_label_transfer@meta.data)

# Cells per sample
dplyr::count(Data_label_transfer@meta.data, orig.ident, sort = TRUE)
dplyr::count(Data_label_transfer@meta.data, predicted.celltype, sort = TRUE)

# The following line of code is needed to merge the umap plots together
Data_label_transfer[["umap"]] <- merge(Data_reference[["umap"]], Data_query[["ref.umap"]])

# Graph both sets of data overlaid with labels
DimPlot(Data_label_transfer, group.by = "Specific_Identities", na.value = "black", shuffle = TRUE, label = TRUE, repel = TRUE, order = orders, cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1)

# Graph both sets of data overlaid without labels
DimPlot(Data_label_transfer, group.by = "Specific_Identities", na.value = "black", shuffle = TRUE, label = FALSE, repel = TRUE, order = NULL, cols = DiscretePalette(17, palette = "parade", shuffle = TRUE), pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1)

# Graph En1 expression on the UMAP
DefaultAssay(Data_label_transfer) <- "RNA"
FeaturePlot(Data_label_transfer, features = c("En1"), label = FALSE, min.cutoff = NA, max.cutoff = NA, slot = "data", order = TRUE) & NoAxes() + theme(aspect.ratio=1)




#### Integration ####




# First merge all of the data together
Data_integrated <- merge(Data_reference, y = Data_query)
Data_reference 
Data_query 
Data_integrated
# even though I called it integrated this data is merged right now

# Set default assay - if you don't do this now downstream stuff will not work
DefaultAssay(object = Data_integrated) <- "RNA"

# Now split the merged object into a list
# split the dataset into a list of seurat objects 
Data_list <- SplitObject(Data_integrated, split.by = "orig.ident")
Data_list

# Sormalize and identify variable features for each dataset independently
# Scale the data individually
# Also note that both datasets were separately normalized but I'm going to do it again to ensure it is all done the same way
# Find individual variable features for each dataset independently
Data_list <- lapply(X = Data_list, FUN = function(x) {
  x <- NormalizeData(x, assay = "RNA")
  all.genes <- rownames(x)
  x <- ScaleData(x, features = all.genes, assay = "RNA")
  # be sure to set assay to RNA
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, assay ="RNA")
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Data_list)

# FindIntegrationAnchors takes a while to run
anchors <- FindIntegrationAnchors(object.list = Data_list, anchor.features = features)

# This command creates an 'integrated' data assay. Also takes a while to run
Data_integrated <- IntegrateData(anchorset = anchors)

# Set default assay
DefaultAssay(Data_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
Data_integrated <- ScaleData(Data_integrated, verbose = TRUE)
Data_integrated <- RunPCA(Data_integrated, npcs = 100, verbose = TRUE)
# Choose PCs to use for downstream analysis
Data_integrated <- RunUMAP(Data_integrated, reduction = "pca", dims = 1:78, n.neighbors = 30L) # Need at least 75 dims, 30 neighbors

# Modify the metadata 
metadata <- Data_integrated@meta.data
colnames(metadata)
# Add source data, bikoff or delile in a new column called source
for (i in 1:nrow(metadata)) {
  #print(i)
  #print(metadata$orig.ident[i])
  #print(startsWith(metadata$orig.ident[i], "D"))
  if(startsWith(as.character(metadata$orig.ident[i]), "D") == TRUE) {
    metadata$source[i] <- "Delile"
  } else {
    metadata$source[i] <- "Bikoff"
  }
}
# Create a new column in which the final identities for all cells reside
# start by making a new column called Final_identities where you copy over the predicted.celltype
metadata$Final_identities <- metadata$predicted.id
# now replace the NA's with the Specific Identities of the Delile data
metadata %>% mutate(Final_identities = coalesce(Final_identities,Specific_Identities))
# Add the updated metadata back to the seurat object
Data_integrated<- AddMetaData(object = Data_integrated, metadata = metadata)
head(Data_integrated@meta.data, 50)
dplyr::count(Data_integrated@meta.data, Final_identities, sort = TRUE)




#### Generate plots using the integrated data ####




# Specify the order of the levels, slightly different than before
# Specify the colors
orders = rev(c("Query", "Null", "Progenitors", "Mesoderm", "Blood", "NC","DRG", "dI1", "dI2","dI3", "dI4", "dI5", "dI6","V0","V1", "V2", "MN", "V3"))
color_list = c("Query" = "black", "Null" = "#cdcd69", "Progenitors" = "#ff37cd", 
                "Mesoderm" = "#379bcd", "Blood" = "#ff6969", 
                "NC" = "#cd69cd","DRG" = "#cd6937", 
                "dI1" = "#ff3769", "dI2" = "#9b9bff","dI3" = "#ffcd69", "dI4" = "#ff6937", 
                "dI5" = "#9b69ff", "dI6" = "#69ffff","V0" = "#37ff69",
                "V1" = "#ff9b37", "V2" = "#37ffcd", 
                "MN" = "#379b37", "V3" = "#37379b")

# Visualization of both datasets together
DimPlot(Data_integrated, group.by = "Specific_Identities", na.value = "black", shuffle = TRUE, label = TRUE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Reference and Query")
ggsave("Graph_UMAP_BothDatasets_Labels.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

DimPlot(Data_integrated, group.by = "Specific_Identities", na.value = "black", shuffle = TRUE, label = FALSE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Reference and Query")
ggsave("Graph_UMAP_BothDatasets_Legend.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

# Visualize En1 expression on both datasets
DefaultAssay(Data_integrated) <- "RNA"
FeaturePlot(Data_integrated, features = c("En1"), label = FALSE, min.cutoff = NA, max.cutoff = NA, slot = "data", order = TRUE) & 
  NoAxes() + theme(aspect.ratio=1)
ggsave("Graph_UMAP_BothDatasets_En1.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

# Separate the datasets to visualize separately
Idents(Data_integrated) <- "source"
levels(Data_integrated)
# subset
Data_Delile_integrated <- subset(Data_integrated, idents = "Delile")
Data_Delile_integrated
Data_Bikoff_integrated <- subset(Data_integrated, idents = "Bikoff")
Data_Bikoff_integrated

# Visualize Delile data alone
DimPlot(Data_Delile_integrated, group.by = "Specific_Identities", na.value = "black", shuffle = FALSE, label = TRUE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Reference data")
ggsave("Graph_UMAP_DelileReference_Labels.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

DimPlot(Data_Delile_integrated, group.by = "Specific_Identities", na.value = "black", shuffle = FALSE, label = FALSE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Reference data")
ggsave("Graph_UMAP_DelileReference_Legend.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

# Visualize Bikoff data alone
DimPlot(Data_Bikoff_integrated, group.by = "Specific_Identities", na.value = "black", shuffle = FALSE, label = FALSE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Query data")
ggsave("Graph_UMAP_Bikoff_Query.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

# Visualize assigned Bikoff identities with labels
DimPlot(Data_Bikoff_integrated, group.by = "predicted.id", na.value = "black", shuffle = FALSE, label = TRUE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Query label transfer")
ggsave("Graph_UMAP_Bikoff_Labels.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

# Visualize assigned Bikoff identities without labels
DimPlot(Data_Bikoff_integrated, group.by = "predicted.id", na.value = "black", shuffle = FALSE, label = FALSE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Query label transfer")
ggsave("Graph_UMAP_Bikoff_Legend.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

# Now graph the confidence of the label transfer score
# FeaturePlot(Data_Bikoff_integrated, features = "prediction.score.max", label = FALSE, repel = TRUE, pt.size = 0.5) + ggtitle("Maximum confidence Score") + NoAxes() + theme(aspect.ratio=1)
FeaturePlot(Data_Bikoff_integrated, features = "prediction.score.V1", label = FALSE, repel = TRUE, pt.size = 0.5) + ggtitle("V1 prediction score") + NoAxes() + theme(aspect.ratio=1)
ggsave("Graph_UMAP_Bikoff_Scores.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

# Visualize all of the data together that is labeled
DimPlot(Data_integrated, group.by = "Final_identities", na.value = "black", shuffle = TRUE, label = TRUE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Reference and Query")
ggsave("Graph_UMAP_BothDatasets_Labels_PostTransfer_Labels.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

# Visualize all of the data together that is labeled
DimPlot(Data_integrated, group.by = "Final_identities", na.value = "black", shuffle = TRUE, label = FALSE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Reference and Query")
ggsave("Graph_UMAP_BothDatasets_Labels_PostTransfer_Legend.pdf", bg = "transparent", device = "pdf", width = 6, height = 6)

# save and clean up
saveRDS(Data_integrated, file = "Data_integrated.rds")
remove(anchors, Data_list, metadata, i, features, metadata)
# read data back in if necessary
Data_integrated <- readRDS("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/3_PositiveSelectionDelile/2_Delile_Bikoff_Embryonic_together/Data_integrated.rds")




#### Label transer after integration #####




# Separate the datasets to visualize separately
Idents(Data_integrated) <- "source"
levels(Data_integrated)
# subset
Data_Delile_integrated <- subset(Data_integrated, idents = "Delile")
Data_Delile_integrated
Data_Bikoff_integrated <- subset(Data_integrated, idents = "Bikoff")
Data_Bikoff_integrated

# Use the integrated data for reference and query
Data_query.anchors <- FindTransferAnchors(reference = Data_Delile_integrated, 
                                          reference.assay = "integrated", 
                                          query = Data_Bikoff_integrated, 
                                          query.assay = "integrated", 
                                          k.anchor = 5,
                                          dims = 1:100, 
                                          npcs = 100)


# Predict cell types using the anchors from above
Data_predictions <- TransferData(anchorset = Data_query.anchors, refdata = Data_Delile_integrated$Specific_Identities)

# Add the predicted cell types back to the Bikoff metadata
Data_Bikoff_integrated <- AddMetaData(Data_Bikoff_integrated, metadata = Data_predictions)

# examine the results
table(Data_Bikoff_integrated$predicted.id)

# Visualize assigned Bikoff identities with labels
DimPlot(Data_Bikoff_integrated, group.by = "predicted.id", na.value = "black", shuffle = FALSE, label = TRUE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Query label transfer")
ggsave("Graph_UMAP_Bikoff_Labels_FIGURE.pdf", bg = "transparent", device = "pdf", width = 7, height = 7)

DimPlot(Data_Bikoff_integrated, group.by = "predicted.id", na.value = "black", shuffle = TRUE, label = FALSE, repel = TRUE, order = orders, 
        cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Query label transfer")
ggsave("Graph_UMAP_Bikoff_Legend_FIGURE.pdf", bg = "transparent", device = "pdf", width = 7, height = 7)




