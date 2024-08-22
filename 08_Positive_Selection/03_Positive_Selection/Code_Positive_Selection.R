# Label Transfer to soupx corrected postnatal data
# AT Updated 6/11/24
# Positive selection on soupx data that has already undergone negative selection

#### Setup #####




# Set up default directory
setwd("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/03_PositiveSelectionDelile/3_PostnatalV1Selection/soupx/LabelTransfer/")

# Install necessary packages if not already installed
# source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/0_IstallPackages.R")

# load packages if necessary
source("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/1_LoadLibraries.R")

# Import custom functions
sourceDirectory("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/Alex_Seurat_functions/", modifiedOnly=FALSE)




#### Load data ####





# Load the Bikoff embryonic data analyzed alone that have had cell identities labeled by the label transfer method
Data_reference <- readRDS("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/03_PositiveSelectionDelile/2_Delile_Bikoff_Embryonic_together/Data_Bikoff_LabelTransfer.rds")

# VIEW REFERENCE DATA QUICKLY
#view colnames
colnames(Data_reference@meta.data)
# Specifcy the order you want the labels to appear on the graph
orders = rev(c("Query", "Null", "Progenitors", "Mesoderm", "Blood", "NC","DRG", "dI1", "dI2","dI3", "dI4", "dI5", "dI6","V0","V1", "V2", "MN", "V3"))
# Specifcy the colors for each cluster, from parade pallette
color_list = c("Query" = "black", "Null" = "#cdcd69", "Progenitors" = "#ff37cd", "Mesoderm" = "#379bcd", "Blood" = "#ff6969", "NC" = "#cd69cd","DRG" = "#cd6937", "dI1" = "#ff3769", "dI2" = "#9b9bff","dI3" = "#69ff37", "dI4" = "#ff6937", "dI5" = "#9b69ff", "dI6" = "#69ffff","V0" = "#37ff69", "V1" = "#ff9b37", "V2" = "#37ffcd", "MN" = "#379b37", "V3" = "#37379b")
# Quickly visualize the reference dataset
DimPlot(Data_reference, group.by = "predicted.id", na.value = "black", shuffle = FALSE, label = FALSE, repel = TRUE, order = orders, cols = color_list, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Reference")
DimPlot(Data_reference, group.by = "source", shuffle = FALSE, label = FALSE, repel = TRUE, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Reference")
DimPlot(Data_reference, group.by = "Specific_Identities", shuffle = FALSE, label = FALSE, repel = TRUE, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Reference")
DimPlot(Data_reference, group.by = "Final_identities", shuffle = FALSE, label = FALSE, repel = TRUE, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Reference")

# Load soupX data that has had QC run on it, been log normalized, and undergone negative selection with and without integration
Data_query <- readRDS("/media/ResearchHome/bikoffgrp/home/atrevisa/RNA/Analysis/02_NegativeSelection/soupx/4_Integration/Data_V1_negSelection3.rds")
Data_query # 89545 (should match the number of cells in the negative selection analysis)
DimPlot(Data_query, group.by = "seurat_clusters", shuffle = FALSE, label = FALSE, repel = TRUE, pt.size = 0.5, label.size = 6) + NoAxes() + theme(aspect.ratio=1) + ggtitle("Query")




#### Label Transfer ####




# Find the query anchors
Data_query.anchors <- FindTransferAnchors(reference = Data_reference, 
                                          reference.assay = "RNA", 
                                          query = Data_query,
                                          query.assay = "RNA")     

# predict cell types using the anchors above
Data_predictions <- TransferData(anchorset = Data_query.anchors, 
                                 refdata = Data_reference$predicted.id)

# Add the predicted cell types back to the Bikoff metadata
Data_query <- AddMetaData(Data_query, metadata = Data_predictions)

# examine the results
table(Data_query$predicted.id, useNA = "always")
# Most of the cells, 89494, were labeled as V1

# Export the V1 cell identities to a csv
Idents(Data_query) <- "predicted.id"
levels(Data_query)
Cells_V1_pos <- WhichCells(Data_query, idents = "V1")
length(Cells_V1_pos)
Cells_V1_pos
write.csv(Cells_V1_pos, "Cells_pos.csv", row.names = TRUE)

