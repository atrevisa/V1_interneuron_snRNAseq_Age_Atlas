# 6_Norm
# Normalize and scale your data
# AT Updated 11/16/23

# This is a function that will normalize and scale your data
# All of these changes will be stored in the "RNA" assay
# It will log-normalize the slot = "counts" (raw counts) data and put the log normalized values into the "data" slot of the "RNA" assay
# Log normalization is done on a per cell basis
# This function will also scale all of the genes and put the scaled values in the scale.data slot of the RNA assay
# Note that you can optionally regress out some variables with scale.data, but the default is NULL
# Regression usually takes a long time

# Example

# Check the data before normalizing for comparison
# head(Data_QC[["RNA"]]@counts['Foxp2',]) # look at raw counts, will be integers
# head(Data_QC[["RNA"]]@data['Foxp2',]) # notice that the log normalized data looks exactly the same as the counts right now, because we haven't done the normalization yet
# head(Data_QC[["RNA"]]@scale.data['Foxp2',]) # this will throw an error because we haven't scaled the data yet either

# Now run the function
# Data_norm <- Norm(Data_QC)

# Now see how the values in each slot of the RNA assay have changed
# head(Data_QC[["RNA"]]@counts['Foxp2',]) # These should look exactly the same as before, as they are untouched by the function
# head(Data_QC[["RNA"]]@data['Foxp2',]) # These numbers should now be floats >0
# head(Data_QC[["RNA"]]@scale.data['Foxp2',]) # You will now see values here and they may be negative or positive

Norm <- function(seurat_object, regression = NULL){
  
  # Set default Ident to orig.ident
  Idents(seurat_object) <- "orig.ident"
  
  # Set default assay to RNA
  DefaultAssay(seurat_object) <- "RNA"
  
  # Log Normalize - according to what I have rerad this is on a per cell basis so you do not need to split a merged object
  seurat_object <- NormalizeData(seurat_object, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
  
  # Scale the data on all genes, optionally you may regress out things here
  DefaultAssay(seurat_object) <- "RNA"
  all.genes <- rownames(seurat_object)
  seurat_object <- ScaleData(seurat_object, features = all.genes, assay = "RNA", vars.to.regress = regression)
  
  # Return the modified seurat object
  return(seurat_object)
  
}

