# This is a function that will automatically integrate a bunch of seurat samples for you
# AT Updated 3/18/24

# This function normalizes with SCT
# Default is to use 50 PCs for integration, unless the user specifies otherwise
# Must provide a merged seurat object as input
# Must provide reference datasets. If you do not wish to do so, set references = None. References should be in vector format using the index of the seurat object list generated in the first line of the function

# Example:
# Integrated <- Parallel_integrate(QC_Data, sreferences = c(1,2,3))

# You need to increase the default memory to run this
# options(future.globals.maxSize = 50000 * 1024^2)

Parallel_integrate <- function(seurat_object, ndim = 30, anchors = 10, sreferences = NULL, weight = 100) {
  
  # Split the seurat object into a list of individual experiments
  seurat_object.list <- SplitObject(seurat_object, split.by = "orig.ident")
  
  # SCT
  print("SCT")
  seurat_object.list <- lapply(X = seurat_object.list, FUN = SCTransform, method = "glmGamPoi")

  # Third, select features that are repeatedly variable across datasets for integration
  print("features")
  features <- SelectIntegrationFeatures(object.list = seurat_object.list, nfeatures = 3000)
  
  # Fourth, prepare the SCT list object for integration
  print("Prep")
  seurat_object.list <- PrepSCTIntegration(object.list = seurat_object.list, anchor.features = features)

  # Fifth, PCA
  print("PCA")
  seurat_object.list <- lapply(X = seurat_object.list, FUN = RunPCA, features = features, verbose = FALSE)
 
  # Sixth, find best buddies i.e. integration anchors. This step requires references for large datasets, an optional argumnet for this function
  print("anchors")
  #print(as.character(sreferences))
  seurat_object.anchors <- FindIntegrationAnchors(object.list = seurat_object.list, normalization.method = "SCT",anchor.features = features, dims = 1:ndim, k.anchor = anchors, reference = sreferences, verbose = FALSE)
  
  # Seventh, integrate across conditions
  # The integrated data will be stored in a new Assay "integrate"
  print("integrate")
  seurat_object <- IntegrateData(anchorset = seurat_object.anchors, normalization.method = "SCT", dims = 1:ndim, k.weight = weight, verbose = FALSE)
  
  return(seurat_object)
  
}