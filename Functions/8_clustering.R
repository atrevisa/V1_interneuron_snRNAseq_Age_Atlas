# Run clustering on the code
# AT Updated 3/31/24

# Arguments
# resolution = c(0.5)
# example
# Data_V1 <- clustering(Data_V1, resolution = c(1))

# This is intended to be used on integrated data, so that will be the data used here

clustering <- function(seurat_object, resolution){
  
  # Set default assay to integrated
  DefaultAssay(seurat_object) <- "integrated"
  
  # do not run scale.data if using SCT based normalization methods, which we did for our integration
  
  # First run PCA
  seurat_object <- RunPCA(seurat_object, verbose = TRUE, assay = 'integrated')

  # Function to determine the number of PC's to use in downstream analysis
  find_PC <- function(seurat_object){
    # Determine the percent of variation associated with each PC
    pct1 <- seurat_object[["pca"]]@stdev / sum(seurat_object[["pca"]]@stdev) * 100
    # Calculate the cumulative percents for each PC
    cumu1 <- cumsum(pct1)
    # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
    co11 <- which(cumu1 > 90 & pct1 < 5)[1]
    # Determine the difference between variation of PC and subsequent PC
    co21 <- sort(which((pct1[1:length(pct1) - 1] - pct1[2:length(pct1)]) > 0.1), decreasing = T)[1] + 1
    # Last point where the change of % of variation is more than 0.1%
    # Choose the minimum of these two metrics as the PCs covering the majority of the variation in the data
    pcs1 <- min(co11, co21)
    # Create a df with values
    #plot_df1 <- data.frame(pct = pct1, cumu = cumu1, rank = 1:length(pct1))
    #Elbow plot to visualize
    #ggplot(plot_df1, aes(cumu1, pct1, label = rank, color = rank > pcs1)) + geom_text() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    return(pcs1)
  }
  
  # Actually run the find_PC function on the data and return the number
  pcs = find_PC(seurat_object)
  
  # Run UMAP
  seurat_object <- RunUMAP(seurat_object, dims = 1:pcs, reduction = "pca")
  
  # Determine the K-nearest neighbor graph
  seurat_object <- FindNeighbors(object = seurat_object,dims = 1:pcs, k.param = 50)
  
  # Determine the clusters for various resolutions                                
  seurat_object <- FindClusters(object = seurat_object, resolution = resolution, algorithm = 3)
  
  return(seurat_object)
  
}