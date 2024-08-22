# Find the # of PC's to use
# AT Updated 9/15/23

# EXAMPLE
# pcs <- find_PC(seurat_object)
# returns an integer

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
  
  # Return the number of PC's to use
  return(pcs1)
  
}