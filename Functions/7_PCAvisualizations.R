# Plot PCA plots 
# AT Updated 5/22/23

# You must run parallel integrate before running this function
# Parallel integrated with run RunPCA on your sample
# PCA_Plots will graph the PCA plots for you and output them in a list

# EXAMPLE: PCAplot_list <- PCA_plots(seuratobject, color = "age")
# View plots using PCAplot_list[[1]], PCAplot_list[[2]], PCAplot_list[[3]], PCAplot_list[[4]], etc
# There should be __ plots total
# The color argument can be any metadata variable

PCA_plots <- function(seurat_object, color = "orig.ident", filename = NULL){
  
  # Plot PCA to check
  plot1 = PCAPlot(seurat_object,group.by = color, shuffle = TRUE)

  # Plot PCA
  plot2 = DimPlot(seurat_object, reduction = "pca", shuffle = TRUE, group.by = color)
  
  # Look at the main genes contributing to the first 2 PC's
  plot3 = VizDimLoadings(seurat_object, dims = 1:2, reduction = "pca")
  
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
    plot_df1 <- data.frame(pct = pct1, cumu = cumu1, rank = 1:length(pct1))
    #Elbow plot to visualize
    plot5 = ggplot(plot_df1, aes(cumu1, pct1, label = rank, color = rank > pcs1)) + geom_text() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    return(plot5)
  }
  
  # Elbow plot
  plot4 = find_PC(seurat_object)
  
  # Put all the PCA plots in a list
  plot_list = list(plot1, plot2, plot3, plot4)
  
  # Export each plot on a separate page of a pdf file
  if(!is.null(filename)) {
    pdf(filename)
    for (i in 1:length(plot_list)){
      #pdf(filename)
      print(plot_list[[i]])
    }
    dev.off()
  }
  
  # Return lists of plots
  return(plot_list)
  
}