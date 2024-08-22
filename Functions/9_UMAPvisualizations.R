# Standard UMAP visualization output functions
# AT Updated 5/20/23

# Just some common plots you would want to see after clustering
# Plots output into a list

# EXAMPLE: umap_plot_list <- CommonPlots(seurat_object)
# umap_plot_list[[1]], umap_plot_list[[2]], ...


CommonPlots <- function(seurat_object, filename = NULL, Identby = "seurat_clusters", highlights = NULL){
  
  # Plot the clusters on the umap
  Idents(seurat_object) <- Identby
  plot1 = DimPlot(seurat_object,reduction = "umap",label = TRUE,repel = TRUE, label.size = 6, raster = FALSE, cells.highlight = highlights) & NoAxes()
  
  # Plot standard QC metrics over UMAP plot
  Idents(seurat_object) <- Identby
  metrics <-  c("nFeature_RNA", "nCount_RNA", "percent_mt")
  plot2 = FeaturePlot(seurat_object, reduction = "umap", features = metrics, pt.size = 1, label = TRUE, raster = FALSE)
  
  # Violin plot of QC metrics
  Idents(seurat_object) <- Identby
  plot3 = VlnPlot(seurat_object, features = metrics, pt.size = 0)
  
  # Put all the plots in one list
  plot_list = list(plot1, plot2, plot3)
  
  # Export each plot on a separate page of a pdf file
  pdf(filename)
  pdf.options(width = 9, height = 7)
  for (i in 1:length(plot_list)){
    print(plot_list[[i]])
  }
  dev.off()
  
  # Return the list of plots
  return(plot_list)
  
}
