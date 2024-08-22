# QC graphs
# AT Updated 11/14/23

# qc_visuals outputs standard QC metrics of a seurat object
# Cells per sample, nFeature, nCount, percent_mt, and UMI/gene relationship
# this function will output the graphs into a PDF
# EXAMPLE: qc_list = qc_visuals(seurat_object_name_here)
# use qc_list[[1]], qc_list[[2]], etc. to access each graph individually

qc_visuals <- function(seurat_object, filename = NULL, spliton = "orig.ident") {
  
  # Set identity level of the Seurat object, default is orig.ident
  Idents(seurat_object) <- spliton
  
  # QC stats, broken down by orig.ident
  QC_plot <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3, pt.size = 0)
  
  #Extract metadata from seurat object
  metadata <- seurat_object@meta.data
  
  # Visualize the number of cell counts per sample
  cellnumber <- metadata %>% ggplot(aes_string(x=spliton, fill=spliton)) + geom_bar() + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(plot.title = element_text(hjust=0.5, face="bold")) + ggtitle("Cells per sample")
  
  # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
  # UMI_gene_cor <- metadata %>% ggplot(aes(x=lognCount_RNA, y=lognFeature_RNA, color=percent_mt)) + geom_point() + scale_colour_gradient(low = "gray90", high = "black") + stat_smooth(method=lm) + scale_x_log10() + scale_y_log10() + theme_classic() + facet_wrap(spliton)
  UMI_gene_cor <- metadata %>% ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=novelty_score)) + geom_point() + scale_colour_gradient(low = "gray90", high = "black") + theme_classic() + facet_wrap(spliton)

  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  count_gene_cor <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  # View novelty, the log10GenesPerUMI
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
  Nov_score <- metadata %>%
    ggplot(aes_string(x="novelty_score", color = spliton)) +
    geom_density(alpha = 0.2) +
    theme_classic() 
  
  # Create an empty list to put all the QC plots in
  plot_list = list(QC_plot, cellnumber, UMI_gene_cor, count_gene_cor, Nov_score)
  
  # Export each plot on a separate page of a pdf file
  if(!is.null(filename)) {
    pdf(filename)
    pdf.options(width = 9, height = 7)
    for (i in 1:length(plot_list)){
      print(plot_list[[i]])
    }
    dev.off()
  }
  
  # Return list of plots
  return(plot_list)
  
}
