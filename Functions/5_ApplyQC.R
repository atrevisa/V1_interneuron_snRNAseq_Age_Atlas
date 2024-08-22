# Apply QC to data
# AT Updated 11/14/23

# Description: This function will run QC on a seurat object containing many merged samples. 
# QC parameters applied:
# Any cell with % mt reads > than max_mt (default is 2.5)
# Will only use top 5000 cells when sorted by # of Features, may change this number if you want using the featres_cutoff argument but the default is 5000
# Will remove any cells outside of the median +- 3 MADs for features and counts, may change cutoff
# Will remove cells with novelty scores < 3 MADs
# Only run once on your data or it will continue cutting off cells that fall within these parameters each time you do so

# Example
# QC_Data <- ApplyQC(Raw_Data)


ApplyQC <- function(seurat_object, max_mt = 2.5, low = 3, high = 3, features_cutoff = 5000) {
    
  # Remove cells with more than max_mt mt reads - need higher cutoff for single cell 
  seurat_object_filtered <- subset(seurat_object, subset = percent_mt < max_mt)
  
  # Split the dataset into a list of seurat objects
  seurat_list <- SplitObject(seurat_object_filtered, split.by = "orig.ident")
  
  # Create function to select the top 5,000 cells with the greatest number of features (genes)
  topFeatures <- function(a_seurat_object){
    qc = a_seurat_object@meta.data
    qc <- qc %>% arrange(desc(nFeature_RNA))
    qc_top = head(qc,features_cutoff)
    qc_cells = rownames(qc_top)
    a_seurat_object = subset(x = a_seurat_object, cells = c(qc_cells))
    return(a_seurat_object)
  }
  
  # Iterate through each individual sample and get rid of high and low outliers within each one to account for batch variation
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    
    # Top n cells
    x <- topFeatures(x)
    
    # Remove cells that are +-2 MAD's above or  below the median number of features and counts
    upper_Counts <- median(x@meta.data$nCount_RNA) + high*mad(x@meta.data$nCount_RNA) # upper limit
    lower_Counts <- median(x@meta.data$nCount_RNA) - low*mad(x@meta.data$nCount_RNA) # lower limit
    upper_Features <- median(x@meta.data$nFeature_RNA) + high*mad(x@meta.data$nFeature_RNA) # upper limit
    lower_Features <- median(x@meta.data$nFeature_RNA) - low*mad(x@meta.data$nFeature_RNA) # upper limit
    lower_novelty <- median(x@meta.data$novelty_score) - low*mad(x@meta.data$novelty_score)
    
    # actually subset
    x <- subset(x, subset = nFeature_RNA > lower_Features & nFeature_RNA < upper_Features & nCount_RNA > lower_Counts & nCount_RNA < upper_Counts & novelty_score > lower_novelty)
    
  })
  
  # Re-merge big.list back together again
  # Create a merged Seurat object.
  if(length(seurat_list) == 1) { 
    seurat_object_filtered <- seurat_list[[1]]
  } else {
    seurat_object_filtered <- merge(seurat_list[[1]], y = c(seurat_list[2:length(seurat_list)]))
  }
  
  # Retrun modified seurat object
  return(seurat_object_filtered)
  
}
