# Add metadata about cells for Bikoff samples only 
# AT Updated 5/19/23

# Function to add metadata to Bikoff samples only
# Will only add the metadata specified by the user i.e. region, replicate, and/or age
# Will also add a dummy variable - useful for setting the Ident to all cells
# example: P0LR1 <- AddMetadata(P0LR1)
# If metadata is already there it will not hurt anything
# Metadata added = "region", "replicate", "age", "dummy_variable"

# example: QC_Data = Add_metadata(QC_Data)

Add_metadata <- function(Seurat_object) {
  
  # The format for my data is to have a column called "orig.ident" that has the structure AGEREGIONREPLICATE like this: P0LR1
  
  # First save metadata in a dataframe
  metadata <- Seurat_object@meta.data
  
  # Extract 3rd character - this is the region
  metadata$region = str_sub(metadata$orig.ident,-3,-3)
  
  # Extract last character - this is the replicate
  metadata$replicate = str_sub(metadata$orig.ident,-1, -1)
  
  # Extract first numerical data - this is the age
  metadata$age <- as.numeric(gsub(".*?([0-9]+).*", "\\1", metadata$orig.ident))
  
  # Add a dummy variable - if you want to look at all cells for a graph or something I have found this helpful
  metadata$dummy_variable = "dummy_variable"
  
  # Distinguish between Bikoff and Delile - only relevant for positive selection method
  # Adds source data, bikoff or delile in a new column called source
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
  
  # Add new metadata back to the seurat object
  Seurat_object<- AddMetaData(object = Seurat_object, metadata = metadata)
  
  # Return modified Seurat object
  return(Seurat_object)
  
}
