# Automatically load data and name Seurat objects, calculate some stats 
# AT Updated 12/14/23

# Description: This function automatically detects folders in a directory ending in *_cellranger
# Then it makes individual seurat objects out of them. 
# You can choose to use the orginal filtered_feature_barcodes or choose to use the soupX corrected data (strainedCounts), but the original one is the default
# This choice is dictated by the optional argument data_type which is set to a string which may either be "filtered_feature_bc_matrix" or "strainedCounts"
# The individual Seurat objects are then merged into one large Seurat object
# The name of the file is added to a metadata column called orig.ident
# This function calculates the percent_mt for each cell and adds that to a column in the metadata
# This function calculates the "novelty score" and adds it in a column in the metadata
# Novelty score is the log10 of the number of genes detected per cell and the log10 of the number of UMIs per cell, then divide the log10 number of genes by the log10 number of UMIs
# This function also caluclates the percent_ribo for each cell and adds it to a column of the same name
# Note that genes that appear in less than 10 cells will be removed and cells with less than 1000 genes will be removed by default. These parameters can be customized

# Arguments: folder_name (required), min_cells (optional, default is 10), min_features (optional, default = 1000), data_type (required, default is "filtered_feature_bc_matrix")

# EXAMPLE: note that this function will not work unless you set it equal to a variable
# folder_name = "/mnt/storage1/AlexTrev/Alex_Seurat_functions/"
# Please note that the folder name must end in a "/"
# Raw_data = Load_and_merge("folder_name")

Load_and_merge <- function(folder_name, min_cells = 10, min_features = 1000, data_type = "filtered_feature_bc_matrix") {
  
  # filepaths is a string used to get the filepaths for all cellranger files in the given folder
  filepaths = paste0(folder_name, "*_cellranger")

  # Filepaths are the actual paths to all the cellranger files to be analyzed
  files <- Sys.glob(filepaths)

  # Create an empty list to put seurat objects in
  SeuratList <- vector("list")

  # Which form of the data are you using - choose the original filtered_feature_bc_matrix or soupX which is "strainedCounts"
  directory = paste0("/outs/", data_type)

  # For every cellranger file do the following:
  for (file in c(files)) {

    # Read in 10X filtered feature barcode matrix
    seurat_data <- Read10X(data.dir = paste0(file, directory), unique.features = TRUE, strip.suffix = TRUE)

    # Extract the "orig.ident" from the filename
    file_names_full = basename(file)
    file_names_final = strsplit(file_names_full, "_")[[1]][1]
	
	# Create seurat object
    SeuratList[[file_names_final]] <- CreateSeuratObject(counts = seurat_data, min.cells = min_cells, min.features = min_features, project = file_names_final)

  }

  # Create a merged Seurat object.
  if(length(SeuratList) == 1) {
    Seurat_raw_data <- SeuratList[[1]]
    } else {
      Seurat_raw_data <- merge(SeuratList[[1]], y = c(SeuratList[2:length(SeuratList)]), add.cell.ids = c(names(SeuratList)))
    }

  # Calculate percent_mt for all the cells
  Seurat_raw_data[["percent_mt"]] <- PercentageFeatureSet(object = Seurat_raw_data, pattern = "^mt-")
  
  # Calculate the novelty score for all the cells
  Seurat_raw_data[["novelty_score"]] <- log10(Seurat_raw_data$nFeature_RNA) / log10(Seurat_raw_data$nCount_RNA)

  # Calculate the log of nFeatures and nCounts
  Seurat_raw_data[["lognFeature_RNA"]] <- log10(Seurat_raw_data$nFeature_RNA)
  Seurat_raw_data[["lognCount_RNA"]] <- log10(Seurat_raw_data$nCount_RNA)
  
  # Calculate the percent of genes coming from the ribosome
  Seurat_raw_data[["percent_ribo"]] <- PercentageFeatureSet(object = Seurat_raw_data, pattern = "^Rp")
  # grep("^Rp",rownames(Data_norm@assays$RNA@counts),value = TRUE) # check all the ribo gene names
  
  # Return the Seurat object
  return(Seurat_raw_data)

}