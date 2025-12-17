suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
})


## Dataset
## Take from http://livercellatlas.org/ human single-nuclei and spatial datasetts


## Look at spatial data

object.names <- c("H35_1", "H35_2", "H36", "H37", "H38")

spatial.rna.objects <- list()
for (object.name in object.names) {
  object <- readRDS(paste0("/NAS/qscb2025/Guilliams_2022_livercellatlas/data/liver/spatial/", object.name, ".rds"))
  object2 <- CreateSeuratObject(
    counts = LayerData(object, assay = "Spatial", layer = "counts"),
    meta.data = object@meta.data
  )
  object2 <- NormalizeData(object2)
  spatial.rna.objects[[object.name]] <- object2
}

# Check metadata
head(spatial.rna.merged)

spatial.rna.merged <- merge(spatial.rna.objects[[1]], spatial.rna.objects[2:length(spatial.rna.objects)])
spatial.rna.merged <- SCTransform(spatial.rna.merged, verbose = FALSE)
spatial.rna.merged <- RunPCA(spatial.rna.merged, verbose = FALSE)
spatial.rna.merged <- FindNeighbors(spatial.rna.merged, dims = 1:30)
spatial.rna.merged <- RunUMAP(spatial.rna.merged, dims = 1:30)

# Exercise 1: Does this require integration?
# If so show implement it
DimPlot(spatial.rna.merged, group.by = "sample_name")
# YOUR CODE HERE


# Exercise 2: Does this require integration?
# If so show implement it
rna.merged <- readRDS("/NAS/qscb2025/Guilliams_2022_livercellatlas/data/liver/rna/GSE192740_rna.rds")

# check metadata
head(rna.merged)

# Perform steps
# YOUR CODE HERE

DimPlot(rna.merged, group.by = "sample")
# YOUR CODE HERE


# Exercise 3: Subset H37 patient and predict the celltypes in H38 patient
# YOUR CODE HERE


# Exercise 4: You were given the objects pre-processed (to some extent)
# Go to http://livercellatlas.org/ or search GSE192740 to download the raw
# counts file and metadata and create the seurat objects

# Exercise 5: Integrate the mouse and the human datasets
