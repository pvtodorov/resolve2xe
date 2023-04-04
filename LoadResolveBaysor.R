library(Seurat)

LoadResolveBaysor <- function(data.dir, fov = 'fov', assay = 'Xenium') {
  data <- ReadXenium(
    data.dir = data.dir,
    type = c("centroids", "segmentations"),
  )

  segmentations.data <- list(
    "centroids" = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations)
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = assay
  )
  counts = data$matrix
  zeros <- Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = dim(counts), 
                          dimnames = dimnames(counts))
  zeros <- as(zeros, "dgCMatrix")
  data$matrix = list()
  data$matrix[["Gene Expression"]] = counts
  data$matrix[["BlankCodeword"]] = zeros
  data$matrix[["ControlCodeword"]] = zeros
  data$matrix[["ControlProbe"]] = zeros
  
  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
  
  xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = zeros)
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = zeros)
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = zeros)

  xenium.obj[[fov]] <- coords
  return(xenium.obj)
}
