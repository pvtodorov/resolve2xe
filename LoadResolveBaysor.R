LoadResolveBaysor <- function(data.dir, fov = 'fov', assay = 'Xenium') {
  data <- Seurat::ReadXenium(
    data.dir = data.dir,
    type = c("centroids", "segmentations"),
  )

  segmentations.data <- list(
    "centroids" = SeuratObject::CreateCentroids(data$centroids),
    "segmentation" = SeuratObject::CreateSegmentation(data$segmentations)
  )
  coords <- SeuratObject::CreateFOV(
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
  
  xenium.obj <- SeuratObject::CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
  
  xenium.obj[["BlankCodeword"]] <- SeuratObject::CreateAssayObject(counts = zeros)
  xenium.obj[["ControlCodeword"]] <- SeuratObject::CreateAssayObject(counts = zeros)
  xenium.obj[["ControlProbe"]] <- SeuratObject::CreateAssayObject(counts = zeros)

  cells = read_csv(paste0(data.dir, 'cells.csv.gz')) %>% 
    distinct(cell_id, .keep_all = TRUE) %>%
    column_to_rownames("cell_id") %>%
    select(cell_area, density, elongation, avg_confidence, nucleus_area)
  xenium.obj = Seurat::AddMetaData(xenium.obj, metadata = cells)

  xenium.obj[[fov]] <- coords
  return(xenium.obj)
}
