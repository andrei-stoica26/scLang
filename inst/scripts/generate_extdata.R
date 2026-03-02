# A self-contained script to generate the data in inst/extdata
#'
createExtData <- function(){
    if (requireNamespace(c('Matrix', 'qs2', 'Seurat', 'withr'),
                         quietly=TRUE)){

        nGenes <- 500
        nCells <- 200

        countsMat <- withr::with_seed(1, matrix(
            rnbinom(nGenes * nCells, size=5, mu=5),
            nrow = nGenes))

        countsMat[withr::with_seed(1, sample(length(countsMat),
                             0.7 * length(countsMat)))] <- 0
        countsMat <- Matrix::Matrix(countsMat, sparse=TRUE)

        rownames(countsMat) <- paste0("Gene", seq(nGenes))
        colnames(countsMat) <- paste0("Cell", seq(nCells))

        seuratObj <- CreateSeuratObject(countsMat)
        seuratObj <- NormalizeData(seuratObj)
        seuratObj <- FindVariableFeatures(seuratObj)
        seuratObj <- ScaleData(seuratObj)
        seuratObj <- withr::with_seed(1, RunPCA(seuratObj))
        seuratObj <- withr::with_seed(1, RunUMAP(seuratObj, dims=1:30))

        scCol(seuratObj, 'Cluster') <- withr::with_seed(1,
                                                     sample(paste0('Cluster',
                                                                   seq(5)),
                                                            dim(seuratObj)[2],
                                                            replace=TRUE))
        scCol(seuratObj, 'Donor') <- rep('Donor1', dim(seuratObj)[2])
        for (i in seq(5)){
            scCol(seuratObj, 'Donor')[withr::with_seed(1,
                                                    sample(which(scCol(seuratObj, 'Cluster') == paste0('Cluster', i))
                                                           , 15))]<- paste0('Donor', i)
            scCol(seuratObj, 'Donor')[withr::with_seed(1,
                                                    sample(which(scCol(seuratObj,
                                                                       'Cluster') ==
                                                                     paste0('Cluster', i))
                                                           , 15))]<- paste0('Donor', i + 1)
        }
        qs2::qs_save(seuratObj, 'inst/extdata/seuratObj.qs2')
        sceObj <- as.SingleCellExperiment(seuratObj)
        qs2::qs_save(sceObj, 'inst/extdata/sceObj.qs2')
    }
}

createExtData()
