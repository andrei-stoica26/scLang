#' Extract the names of available dimensionality reductions from object
#'
#' This function extracts the names of available dimensionality reductions
#' from a Seurat or SingleCellExperiment object.
#'
#' @inheritParams metadataNames
#'
#' @return A character vector.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' dimredNames(sceObj)
#'
#' @export
#'
dimredNames <- function(scObj)
    UseMethod(generic='dimredNames', object=scObj)

#' Extract metadata from object as a data frame
#'
#' This function extracts the metadata from a \code{Seurat} or
#' \code{SingleCellExperiment} object as a data frame.
#'
#' @param scObj A \code{Seurat} or \code{SingleCellExperiment} object.
#' @param value A data frame to replace metadata with.
#'
#' @return A metadata data frame.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' df <- metadataDF(sceObj)
#'
#' @export
#'
metadataDF <- function(scObj)
    UseMethod(generic='metadataDF', object=scObj)

#' @rdname metadataDF
#' @export
#'
`metadataDF<-` <- function(scObj, value)
    UseMethod("metadataDF<-")

#' Extract a metadata/coldata column from object.
#'
#' This function extracts a metadata/coldata column from a \code{Seurat} or
#' \code{SingleCellExperiment} object.
#'
#' @inheritParams metadataNames
#' @param col Column name.
#'
#' @return A vector.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' v <- scCol(sceObj, 'Mutation_Status')
#'
#' @export
#'
scCol <- function(scObj, col)
    UseMethod(generic='scCol', object=scObj)

#' @rdname scCol
#' @export
#'
`scCol<-` <- function(scObj, col, value)
    UseMethod("scCol<-")

#' Extracts a dimensionality reduction matrix from object.
#'
#' This function extracts a dimensionality reduction matrix from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataNames
#' @param dimred Dimensionality reduction.
#'
#' @return A dimensionality reduction matrix.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' pcaMat <- scDimredMat(sceObj, 'pca')
#'
#'
#' @export
#'
scDimredMat <- function(scObj, dimred)
    UseMethod(generic='scDimredMat', object=scObj)

#' Extracts the expression matrix from object.
#'
#' This function extracts an expression matrix from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams scGeneExp
#' @param genes Selected genes. If \code{NULL}, all genes will be retained.
#' @param densify Whether to convert to dense matrix.
#'
#' @return An expression matrix.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' mat <- scExpMat(sceObj, 'counts')
#'
#' @export
#'
scExpMat <- function(scObj, slot, genes = NULL, densify = TRUE)
    UseMethod(generic='scExpMat', object=scObj)

#' Extracts the expression of a single gene
#'
#' This function extracts the expression of a single gene from a Seurat,
#' SingleCellExperiment, dgCMatrix or matrix object.
#
#' @param scObj A \code{Seurat}, \code{SingleCellExperiment},
#' \code{dgCMatrix} or \code{matrix} object.
#' @param gene Selected gene.
#' @param slot Gene expression slot. Ignored if \code{scObj} is of class
#' \code{dgCMatrix} or \code{matrix}.
#'
#' @return A gene expression vector.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' v <- scGeneExp(sceObj, 'Gene_0491')
#'
#' @export
#'
scGeneExp <- function(scObj, gene, slot)
    UseMethod(generic='scGeneExp', object=scObj)
