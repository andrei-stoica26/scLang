test_that("dimredNames works", {
    expect_identical(dimredNames(sceObj), c('PCA', 'UMAP'))
    expect_identical(dimredNames(seuratObj), c('pca', 'umap'))
    expect_error(dimredNames(c(1, 2, 3)))
})

test_that("metadataDF works", {
    v <- metadataDF(seuratObj)
    w <- metadataDF(sceObj)
    expect_equal(dim(v), c(200, 5))
    expect_equal(dim(w), c(200, 6))
    metadataDF(sceObj) <- v
    expect_identical(v, metadataDF(sceObj))
    expect_error(metadataDF(c(1, 2, 3)))
})

test_that("metadataNames works", {
    cols <- c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Cluster', 'Donor')
    expect_equal(metadataNames(seuratObj), cols)
    expect_equal(metadataNames(sceObj), c(cols, 'ident'))
    expect_error(metadataNames(c(1, 2, 3)))
})

test_that("scCol works", {
    scCol(seuratObj, 'featCountRatio') <-
        scCol(seuratObj, 'nFeature_RNA') / scCol(seuratObj, 'nCount_RNA')
    expect_equal(mean(scCol(seuratObj, 'featCountRatio')), 0.1934944,
                 tolerance=0.0001)
    scCol(sceObj, 'featCountRatio') <-
        scCol(sceObj, 'nFeature_RNA') / scCol(sceObj, 'nCount_RNA')
    expect_equal(mean(scCol(sceObj, 'featCountRatio')), 0.1934944,
                 tolerance=0.0001)
    expect_error(scCol(c(1, 2, 3)))
    expect_error(scCol(seuratObj, 'noCol'))
    expect_error(scCol(sceObj, 'noCol'))
})

test_that("counting functions work", {
    expect_equal(scColCounts(seuratObj, 'Donor')['Donor1'],
                 setNames(125, 'Donor1'))
    expect_equal(scColCounts(sceObj, 'Donor')['Donor1'],
                 setNames(125, 'Donor1'))
    expect_error(scColCounts(c(1, 2, 3)))

    expect_equal(scColPairCounts(seuratObj, 'Cluster', 'Donor')[1, 3],
                 32)
    expect_equal(scColPairCounts(sceObj, 'Cluster', 'Donor')[1, 3],
                 32)
    expect_error(scColPairCounts(c(1, 2, 3)))

    expect_equal(scColPairPercs(seuratObj, 'Cluster', 'Donor')[1, 4],
                 25.6)
    expect_equal(scColPairPercs(sceObj, 'Cluster', 'Donor')[1, 4],
                 25.6)
    expect_error(scColPairPercs(c(1, 2, 3)))
})

test_that("scGeneExp works", {
    expect_equal(length(scGeneExp(sceObj, 'Gene474')), 200)
    expect_equal(length(scGeneExp(seuratObj, 'Gene474')), 200)
    expect_equal(mean(scGeneExp(sceObj, 'Gene474')), 1.311307,
                 tolerance=0.0001)
    expect_equal(mean(scGeneExp(seuratObj, 'Gene474')), 1.311307,
                 tolerance=0.0001)
    mat <- scExpMat(sceObj)
    expect_equal(mean(scGeneExp(mat, 'Gene474')), 1.311307,
                 tolerance=0.0001)
    expect_error(scGeneExp(sceObj, 'missingGene'))
    expect_error(scGeneExp(seuratObj, 'missingGene'))
    expect_error(scGeneExp(c(1, 2, 3)))
})

test_that("scExpMat works", {
    v <- scExpMat(sceObj)
    w <- scExpMat(seuratObj)
    expect_identical(v, w)
    mat <- scExpMat(v)
    expect_identical(mat, v)
    expect_equal(dim(v), c(500, 200))

    v <- scExpMat(sceObj, slot='counts')
    w <- scExpMat(seuratObj, slot='counts')
    expect_identical(v, w)
    expect_error(scExpMat(sceObj, slot='missingSlot'))
    expect_error(scExpMat(seuratObj, slot='missingSlot'))
    expect_equal(dim(v), c(500, 200))

    v <- scExpMat(sceObj, genes=c('Gene285', 'Gene287', 'Gene291'))
    w <- scExpMat(seuratObj, genes=c('Gene285', 'Gene287', 'Gene291'))
    expect_equal(dim(v), c(3, 200))
    expect_identical(v, w)
    mat <- scExpMat(v)
    expect_identical(mat, v)
    expect_error(scExpMat(seuratObj, genes=c('Gene285', 'noGene')))
    expect_error(scExpMat(c(1, 2, 3)))
})

test_that("scExpMat works for sparse matrices", {
    v <- scExpMat(sceObj, densify=FALSE)
    w <- scExpMat(seuratObj, densify=FALSE)
    expect_identical(v, w)
    expect_equal(dim(v), c(500, 200))
})

test_that("scPCAMat works", {
    v <- scPCAMat(sceObj)
    w <- scPCAMat(seuratObj)
    expect_identical(v, w)
    expect_equal(dim(v), c(200, 50))
    expect_error(scPCAMat(c(1, 2, 3)))
})

test_that("scUMAPMat works", {
    v <- scUMAPMat(sceObj)
    w <- scUMAPMat(seuratObj)
    expect_identical(v, w)
    expect_equal(dim(v), c(200, 2))
    expect_error(scUMAPMat(c(1, 2, 3)))
})

test_that("visualization functions work", {
    p <- dimPlot(sceObj, groupBy='Donor')
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)

    p <- featurePlot(sceObj, 'Gene289')
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)

    p <- featurePlot(sceObj, 'nCount_RNA')
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)

    p <- violinPlot(sceObj, 'Gene289')
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)

    p <- violinPlot(sceObj, 'nCount_RNA')
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)

    expect_error(featurePlot(sceObj, 'missingFeat'))
    expect_error(violinPlot(sceObj, 'missingFeat'))
})
