test_that("compatibility functions work", {
    expect_equal(metadataNames(sceObj), c('Mutation_Status',
                                          'Cell_Cycle',
                                          'Treatment',
                                          'sizeFactor',
                                          'Cluster',
                                          'Donor'))

    expect_error(metadataNames(c(1, 2, 3)))
    expect_equal(length(scCol(sceObj, 'Cell_Cycle')), 200)

    expect_equal(scColCounts(sceObj, 'Cell_Cycle')['G2M'],
                 setNames(56, 'G2M'))
    expect_equal(scColPairCounts(sceObj, 'Mutation_Status', 'Cell_Cycle')[1, 3],
                 24)

    expect_equal(length(scGeneExp(sceObj, 'Gene_0480')), 200)
    expect_equal(dim(scExpMat(sceObj)), c(500, 200))

    v <- scPCAMat(sceObj)
    w <- scPCAMat(seuratObj)
    colnames(v) <- paste0('PC_', seq(50))
    expect_equal(v, w)

    v <- scUMAPMat(sceObj)
    w <- scUMAPMat(seuratObj)
    colnames(v) <- paste0('UMAP_', seq(2))
    expect_equal(v, w)
})

test_that("visualization functions work", {
    p <- dimPlot(sceObj, groupBy='Donor')
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
    p <- featurePlot(sceObj, 'Gene_0289')
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
    p <- violinPlot(sceObj, 'Gene_0289')
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
})
