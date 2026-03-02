#' Internal function used for documenting other functions
#'
#' This function is used internally to help document other functions.
#'
#' @inheritParams scDimredMat
#' @param title Plot title.
#' @param groupBy Grouping variable. Must exist in the metadata/coldata of
#' the single-cell expression object.
#' @param dims A numeric vector of size 2 representing the dimensions selected
#' for the plot.
#' @param xLab x axis label.
#' @param yLab y axis label.
#' @param legendTitle Legend title.
#' @param legendLabs Legend labels.
#' @param legendPos Legend position.
#' @param palette Color palette.
#' @param alpha Opaqueness level.
#' @param pointSize Point size.
#' @param legendTitleSize Legend title size.
#' @param legendTextSize Legend text size.
#' @param axisTextSize Axis text size.
#' @param axisTitleSize Axis title size.
#' @param ... Additional arguments passed to \code{henna::centerTitle}.
#'
#' @return \code{NULL}. This function is only used internally for
#' documentation.
#'
#' @keywords internal

documentFun <- function(scObj = NULL,
                        title = NULL,
                        groupBy = NULL,
                        dimred = 'umap',
                        dims = c(1, 2),
                        xLab = 'x',
                        yLab = 'y',
                        legendTitle = 'Legend',
                        legendLabs = c('a', 'b'),
                        legendPos = 'right',
                        palette = 'Spectral',
                        alpha = 1,
                        pointSize = 0.8,
                        legendTitleSize = 10,
                        legendTextSize = 10,
                        axisTextSize = 12,
                        axisTitleSize = 12,
                        ...){}
