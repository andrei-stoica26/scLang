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
#' @param xLab Label of x axis.
#' @param yLab Label of y axis.
#' @param legendTitle Legend title.
#' @param legendLabs Legend labels.
#' @param legendPos Legend position.
#' @param palette Color palette.
#' @param alpha Opaqueness level.
#' @param labelDF Label data frame.
#' @param labelType Whether to draw a box around labels (option 'boxed') or not
#' (option 'free'). Default is 'free'.
#' @param labelSize Label size.
#' @param labelColor Label color.
#' @param labelRepulsion Repulsion strength between labels.
#' @param labelPull Attraction strength between a text label
#' and its data point.
#' @param maxOverlaps Maximum overlaps.
#' @param boxPadding Amount of padding around box.
#' @param labelPadding Amount of padding around label.
#' @param pointSize Point size.
#' @param pointShape Point shape.
#' @param legendTitleSize Legend title size.
#' @param legendTextSize Legend text size.
#' @param axisTextSize Axis text size.
#' @param axisTitleSize Axis title size.
#' @param xAngle Angle of x axis text.
#' @param vJust Vertical justification in [0, 1].
#' @param theme Plot theme.
#' @param ... Additional arguments passed to \code{centerTitle}.
#'
#' @return \code{NULL}. This function is only used internally for
#' documentation.
#'
#' @keywords internal

documentFun <- function(scObj = NULL,
                        title = NULL,
                        dimred = 'umap',
                        dims = c(1, 2),
                        xLab = 'x',
                        yLab = 'y',
                        legendTitle = 'Legend',
                        legendLabs = c('a', 'b'),
                        legendPos = 'right',
                        palette = 'Spectral',
                        alpha = 1,
                        labelSize = 2.5,
                        labelColor ='black',
                        labelRepulsion = 1,
                        labelPull = 1,
                        maxOverlaps = 50,
                        boxPadding = 0.2,
                        labelPadding = 0.1,
                        pointSize = 0.8,
                        pointShape = 1,
                        legendTitleSize = 10,
                        legendTextSize = 10,
                        axisTextSize = 12,
                        axisTitleSize = 12,
                        xAngle = 45,
                        vJust = 0.6,
                        theme = 'linedraw',
                        ...){}
