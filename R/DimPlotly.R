#seuratPlotly.R: Functions replacing the ggplot2-based ploting functions of Seurat with those of Plot.ly
#' Plot dimensional reduction for a Seurat object
#'
#' Create a scatterplot of a given dimensional reduction set for a scRNA-seq data object,
#' coloring points by the given grouping variable
#'
#' @param object scRNA-seq object
#' @param group.by Variable by which to group cells. Currently only works with the current ident and column names from meta.data (default: ident)
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param do.label Add a label showing thr group name to the graph (default: FALSE)
#' @param label.size Label font size (default: 12)
#' @param show.arrow Offset the position of the labels and instead point to each group with an arrow (default: FALSE)
#' @param label.color Color for label border and arrow.  Need hex value. (default = '000000')
#' @param pt.size Size of the points in pixels (default: 2)
#' @param pt.shape Shape to use for the points (default: circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param palette.use Color palette to use.  Must be a palette available in the Paletteer package
#' @param plot.height Plot height in pixels (default: 900)
#' @param plot.width Plot width in pixels (default: 900)
#' @param legend Display legend? (default: TRUE)
#' @param legend.font.size Legend font size (default: 12)
#' @param pt.info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param do.return Return the plot object instead of displaying it (default: FALSE)
#'
#' @import dplyr
#' @importFrom plotly plot_ly layout
#'
#' @return plotly object
#' @export
#'
#' @examples
#' object <- RunTSNE(object)
#' DimPlotly(object, group.by = 'ident', pt.size = 4, opacity = 0.5, plot.title = "Test Plot", reduction.use = "tsne")
DimPlotly <- function(object,
                      group.by = "ident",
                      do.label = FALSE,
                      label.size = 12,
                      show.arrow = FALSE,
                      label.color = "000000",
                      do.return = FALSE,
                      pt.size = 4,
                      pt.shape = "circle",
                      opacity = 0.75,
                      reduction.use = "tsne",
                      dim.1 = 1,
                      dim.2 = 2,
                      palette.use = "default_ucscgb",
                      plot.height = 900,
                      plot.width = 900,
                      plot.title = NULL,
                      pt.info = NULL,
                      legend = TRUE,
                      legend.font.size = 12){

  df <- PrepDf(object,
               reduction.use,
               dim.1 = dim.1,
               dim.2 = dim.2,
               group.by = group.by)

  df <- PrepInfo(object = object,
                 pt.info = pt.info,
                 df = df)

  pal <- PrepPalette(df = df,
                     palette.use = palette.use)

  if (do.label) {
    df %>%
      dplyr::group_by(ident) %>%
      centers <- summarize(
        x = median(x = as.double(x)),
        y = median(x = as.double(y))
      )
      labels <- list(x = centers$x,
                     y = centers$y,
                     text = centers$ident,
                     font = list(size = label.size),
                     showarrow = show.arrow,
                     bordercolor = label.color,
                     bgcolor = "FFFFFF",
                     opacity = 1
      )
  } else {
    labels <- NULL
  }

  if (is.null(plot.title)){
    plot.title <- reduction.use
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               color = ~ident,
               colors = pal,
               marker = list(
                 symbol = pt.shape,
                 size = pt.size,
                 opacity = opacity,
                 mode = "markers"
               ),
               width = plot.width,
               height = plot.height,
               type = "scattergl",
               mode = "markers",
               showlegend = legend,
               hoverinfo = "text",
               text = ~meta.info
  ) %>%
    layout(
      title = plot.title,
      xaxis = list(title = glue("{reduction.use}_{dim.1}")),
      yaxis = list(title = glue("{reduction.use}_{dim.2}")),
      annotations = labels
    )

  p <- p %>% layout(legend = list(
    font = list(
      size = legend.font.size)
    )
  )

  if (isTRUE(do.return)){
    return(p)
  } else {
    p
  }
}
