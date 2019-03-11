#seuratPlotly.R: Functions replacing the ggplot2-based ploting functions of Seurat with those of Plot.ly
#' Plot dimensional reduction for a Seurat object
#'
#' Create a scatterplot of a given dimensional reduction set for a scRNA-seq data object,
#' coloring points by the given grouping variable
#'
#' @param object scRNA-seq object
#' @param grouping Variable by which to group cells. Currently only works with the current ident and column names from meta.data (default: ident)
#' @param reduction_use Dimensional reduction to display (default: tsne)
#' @param dim_1 Dimension to display on the x-axis (default: 1)
#' @param dim_2 Dimension to display on the y-axis (default: 2)
#' @param do.label Add a label showing thr group name to the graph (default: FALSE)
#' @param label.size Label font size (default: 12)
#' @param show.arrow Offset the position of the labels and instead point to each group with an arrow (default: FALSE)
#' @param label.color Color for label border and arrow.  Need hex value. (default = '000000')
#' @param pt.size Size of the points in pixels (default: 2)
#' @param pt_shape Shape to use for the points (default: circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param palette_use Color palette to use.  Must be a palette available in the Paletteer package
#' @param plot_height Plot height in pixels (default: 900)
#' @param plot_width Plot width in pixels (default: 900)
#' @param legend Display legend? (default: TRUE)
#' @param legend_font_size Legend font size (default: 12)
#' @param pt_info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param return Return the plot object instead of displaying it (default: FALSE)
#'
#' @import dplyr group_by summarise
#' @importFrom plotly plot_ly layout
#'
#' @return plotly object
#' @export
#'
#' @examples
#' object <- RunTSNE(object)
#' DimPlotly(object, grouping = 'ident', pt.size = 4, opacity = 0.5, plot_title = "Test Plot", reduction_use = "tsne")
DimPlotly <- function(object,
                      grouping = "ident",
                      do.label = FALSE,
                      label.size = 12,
                      show.arrow = FALSE,
                      label.color = "000000",
                      return = FALSE,
                      pt.size = 4,
                      pt_shape = "circle",
                      opacity = 0.75,
                      reduction_use = "tsne",
                      dim_1 = 1,
                      dim_2 = 2,
                      palette_use = "default_ucscgb",
                      plot_height = 900,
                      plot_width = 900,
                      plot_title = NULL,
                      pt_info = NULL,
                      legend = TRUE,
                      legend_font_size = 12){

  df <- PrepDf(object,
               reduction_use,
               dim_1 = dim_1,
               dim_2 = dim_2,
               grouping = grouping)

  df <- PrepInfo(object = object,
                 pt_info = pt_info,
                 df = df)

  pal <- PrepPalette(df = df,
                     palette_use = palette_use)

  if (do.label) {
    df %>%
      group_by(ident) %>%
      centers <- summarise(
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

  if (is.null(plot_title)){
    plot_title <- reduction_use
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               color = ~ident,
               colors = pal,
               marker = list(
                 symbol = pt_shape,
                 size = pt.size,
                 opacity = opacity,
                 mode = "markers"
               ),
               width = plot_width,
               height = plot_height,
               type = "scattergl",
               mode = "markers",
               showlegend = legend,
               hoverinfo = "text",
               text = ~meta.info
  ) %>%
    layout(
      title = plot_title,
      xaxis = list(title = glue("{reduction_use}_{dim_1}")),
      yaxis = list(title = glue("{reduction_use}_{dim_2}")),
      annotations = labels
    )

  p <- p %>% layout(legend = list(
    font = list(
      size = legend_font_size)
    )
  )

  if (isTRUE(return)){
    return(p)
  } else {
    p
  }
}
