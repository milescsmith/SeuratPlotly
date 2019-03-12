#seuratPlotly.R: Functions replacing the ggplot2-based ploting functions of Seurat with those of Plot.ly
#' Plot dimensional reduction for a Seurat object
#'
#' Create a scatterplot of a given dimensional reduction set for a scRNA-seq data object,
#' coloring points by the given grouping variable
#'
#' @param object scRNA-seq object
#' @param grouping Variable by which to group cells. Currently only works with the current ident and column names from meta.data. Default: ident
#' @param reduction_use Dimensional reduction to display. Default: tsne
#' @param dim_1 Dimension to display on the x-axis. Default: 1
#' @param dim_2 Dimension to display on the y-axis. Default: 2
#' @param label Add a label showing thr group name to the graph. Default: FALSE
#' @param label_size Label font size. Default: 12
#' @param show_arrow Offset the position of the labels and instead point to each group with an arrow. Default: FALSE
#' @param label_color Color for label border and arrow.  Need hex value.. Default = '000000'
#' @param pt_size Size of the points in pixels. Default: 2
#' @param pt_shape Shape to use for the points. Default: circle
#' @param opacity Transparency level to use for the points, on a 0-1 scale. Default: 1
#' @param palette_use Color palette to use.  Must be a palette available in the Paletteer package
#' @param plot_height Plot height in pixels. Default: 900
#' @param plot_width Plot width in pixels. Default: 900
#' @param legend Display legend?. Default: TRUE
#' @param legend_font_size Legend font size. Default: 12
#' @param pt_info Meta.data columns to add to the hoverinfo popup.. Default: ident
#' @param return Return the plot object instead of displaying it. Default: FALSE
#' @param plot_title
#'
#' @importFrom dplyr group_by summarise
#' @importFrom plotly plot_ly layout
#'
#' @return plotly object
#' @export
#'
#' @examples
#' object <- RunTSNE(object
#' DimPlotly(object, grouping = 'ident', pt_size = 4, opacity = 0.5, plot_title = "Test Plot", reduction_use = "tsne")
DimPlotly <- function(object,
                      grouping = "ident",
                      label = FALSE,
                      label_size = 12,
                      show_arrow = FALSE,
                      label_color = "000000",
                      return = FALSE,
                      pt_size = 4,
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

  df <- PrepDr(object,
               reduction_use,
               dim_1 = dim_1,
               dim_2 = dim_2,
               grouping = grouping)

  pal <- PrepPalette(df = df,
                     palette_use = palette_use)

  if (label) {
    df %>%
      group_by(ident) %>%
      centers <- summarise(
        x = median(x = as.double(x)),
        y = median(x = as.double(y))
      )
      labels <- list(x = centers$x,
                     y = centers$y,
                     text = centers$ident,
                     font = list(size = label_size),
                     showarrow = show_arrow,
                     bordercolor = label_color,
                     bgcolor = "FFFFFF",
                     opacity = 1
      )
  } else {
    labels <- NULL
  }

  md <- GetFeatureValues(object = object,
                         features = c(pt_info, "ident")) %>%
    mutate_at(vars(-cell),
              list(~paste0('</br> ', substitute(.), ": ", .))) %>%
    unite(info, -cell)

  df %<>% inner_join(md)

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
                 size = pt_size,
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
