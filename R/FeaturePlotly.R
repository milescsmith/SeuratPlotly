#' FeaturePlotly
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' object@@dr slot
#'
#' @param object Seurat object
#' @param feature Variable to display. Currently only works with gene names
#' @param reduction Dimensional reduction to display. Default: tsne
#' @param dim_1 Dimension to display on the x-axis. Default: 1
#' @param dim_2 Dimension to display on the y-axis. Default: 2
#' @param pt_scale Factor by which to multiply the size of the points. Default: 5
#' @param pt_shape Shape to use for the points. Default = circle
#' @param opacity Transparency level to use for the points, on a 0-1 scale. Default: 1
#' @param colors_use Color palette to use.  Palettes from RColorBrewer and viridis or a list of colors.. Default: Reds
#' @param bins Number of bins to use in dividing expression levels.. Default: 10
#' @param plot_height Plot height in pixels. Default: 900
#' @param plot_width Plot width in pixels. Default: 900
#' @param plot_title  Display title with the name of the feature?. Default TRUE
#' @param pt_info Meta.data columns to add to the hoverinfo popup.. Default: ident
#' @param legend Display legend?. Default: TRUE
#' @param legend_font_size Legend font size. Default: 12
#' @param return Return the plot dataframe instead of displaying it. Default: FALSE
#' @param assay_use
#' @param slot_use
#' @param reverse_color_scale
#'
#' @importFrom plotly plot_ly layout
#'
#' @return plotly object
#' @export
#'
#' @examples
FeaturePlotly <- function(object,
                          feature = NULL,
                          assay_use = "RNA",
                          slot_use = "data",
                          return = FALSE,
                          pt_scale = 5,
                          pt_shape = "circle",
                          opacity = 1,
                          reduction = "umap",
                          dim_1 = 1,
                          dim_2 = 2,
                          colors_use = "Red",
                          reverse_color_scale = FALSE,
                          bins = 10,
                          plot_height = 900,
                          plot_width = 900,
                          plot_title = FALSE,
                          pt_info = NULL,
                          legend = TRUE,
                          legend_font_size = 12){

  df <- PrepDr(object,
               reduction,
               dim_1 = dim_1,
               dim_2 = dim_2)

  df <- PrepInfo(object = object,
                 pt_info = pt_info,
                 df = df)

  if (is.null(feature)){
    stop("No gene or feature given")
    }

  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature,
                         bins = bins,
                         use.scaled = TRUE,
                         assay_use = assay_use)
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature,
                         use.scaled = FALSE,
                         assay_use = assay_use,
                         suffix = "size")
  df[,ncol(df)] <- df[,ncol(df)] * pt_scale

  md <- GetFeatureValues(object = object,
                         features = c(pt_info, "ident")) %>%
    mutate_at(vars(-cell),
              list(~paste0('</br> ', substitute(.), ": ", .))) %>%
    unite(info, -cell)

  df %<>% inner_join(md)

  if (isTRUE(plot_title)){
    plot_title = feature
  } else {
    plot_title = NULL
  }

  pal <- PrepQualitativePalette(bins = binds,
                                palette_use = colors_use)

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               mode = 'markers',
               type = "scattergl",
               size = ~get(str_glue("{feature}_size")),
               sizes = c(0,max(df[[str_glue("{feature}_size")]])),
               marker = list(symbol = pt_shape,
                             opacity = opacity,
                             color = ~feature,
                             line = list(width = 0),
                             colorscale=pal,
                             reversescale = reverse_color_scale,
                             cmin = 0,
                             cmax = 1,
                             sizemode = "diameter",
                             colorbar = list(title = feature)),
               width = plot_width,
               height = plot_height,
               showlegend = legend,
               hoverinfo = "text",
               text = ~meta.info) %>%
    layout(
      title = plot_title,
      xaxis = list(title = str_glue("{reduction}_{dim_1}")),
      yaxis = list(title = str_glue("{reduction}_{dim_2}"))
    )

  p <- p %>% layout(legend = list(
    font = list(
      size = legend_font_size)
  )
  )

  if (isTRUE(return)){
    return(df)
  } else {
    return(p)
  }
}
