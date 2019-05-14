#' FeaturePlotly
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' object@@dr slot
#'
#' @param object Seurat object
#' @param feature Feature values to display. Works with anything \code{Seurat::\link[Seurat]{FetchData}} can retrieve.
#' @param reduction Dimensional reduction to display. Default: umap
#' @param assay Assay to pull values from for feature. Default: NULL
#' @param slot Slot to pull values from for feature. Default: NULL
#' @param filter_zeros Remove points with no expression data. Default: TRUE
#' @param dim_1 Dimension to display on the x-axis. Default: 1
#' @param dim_2 Dimension to display on the y-axis. Default: 2
#' @param pt_scale Factor by which to multiply the size of the points. Default: 5
#' @param pt_shape Shape to use for the points. Default = circle
#' @param opacity Transparency level to use for the points, on a 0-1 scale. Default: 1
#' @param colors_use Color palette to use. Default: Reds
#' @param bins Number of bins to use in dividing expression levels.. Default: 10
#' @param plot_height Plot height in pixels. Default: 900
#' @param plot_width Plot width in pixels. Default: 900
#' @param plot_title  Display title with the name of the feature?. Default TRUE
#' @param pt_info Meta.data columns to add to the hoverinfo popup. Default: ident
#' @param legend Display legend?. Default: TRUE
#' @param legend_font_size Legend font size. Default: 12
#' @param return Return the plot dataframe instead of displaying it. Default: FALSE
#' @param reverse_color_scale Reverse the color scale. Default = FALSE
#'
#' @importFrom plotly plot_ly layout
#' @importFrom dplyr mutate_at vars inner_join
#' @importFrom tidyr unite

#' @return plotly object
#' @export
#'
FeaturePlotly <- function(object,
                          feature = NULL,
                          reduction = "umap",
                          assay = NULL,
                          slot = NULL,
                          filter_zeros = TRUE,
                          dim_1 = 1,
                          dim_2 = 2,
                          pt_scale = 5,
                          pt_shape = "circle",
                          opacity = 1,
                          colors_use = "Red",
                          reverse_color_scale = FALSE,
                          bins = 10,
                          plot_height = 900,
                          plot_width = 900,
                          plot_title = FALSE,
                          pt_info = NULL,
                          legend = TRUE,
                          legend_font_size = 12,
                          return = FALSE){

  # global variable binding hack
  cell <- NULL
  info <- NULL

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
                         assay = assay,
                         suffix = "expr")
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature,
                         use.scaled = TRUE,
                         assay = assay,
                         suffix = "size")
  df[[str_glue("{feature}_size")]] <- df[[str_glue("{feature}_size")]] * pt_scale

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

  pal <- PrepQuantitativePalette(bins = bins,
                                palette = colors_use)

  if (filter_zeros){
    df %<>% filter(get(str_glue("{feature}_size")) > 0)
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               mode = "markers",
               type = "scattergl",
               marker = list(symbol = pt_shape,
                             opacity = opacity,
                             size = ~get(str_glue("{feature}_size")),
                             sizes = c(0,max(df[[str_glue("{feature}_size")]])),
                             color = ~get(str_glue("{feature}_expr")),
                             colors=pal,
                             line = list(width = 0)),
               width = plot_width,
               height = plot_height,
               showlegend = legend,
               hoverinfo = "text",
               text = ~meta_info) %>%
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
