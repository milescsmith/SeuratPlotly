#' Feature2Plotly
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of two different features.
#'
#' @param object Seurat object
#' @param reduction Dimensional reduction to display. Default: 'umap'
#' @param feature_1 First feature values to display. Works with anything \code{Seurat::\link[Seurat]{FetchData}} can retrieve.
#' @param feature_2 Second feature values to display.
#' @param assay_1 Assay to pull values from for feature_1. Default: NULL
#' @param assay_2 Assay to pull values from for feature_2. Default: NULL
#' @param colors_1 Colors to use to display values for feature_1. Default: "Reds"
#' @param colors_2 Colors to use to display values for feature_2. Default: "Blues"
#' @param dim_1 Dimension to display on the x-axis. Default: 1
#' @param dim_2 Dimension to display on the y-axis. Default: 2
#' @param pt_scale Factor by which to multiply the size of the points. Default: 5
#' @param pt_shape Shape to use for the points. Default = circle
#' @param opacity Transparency level to use for the points, on a 0-1 scale. Default: 1
#' @param bins Number of bins to use in dividing expression levels.. Default: 10
#' @param plot_height Plot height in pixels. Default: 900
#' @param plot_width Plot width in pixels. Default: 900
#' @param pt_info Meta.data columns to add to the hoverinfo popup.. Default: ident
#' @param legend Display legend?. Default: TRUE
#' @param legend_font_size Legend font size. Default: 12
#' @param return Return the plot dataframe instead of displaying it. Default: FALSE
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom plotly plot_ly layout
#' @importFrom grDevices colorRampPalette
#' @importFrom stringr str_glue
#'
#' @return
#' @export
#'
Feature2Plotly <- function(object,
                           feature_1 = NULL,
                           feature_2 = NULL,
                           reduction = "umap",
                           assay_1 = NULL,
                           assay_2 = NULL,
                           return = FALSE,
                           pt_scale = 1.5,
                           pt_shape = "circle",
                           opacity = 0.5,
                           dim_1 = 1,
                           dim_2 = 2,
                           colors_1 = "Reds",
                           colors_2 = "Blues",
                           bins = 10,
                           plot_height = "750",
                           plot_width = "750",
                           pt_info = NULL,
                           legend = TRUE,
                           legend_font_size = 12){

  cell <- NULL
  info <- NULL

  df <- PrepDr(object = object,
               reduction = reduction,
               dim_1 = dim_1,
               dim_2 = dim_2)

  pal_1 <- PrepQuantitativePalette(bins, colors_1)
  pal_2 <- PrepQuantitativePalette(bins, colors_2)

  #get feature_1 values
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_1,
                         bins = bins,
                         use.scaled = TRUE,
                         assay = assay_1,
                         suffix = "expr")
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_1,
                         use.scaled = TRUE,
                         assay = assay_1,
                         suffix = "size")
  df[[str_glue("{feature_1}_size")]] <- df[[str_glue("{feature_1}_size")]] * pt_scale

  #get feature_2 values
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_2,
                         bins = bins,
                         use.scaled = TRUE,
                         assay = assay_2,
                         suffix = "expr")
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_2,
                         use.scaled = TRUE,
                         assay = assay_2,
                         suffix = "size")
  df[[str_glue("{feature_2}_size")]] <- df[[str_glue("{feature_2}_size")]] * pt_scale


  md <- GetFeatureValues(object = object,
                         features = c(pt_info, "ident")) %>%
    mutate_at(vars(-cell),
              list(~paste0('</br> ', substitute(.), ": ", .))) %>%
    unite(info, -cell)

  df %<>% inner_join(md)

 p <- plot_ly(
  filter(
    df,
    get(str_glue("{feature_1}_expr")) > 0
  ),
  x = ~x,
  y = ~y,
  mode = "markers",
  type = "scattergl",
  color = ~ get(str_glue("{feature_1}_expr")),
  colors = pal_1,
  marker = list(
    symbol = pt_shape,
    opacity = opacity,
    size = ~ get(str_glue("{feature_1}_size")),
    sizes = c(
      0,
      max(df[[str_glue("{feature_1}_size")]])
    )
  ),
  width = plot_width,
  height = plot_height,
  showlegend = legend,
  name = feature_1
  ) %>%
  colorbar(
    title = feature_1,
    which = 1
  ) %>%
  add_trace(
    data = filter(
      df,
      get(str_glue("{feature_2}_expr")) > 0
    ),
    x = ~x,
    y = ~y,
    color = ~ get(str_glue("{feature_2}_expr")),
    colors = pal_2,
    marker = list(
      symbol = pt_shape,
      opacity = opacity,
      size = ~ get(str_glue("{feature_2}_size")),
      sizes = c(0, max(df[[str_glue("{feature_2}_size")]]))
    ),
    name = feature_2
  ) %>%
  colorbar(
    title = feature_2,
    which = 2
  ) %>%
  layout(
    title = str_glue("{feature_1} x {feature_2}"),
    xaxis = list(title = str_glue("{reduction}_{dim_1}")),
    yaxis = list(title = str_glue("{reduction}_{dim_2}")),
    margin = c(100, NA, NA, NA)
  )

  if(!is.null(pt_info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = ~meta_info,
                           showlegend = FALSE
    )
  }

  p <- p %>% layout(legend = list(
    font = list(
      size = legend_font_size)))

  if (isTRUE(return)){
    return(df)
  } else {
    p
  }
}
