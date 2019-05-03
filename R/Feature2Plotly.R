#' Feature2Plotly
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of two different features.
#'
#' @param object Seurat object
#' @param feature_1 First variable to display. Currently only works with gene names
#' @param feature_2 Second variable to display. Currently only works with gene names
#' @param reduction Dimensional reduction to display (default: tsne)
#' @param dim_1 Dimension to display on the x-axis (default: 1)
#' @param dim_2 Dimension to display on the y-axis (default: 2)
#' @param pt_scale Factor by which to multiply the size of the points (default: 5)
#' @param pt_shape Shape to use for the points (default = circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors_1 Color palette to use for feature 1.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param colors_2 Color palette to use for feature 2.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param bins Number of bins to use in dividing expression levels. (default: 10)
#' @param plot_height Plot height in pixels (default: 900)
#' @param plot_width Plot width in pixels (default: 900)
#' @param plot_title  Display title with the name of the feature? (default TRUE)
#' @param pt_info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param legend Display legend? (default: TRUE)
#' @param legend_font_size Legend font size (default: 12)
#' @param return Return the plot dataframe instead of displaying it (default: FALSE)
#'
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly layout
#' @importFrom grDevices colorRampPalette
#' @importFrom stringr str_glue
#'
#' @return
#' @export
#'
#' @examples
Feature2Plotly <- function(object,
                           feature_1 = NULL,
                           feature_2 = NULL,
                           assay_1 = "RNA",
                           assay_2 = "RNA",
                           return = FALSE,
                           pt_scale = 1.5,
                           pt_shape = "circle",
                           opacity = 0.5,
                           reduction = "tsne",
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

  df <- PrepDf(object = object,
               reduction = reduction,
               dim_1 = dim_1,
               dim_2 = dim_2)

  #if (is.null(feature)){ stop("No gene or feature given") }
  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (colors_1 %in% rownames(brewer.pal.info)){
    pal.1 <- colorRampPalette(brewer.pal(brewer.pal.info[colors_1,]$maxcolors,colors_1))(bins)
  } else if (colors_1 %in% viridis_palettes){
    pal.1 <- viridis(n = bins, option = colors_1)
  } else {
    pal.1 <- colors_1
  }

  if (colors_2 %in% rownames(brewer.pal.info)){
    pal.2 <- colorRampPalette(brewer.pal(brewer.pal.info[colors_2,]$maxcolors,colors_2))(bins)
  } else if (colors_2 %in% viridis_palettes){
    pal.2 <- viridis(n = bins, option = colors_2)
  } else {
    pal.2 <- colors_2
  }

  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_1,
                         bins = bins,
                         use.scaled = TRUE,
                         assay = assay_1)
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_1,
                         bins = bins,
                         use.scaled = FALSE,
                         assay = assay_1,
                         suffix = "size")
  df[,ncol(df)] <- df[,ncol(df)] * pt_scale

  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_2,
                         bins = bins,
                         use.scaled = TRUE,
                         assay = assay_2)
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_2,
                         bins = bins,
                         use.scaled = FALSE,
                         assay = assay_2,
                         suffix = "size")
  df[,ncol(df)] <- df[,ncol(df)] * pt_scale

  if(!is.null(pt_info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt_info
      rowinfo = ""
      for(j in 1:length(pt_info)){
        rowinfo <- str_glue("{rowinfo} </br> {pt_info[j]}: {object@meta.data[i, pt_info[j]]}")
      }
      rowinfo <- str_glue("{rowinfo} </br> Expr {feature_1}: {df[i,'feature_1']}")
      rowinfo <- str_glue("{rowinfo} </br> Expr {feature_2}: {df[i,'feature_2']}")
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- meta.info
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               color = ~feature_1,
               mode = 'markers',
               colors = c(pal.1,pal.2),
               size = ~get(tr_glue("{feature_1}_size")),
               sizes = c(0,max(df[[str_glue("{feature_1}_size")]])),
               marker = list(symbol = pt_shape,
                             opacity = opacity,
                             sizemode = "diameter"
               ),
               width = plot_width,
               height = plot_height,
               type = 'scattergl',
               legendgroup = 1,
               showlegend = legend,
               name = feature_1) %>%
    add_trace(df,
              x = ~x,
              y = ~y,
              color = ~feature.2,
              mode = 'markers',
              type = 'scattergl',
              size = ~get(str_glue("{feature_1}_size")),
              sizes = c(0,max(df[[str_glue("{feature_1}_size")]])),
              marker = list(symbol = pt_shape,
                            opacity = opacity,
                            sizemode = "diameter"),
              legendgroup = 1,
              showlegend = legend,
              name = feature_2) %>%
    layout(title = str_glue("{feature_1} x {feature_2}"),
           xaxis = list(title = dim.axes[as.numeric(dim_1)]),
           yaxis = list(title = dim.axes[as.numeric(dim_2)]),
           margin = c(100,NA,NA,NA)
    )

  if(!is.null(pt_info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = ~meta.info,
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
