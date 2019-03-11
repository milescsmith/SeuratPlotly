#' DimPlotly3D
#'
#' @param object Seurat object
#' @param grouping Variable by which to group cells. Currently only works with the current ident and column names from meta.data (default: ident)
#' @param reduction_use Dimensional reduction to display (default: tsne)
#' @param dim_1 Dimension to display on the x-axis (default: 1)
#' @param dim_2 Dimension to display on the y-axis (default: 2)
#' @param dim_3 Dimension to display on the z-axis (default: 3)
#' @param do.label Add a label showing thr group name to the graph (default: FALSE)
#' @param label.size Label font size (default: 12)
#' @param show.arrow Offset the position of the labels and instead point to each group with an arrow (default: FALSE)
#' @param label.color Color for label border and arrow.  Need hex value. (default = '000000')
#' @param return Return the plot object instead of displaying it (default: FALSE)
#' @param pt.size Size of the points in pixels (default: 2)
#' @param pt_shape Shape to use for the points (default: circle)
#' @param opacity Transparency level to use for the points on a 0-1 scale (default: 1)
#' @param palette.use Color palette to use.  Must be a palette available in the Paletteer package.  (default: 'Set1')
#' @param plot_height Plot height in pixels (default: 900)
#' @param plot_width Plot width in pixels (default: 900)
#' @param pt_info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param legend Display legend? (default: TRUE)
#' @param legend_font_size Legend font size (default: 12)
#' @param plot_title Plot title (default: reduction_use)
#' @param plot_axes Display the major x, y, and z axes? (default: FALSE)
#' @param plot_grid Display the major unit tick marks? (default: FALSE)
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#'
#' @return plotly object
#' @export
#'
#' @examples
#' DimPlotly3D(object, grouping = "res.0.6", do.label = TRUE, show.arrow = FALSE)
DimPlotly3d <- DimPlotly3D <- function(object,
                                       grouping = "ident",
                                       do.label = FALSE,
                                       label.size = 12,
                                       label.color = "000000",
                                       show.arrow = FALSE,
                                       return = FALSE,
                                       pt.size = 2,
                                       pt_shape = "circle",
                                       opacity = 1,
                                       reduction_use = "tsne",
                                       dim_1 = 1,
                                       dim_2 = 2,
                                       dim_3 = 3,
                                       palette.use = "Set1",
                                       plot_height = 900,
                                       plot_width = 900,
                                       plot_title = NULL,
                                       pt_info = NULL,
                                       legend = TRUE,
                                       legend_font_size = 12,
                                       plot_grid = FALSE,
                                       plot_axes = FALSE) {
  df <- PrepDf(object,
    reduction_use,
    dim_1 = dim_1,
    dim_2 = dim_2,
    dim_3 = dim_3,
    grouping = grouping
  )

  df <- PrepInfo(
    object = object,
    pt_info = pt_info,
    df = df
  )

  if (do.label) {
    df %>%
      dplyr::group_by(ident) %>%
      centers() <- summarize(
      x = median(x = x),
      y = median(x = y),
      z = median(x = z)
    )
    labels <- list(
      x = centers$x,
      y = centers$y,
      z = centers$z,
      text = centers$ident,
      font = list(size = label.size)
    )
    compiled.labels <- list()

    if (isTRUE(show.arrow)) {
      border.color <- label.color
      bg.color <- "FFFFFF"
    } else {
      border.color <- "FFFFFF"
      bg.color <- "FFFFFF"
    }
    for (k in 1:length(unique(ident))) {
      tmp <- list(
        showarrow = show.arrow,
        x = labels$x[k],
        y = labels$y[k],
        z = labels$z[k],
        text = labels$text[k],
        font = list(size = label.size),
        bordercolor = border.color,
        bgcolor = bg.color,
        opacity = 0.8
      )
      compiled.labels <- c(compiled.labels, list(tmp))
    }
  } else {
    compiled.labels <- NULL
  }

  if (is.null(plot_title)) {
    plot_title <- reduction_use
  }

  pal <- PrepPalette(
    df = df,
    palette.use = palette.use
  )

  p <- plot_ly(df,
    x = ~x,
    y = ~y,
    z = ~z,
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
    type = "scatter3d",
    mode = "markers",
    showlegend = legend
  ) %>%
    layout(
      title = plot_title,
      scene = list(
        aspectratio = list(
          x = 0,
          y = 0,
          z = -1
        ),
        camera = list(
          center = list(
            x = 0,
            y = 0,
            z = 0
          ),
          eye = list(
            x = 2,
            y = -1,
            z = 0.5
          ),
          up = list(
            x = 1,
            y = 0,
            z = 0
          )
        ),
        dragmode = "turnable",
        xaxis = list(
          title = glue("{reduction_use}_{dim_1}"),
          type = "double",
          showgrid = plot_grid,
          visible = plot_axes
        ),
        yaxis = list(
          title = glue("{reduction_use}_{dim_2}"),
          type = "double",
          showgrid = plot_grid,
          visible = plot_axes
        ),
        zaxis = list(
          title = glue("{reduction_use}_{dim_3}"),
          type = "double",
          showgrid = plot_grid,
          visible = plot_axes
        ),
        annotations = compiled.labels
      )
    )

  if (!is.null(pt_info)) {
    p <- p %>% add_markers(
      hoverinfo = "text",
      hovertext = ~meta.info,
      showlegend = FALSE
    )
  }

  p <- p %>% layout(legend = list(
    font = list(
      size = legend_font_size
    )
  ))

  if (isTRUE(return)) {
    return(p)
  } else {
    p
  }
}
