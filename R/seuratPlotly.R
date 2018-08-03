#seuratPlotly.R: Functions replacing the ggplot2-based ploting functions of Seurat with those of Plot.ly

#' Plot dimensional reduction for a Seurat object
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring points by the given grouping variable
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' seuratObj@@dr slot
#'
#' @param seuratObj Seurat object
#' @param group.by Variable by which to group cells. Currently only works with the current ident and column names from meta.data (default: ident)
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param do.label Add a label showing thr group name to the graph (default: FALSE)
#' @param label.size Label font size (default: 12)
#' @param show.arrow Offset the position of the labels and instead point to each group with an arrow (default: FALSE)
#' @param pt.size Size of the points in pixels (default: 2)
#' @param pt.shape Shape to use for the points (default: circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors.use Color palette to use.  Palettes from RColorBrewer and viridis
#' @param plot.height Plot height in pixels (default: 900)
#' @param plot.width Plot width in pixels (default: 900)
#' @param legend Display legend? (default: TRUE)
#' @param legend.font.size Legend font size (default: 12)
#' @param pt.info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param do.return Return the plot object instead of displaying it (default: FALSE)
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom Seurat GetDimReduction
#' @importFrom Seurat FetchData
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom grDevices colorRampPalette
#'
#' @return plotly object
#' @export
#'
#' @examples
#' seuratObj <- RunTSNE(seuratObj)
#' DimPlotly(seuratObj, group.by = 'ident', pt.size = 4, opacity = 0.5, plot.title = "Test Plot", reduction.use = "tsne")
DimPlotly <- function(seuratObj,
                      group.by = "ident",
                      do.label = FALSE,
                      label.size = 12,
                      show.arrow = FALSE,
                      do.return = FALSE,
                      pt.size = 4,
                      pt.shape = "circle",
                      opacity = 0.75,
                      reduction.use = "tsne",
                      dim.1 = 1,
                      dim.2 = 2,
                      colors.use = "Set1",
                      plot.height = 900,
                      plot.width = 900,
                      pt.info = NULL,
                      legend = TRUE,
                      legend.font.size = 12){

  df <- as.data.frame(GetDimReduction(object = seuratObj,
                                      reduction.type = reduction.use,
                                      slot = "cell.embeddings"))
  dim.code <- GetDimReduction(
    object = seuratObj,
    reduction.type = reduction.use,
    slot = "key"
  )

  dim.axes <- colnames(
    GetDimReduction(
      object = seuratObj,
      reduction.type = reduction.use,
      slot = "cell.embeddings"
    )
  )

  dim.code <- c(dim.axes[[dim.1]], dim.axes[[dim.2]])
  df <- df[,dim.code]
  cell_names <- rownames(df)
  ident <- as.factor(x = seuratObj@ident)
  if (group.by != "ident") {
    ident <- as.factor(x = FetchData(
      object = seuratObj,
      vars.all = group.by
    )[, 1])
  }

  df$x <- df[,1]
  df$y <- df[,2]

  df$ident <- ident

  # Using the members of pt.info,
  # build a string for each row that compiles the indicated
  # pt.info from the meta.data slot along with the necessary
  # html tags and column names
  if(!is.null(pt.info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt.info
      rowinfo = ""
      for(j in 1:length(pt.info)){
        rowinfo <- paste0(rowinfo, " </br> ", pt.info[j], ": ", seuratObj@meta.data[i, pt.info[j]])
      }
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- meta.info
  } else {
    df$meta.info <- df$ident
  }

  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")
  bins = length(unique(df[,'ident']))

  if (colors.use %in% rownames(brewer.pal.info)){
    pal <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use,]$maxcolors,colors.use))(bins)
  } else if (colors.use %in% viridis_palettes){
    pal <- viridis(n = bins, option = colors.use)
  } else {
    pal <- colors.use
  }

  if (do.label) {
    df %>%
      dplyr::group_by(ident) %>%
      summarize(x = median(x = x), y = median(x = as.double(y))) -> centers
    labels <- list(x = centers$x,
                   y = centers$y,
                   text = centers$ident,
                   font = list(size = label.size),
                   showarrow = show.arrow,
                   bordercolor='000000',
                   bgcolor='FFFFFF',
                   opacity=1
    )
  } else {
    labels <- NULL
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               color = ~ident,
               colors = pal,
               marker = list(symbol = pt.shape,
                             size = pt.size,
                             opacity = opacity,
                             mode = "markers",
                             type = "scattergl"),
               width = plot.width,
               height = plot.height,
               type = "scattergl",
               mode = "markers",
               showlegend = legend,
               hoverinfo = "text",
               text = ~meta.info) %>%
    layout(title = reduction.use,
           xaxis = list(title = dim.axes[as.numeric(dim.1)]),
           yaxis = list(title = dim.axes[as.numeric(dim.2)]))

  if(isTRUE(do.label)){
    p <- p %>% layout(annotations = labels,
                      legend = list(
                        font = list(
                          size = legend.font.size)))
  } else {
    p <- p %>% layout(legend = list(
                        font = list(
                          size = legend.font.size)))
  }


  if (isTRUE(do.return)){
    return(p)
  } else {
    p
  }
}


#' DimPlotly3D
#'
#' @param seuratObj Seurat object
#' @param group.by Variable by which to group cells. Currently only works with the current ident and column names from meta.data (default: ident)
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param dim.3 Dimension to display on the z-axis (default: 3)
#' @param do.label Add a label showing thr group name to the graph (default: FALSE)
#' @param label.size Label font size (default: 12)
#' @param show.arrow Offset the position of the labels and instead point to each group with an arrow (default: FALSE)
#' @param do.return Return the plot object instead of displaying it (default: FALSE)
#' @param pt.size Size of the points in pixels (default: 2)
#' @param pt.shape Shape to use for the points (default: circle)
#' @param opacity Transparency level to use for the points on a 0-1 scale (default: 1)
#' @param colors.use Color palette to use.  Accepts palettes from RColorBrewer and viridis (default: Set1)
#' @param plot.height Plot height in pixels (default: 900)
#' @param plot.width Plot width in pixels (default: 900)
#' @param pt.info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param legend Display legend? (default: TRUE)
#' @param legend.font.size Legend font size (default: 12)
#' @param plot.title Plot title (default: reduction.use)
#' @param plot.axes Display the major x, y, and z axes? (default: FALSE)
#' @param plot.grid Display the major unit tick marks? (default: FALSE)
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom Seurat GetDimReduction
#' @importFrom Seurat FetchData
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom grDevices colorRampPalette
#'
#' @return plotly object
#' @export
#'
#' @examples DimPlotly3D(seuratObj, group.by = "res.0.6", do.label = TRUE, show.arrow = FALSE)
#'
DimPlotly3D <- function(seuratObj,
                        group.by = "ident",
                        do.label = FALSE,
                        label.size = 12,
                        do.return = FALSE,
                        pt.size = 2,
                        pt.shape = "circle",
                        opacity = 1,
                        reduction.use = "tsne",
                        dim.1 = 1,
                        dim.2 = 2,
                        dim.3 = 3,
                        colors.use = "Set1",
                        plot.height = 900,
                        plot.width = 900,
                        plot.title = NULL,
                        pt.info = NULL,
                        show.arrow = FALSE,
                        legend = TRUE,
                        legend.font.size = 12,
                        plot.grid = FALSE,
                        plot.axes = FALSE){

  df <- as.data.frame(GetDimReduction(object = seuratObj,
                                      reduction.type = reduction.use,
                                      slot = "cell.embeddings"))
  if (length(colnames(df)) < 3){
    stop(paste("Fewer than 3 dimensions have been calculated for ", reduction.use, "."))
  }
  dim.code <- GetDimReduction(
    object = seuratObj,
    reduction.type = reduction.use,
    slot = "key"
  )

  dim.axes <- colnames(
    GetDimReduction(
      object = seuratObj,
      reduction.type = reduction.use,
      slot = "cell.embeddings"
    )
  )

  dim.code <- c(dim.axes[[dim.1]], dim.axes[[dim.2]], dim.axes[[dim.3]])
  df <- df[,dim.code]
  cell_names <- rownames(df)
  ident <- as.factor(x = seuratObj@ident)
  if (group.by != "ident") {
    ident <- as.factor(x = FetchData(
      object = seuratObj,
      vars.all = group.by
    )[, 1])
  }

  df$x <- df[,1]
  df$y <- df[,2]
  df$z <- df[,3]
  df$ident <- ident
  rownames(df) <- cell_names

  if(!is.null(pt.info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt.info
      rowinfo = ""
      for(j in 1:length(pt.info)){
        rowinfo <- paste0(rowinfo, " </br> ", pt.info[j], ": ", seuratObj@meta.data[i, pt.info[j]])
      }
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- meta.info
  } else {
    df$meta.info <- df$ident
  }

  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (do.label) {
    df %>%
      dplyr::group_by(ident) %>%
      summarize(x = median(x = x), y = median(x = y), z = median(x = z)) -> centers
    labels <- list(x = centers$x,
                   y = centers$y,
                   z = centers$z,
                   text = centers$ident,
                   font = list(size = label.size)
    )
    compiled.labels = list()

    if(isTRUE(show.arrow)){
      border.color = '000000'
      bg.color = 'FFFFFF'
    } else {
      border.color = 'FFFFFF'
      bg.color = 'FFFFFF'
    }
    for(k in 1:length(unique(ident))){
      tmp <- list(
        showarrow = show.arrow,
        x = labels$x[k],
        y = labels$y[k],
        z = labels$z[k],
        text = labels$text[k],
        font = list(size = label.size),
        bordercolor=border.color,
        bgcolor=bg.color,
        opacity=0.8)
      compiled.labels <- c(compiled.labels, list(tmp))
    }
  } else {
    compiled.labels <- NULL
  }

  if (is.null(plot.title)){
    plot.title <- reduction.use
  }

  bins = length(unique(df[,'ident']))
  if (as.character(colors.use) %in% rownames(brewer.pal.info)){
    pal <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use,]$maxcolors,colors.use))(bins)
  } else if (colors.use %in% viridis_palettes){
    pal <- viridis(n = bins, option = colors.use)
  } else {
    pal <- colors.use
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               z = ~z,
               color = ~ident,
               colors = pal,
               marker = list(symbol = pt.shape,
                             "size" = pt.size,
                             opacity = opacity,
                             mode = "markers"),
               width = plot.width,
               height = plot.height,
               type = "scatter3d",
               mode = "markers",
               showlegend = legend
  ) %>%
    layout(
      title = plot.title,
      scene = list(
        aspectratio = list(x = 0,y = 0,z = -1),
        camera = list(
          center = list(x = 0,y = 0,z = 0),
          eye = list(x = 2,y = -1,z = 0.5),#eye = list(x = 0,y = 0,z = max(mean(df$z), mean(df$y))),
          up = list(x = 1,y = 0,z = 0)
        ),
        dragmode = "turnable",
        xaxis = list(title = dim.axes[as.numeric(dim.1)], type = "double", showgrid = plot.grid, visible = plot.axes),
        yaxis = list(title = dim.axes[as.numeric(dim.2)], type = "double", showgrid = plot.grid, visible = plot.axes),
        zaxis = list(title = dim.axes[as.numeric(dim.3)], type = "double", showgrid = plot.grid, visible = plot.axes),
        annotations = compiled.labels
      )
    )

  if(!is.null(pt.info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = ~meta.info,
                           showlegend = FALSE
    )
  }

  p <- p %>% layout(legend = list(
                      font = list(
                        size = legend.font.size)))

  if (isTRUE(do.return)){
    return(p)
  } else {
    p
  }
}

#' FeaturePlotly
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' seuratObj@@dr slot
#'
#' @param seuratObj Seurat object
#' @param feature.use Variable to display. Currently only works with gene names
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param pt.scale Factor by which to multiply the size of the points (default: 5)
#' @param pt.shape Shape to use for the points (default = circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors.use Color palette to use.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param bins Number of bins to use in dividing expression levels. (default: 10)
#' @param plot.height Plot height in pixels (default: 900)
#' @param plot.width Plot width in pixels (default: 900)
#' @param plot.title  Display title with the name of the feature? (default TRUE)
#' @param pt.info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param legend Display legend? (default: TRUE)
#' @param legend.font.size Legend font size (default: 12)
#' @param do.return Return the plot object instead of displaying it (default: FALSE)
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom Seurat GetDimReduction
#' @importFrom Seurat FetchData
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom grDevices colorRampPalette
#'
#' @return plotly object
#' @export
#'
#' @examples
FeaturePlotly <- function(seuratObj,
                          feature.use = NULL,
                          do.return = TRUE,
                          pt.scale = 5,
                          pt.shape = "circle",
                          opacity = 1,
                          reduction.use = "tsne",
                          dim.1 = 1,
                          dim.2 = 2,
                          colors.use = "Reds",
                          bins = 10,
                          plot.height = 900,
                          plot.width = 900,
                          plot.title = FALSE,
                          pt.info = NULL,
                          legend = TRUE,
                          legend.font.size = 12){

  df <- as.data.frame(GetDimReduction(object = seuratObj,
                                      reduction.type = reduction.use,
                                      slot = "cell.embeddings"))
  dim.code <- GetDimReduction(
    object = seuratObj,
    reduction.type = reduction.use,
    slot = "key"
  )

  dim.axes <- colnames(
    GetDimReduction(
      object = seuratObj,
      reduction.type = reduction.use,
      slot = "cell.embeddings"
    )
  )

  dim.code <- c(dim.axes[[dim.1]], dim.axes[[dim.2]])
  df <- df[,dim.code]
  cell_names <- rownames(df)
  ident <- as.factor(x = seuratObj@ident)


  df$x <- df[,1]
  df$y <- df[,2]

  if(!is.null(pt.info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt.info
      rowinfo = ""
      for(j in 1:length(pt.info)){
        rowinfo <- paste0(rowinfo, " </br> ", pt.info[j], ": ", seuratObj@meta.data[i, pt.info[j]])
      }
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- meta.info
  } else {
    df$meta.info <- seuratObj@ident
  }

  if (is.null(feature.use)){ stop("No gene or feature given") }

  feature.data <- FetchData(object = seuratObj,
                            vars.all = feature.use,
                            use.scaled = TRUE)
  size.data <- FetchData(object = seuratObj,
                         vars.all = feature.use,
                         use.scaled = TRUE)
  feature.data[,1][feature.data[,1] == 0] <- NA
  feature.data <- as.matrix(feature.data)
  cut.feature.data <- as.numeric(as.factor(x = cut(x = as.numeric(feature.data), breaks = bins)))

  df[,"feature"] <- cut.feature.data
  df[,"size"] <- size.data[,1] * pt.scale

  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (colors.use %in% rownames(brewer.pal.info)){
    pal <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use,]$maxcolors,colors.use))(bins)
  } else if (colors.use %in% viridis_palettes){
    pal <- viridis(n = bins, option = colors.use)
  } else {
    pal <- colors.use
  }

  if (isTRUE(plot.title)){
    plot.title = feature.use
  } else {
    plot.title = NULL
  }

  p <- plot_ly(df,
               x = ~get(dim.code[1]),
               y = ~get(dim.code[2]),
               color = ~feature,
               colors = pal,
               mode = 'markers',
               type = "scattergl",
               size = ~size,
               sizes = c(0,max(df$size)),
               marker = list(symbol = pt.shape,
                             opacity = opacity,
                             sizemode = "diameter"),
               width = plot.width,
               height = plot.height,
               showlegend = legend,
               hoverinfo = "text",
               text = ~meta.info) %>%
    layout(title = plot.title,
           xaxis = list(title = dim.axes[as.numeric(dim.1)]),
           yaxis = list(title = dim.axes[as.numeric(dim.2)]))

  if (isTRUE(do.return)){
    return(p)
  } else {
    p %>% layout(legend = list(
      font = list(
        size = legend.font.size)))
  }
}

#' FeaturePlotly3D
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' seuratObj@@dr slot
#'
#' @param seuratObj Seurat object
#' @param feature.use Variable to display. Currently only works with gene names
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param dim.3 Dimension to display on the z-axis (default: 3)
#' @param pt.scale Factor by which to multiply the size of the points (default: 5)
#' @param pt.shape Shape to use for the points (default = circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors.use Color palette to use.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param bins Number of bins to use in dividing expression levels. (default: 10)
#' @param plot.height Plot height in pixels (default: 900)
#' @param plot.width Plot width in pixels (default: 900)
#' @param plot.title  Display title with the name of the feature? (default TRUE)
#' @param plot.axes Display the major x, y, and z axes? (default: FALSE)
#' @param plot.grid Display the major unit tick marks? (default: FALSE)
#' @param pt.info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param legend Display legend? (default: TRUE)
#' @param legend.font.size Legend font size (default: 12)
#' @param do.return Return the plot object instead of displaying it (default: FALSE)
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom Seurat GetDimReduction
#' @importFrom Seurat FetchData
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom grDevices colorRampPalette
#'
#' @return If do.return is TRUE, a plotly object.
#' @export
#'
#' @examples
FeaturePlotly3D <- function(seuratObj,
                            feature.use = NULL,
                            do.return = FALSE,
                            pt.scale = 0.5,
                            pt.shape = "circle",
                            opacity = 1,
                            reduction.use = "tsne",
                            dim.1 = 1,
                            dim.2 = 2,
                            dim.3 = 3,
                            colors.use = "Reds",
                            bins = 10,
                            plot.height = 900,
                            plot.width = 900,
                            plot.title = FALSE,
                            pt.info = NULL,
                            legend = TRUE,
                            legend.font.size = 12,
                            plot.grid = FALSE,
                            plot.axes = FALSE){

  df <- as.data.frame(GetDimReduction(object = seuratObj,
                                      reduction.type = reduction.use,
                                      slot = "cell.embeddings"))
  dim.code <- GetDimReduction(
    object = seuratObj,
    reduction.type = reduction.use,
    slot = "key"
  )

  dim.axes <- colnames(
    GetDimReduction(
      object = seuratObj,
      reduction.type = reduction.use,
      slot = "cell.embeddings"
    )
  )

  dim.code <- c(dim.axes[[dim.1]], dim.axes[[dim.2]], dim.axes[[dim.3]])
  cell_names <- rownames(df)

  rownames(df) <- cell_names

  df$x <- df[,1]
  df$y <- df[,2]
  df$z <- df[,3]

  if(!is.null(pt.info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt.info
      rowinfo = ""
      for(j in 1:length(pt.info)){
        rowinfo <- paste0(rowinfo, " </br> ", pt.info[j], ": ", seuratObj@meta.data[i, pt.info[j]])
      }
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- seuratObj@ident
  }

  if (is.null(feature.use)){ stop("No gene or feature given") }

  feature.data <- FetchData(object = seuratObj,
                            vars.all = feature.use,
                            use.scaled = TRUE)
  size.data <- FetchData(object = seuratObj,
                         vars.all = feature.use,
                         use.scaled = FALSE)
  feature.data[,1][feature.data[,1] == 0] <- NA
  feature.data <- as.matrix(feature.data)
  cut.feature.data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature.data), breaks = bins)))
  df[,"feature"] <- cut.feature.data
  df[,"size"] <- size.data[,1] * pt.scale

  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (colors.use %in% rownames(brewer.pal.info)){
    pal <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use,]$maxcolors,colors.use))(bins)
  } else if (colors.use %in% viridis_palettes){
    pal <- viridis(n = bins, option = colors.use)
  } else {
    pal <- colors.use
  }

  if (isTRUE(plot.title)){
    plot.title = feature.use
  } else {
    plot.title = NULL
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               z = ~z,
               color = ~feature,
               mode = 'markers',
               type = "scatter3d",
               colors = pal,
               size = ~size,
               sizes = c(0,max(df$size)),
               marker = list(symbol = pt.shape,
                             opacity = opacity,
                             sizemode = "diameter"
               ),
               width = plot.width,
               height = plot.height,
               showlegend = legend) %>%
    layout(title = plot.title,
           scene = list(
             xaxis = list(title = dim.axes[as.numeric(dim.1)], showgrid = plot.grid, visible = plot.axes),
             yaxis = list(title = dim.axes[as.numeric(dim.2)], showgrid = plot.grid, visible = plot.axes),
             zaxis = list(title = dim.axes[as.numeric(dim.3)], showgrid = plot.grid, visible = plot.axes)
           ),
           margin = c(100,NA,NA,NA)
    )

  if(!is.null(pt.info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = paste(~meta.info, ~feature),
                           showlegend = FALSE,

    )
  }

  p <- p %>% layout(legend = list(
    font = list(
      size = legend.font.size)))

  if (isTRUE(do.return)){
    return(p)
  } else {
    p
  }
}

# TODO: add more ways to sort/slice/dice!
#' BubblePlotly
#'
#' Plot average expression levels and proportion expressing a feature arranged by given grouping.
#' Color intensity of the markers indicates expression levels and size of the marker indicates
#' the proportion of the group expressing the gene above a given threshold. By default, genes are arranged
#' along the x-axis and groups along the y-axis.
#'
#' @param seuratObj Seurat object
#' @param genes.plot A list of genes to plot for each group.
#' @param colors.use Color palette to use.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param dot.min Minimium marker size, in pixels. (default: 0)
#' @param dot.scale Factor by which to scale markers. (default: 2)
#' @param group.by Factor by which to group cells.  (default: ident)
#' @param legend Display legend. (currently nonfunctional) (default TRUE)
#' @param do.return Return the plot object instead of displaying it (default: FALSE)
#' @param x.lab.rot Angle by which to rotate the x-axis labels, in degrees relative to horizontally aligned text. (default -45°)
#' @param plot.height Plot height in pixels. (default: 900)
#' @param plot.width Plot width in pixels. (default: 900)
#' @param x.font.size Size of the x-axis titles. (default: 10)
#' @param y.font.size Size of the y-axis titles. (default: 10)
#' @param title.font.size Size of the plot title. (default: 12)
#' @param legend.text.size Size of the legend text. (default: 10)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param plot.title  Display title with the name of the feature? (default TRUE)
#' @param bins Number of bins to use in dividing expression levels. (default: 10)
#' @param flip Swap the x- and y-axes so that genes are on the y-axis and groups along the x-axis. (default: false)
#' @param alphabetize Alphabetize the display order of genes. (default: TRUE)
#' @param pct.expr.thresh Hide a gene if there is no group that expresses at or above this percentage. (default: 0)
#' @param export.df Return the generated data frame underlying the graph. (default: FALSE)
#' @param use.scaled Use scaled data. (default: FALSE)
#' @param use.raw Use raw data. (default: FALSE)
#' @param show.zeros Display comparison groups for which there there was no data? (default: FALSE)
#' @param add.group Add a null comparison group.  Will be initialized with zeros.
#' @param y.label.order List of labels to use for the comparison groups.  Will be displayed in order of their index.
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column
#' @importFrom Seurat GetDimReduction
#' @importFrom Seurat FetchData
#' @importFrom Seurat SetAllIdent
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom viridis viridis
#' @importFrom compositions normalize
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom grDevices colorRampPalette
#'
#' @return if do.return is TRUE, a plotly object.
#' @return if export.df is TRUE, a data frame.
#' @export
#'
#' @examples
BubblePlotly <- function (seuratObj,
                          genes.plot,
                          colors.use = "Blues",
                          dot.min = 0,
                          dot.scale = 2,
                          group.by,
                          plot.legend = FALSE,
                          do.return = FALSE,
                          x.lab.rot = -45,
                          plot.width = 600,
                          plot.height = 600,
                          x.font.size = 10,
                          y.font.size = 10,
                          title.font.size = 12,
                          legend.text.size = 10,
                          opacity = 1,
                          plot.title = NULL,
                          bins = 50,
                          flip = FALSE,
                          alphabetize = TRUE,
                          pct.expr.thresh = NULL,
                          export.df = FALSE,
                          use.scaled = FALSE,
                          use.raw = FALSE,
                          show.zeros = FALSE,
                          add.group = NULL,
                          y.label.order = NULL)
{

  if (!missing(x = group.by)) {
    seuratObj <- SetAllIdent(object = seuratObj, id = group.by)
  }

  #screen out any genes that are not in our dataset and print them
  original_genes_to_plot <- genes.plot
  genes.plot <- (genes.plot %>% as_tibble() %>% dplyr::filter(value %in% rownames(seuratObj@data)))$value
  not_found <- original_genes_to_plot[!original_genes_to_plot %in% genes.plot]
  print(not_found)

  #let the user alphabetize the list of genes.  makes it easier to find a particular gene in a big table
  #TODO: add more sort options
  if (isTRUE(alphabetize)){
    genes.plot <- sort(genes.plot, decreasing = TRUE)
  }

  if(isTRUE(use.scaled)){
    data.to.plot <- data.frame(FetchData(object = seuratObj,
                                                 vars.all = genes.plot,
                                                 use.scaled = TRUE))
  } else if(isTRUE(use.raw)){
    data.to.plot <- data.frame(FetchData(object = seuratObj,
                                                 vars.all = genes.plot,
                                                 use.scaled = FALSE))
  } else {
    data.to.plot <- data.frame(FetchData(object = seuratObj,
                                                 vars.all = genes.plot))
  }

  # Add 0 data for each gene that was not detected or failed to pass QC
  if (isTRUE(show.zeros)){
    for(i in not_found){
      data.to.plot[,i] <- 0
    }
  }

  data.to.plot <- rownames_to_column(df = data.to.plot, var = "cell")

  data.to.plot$id <- seuratObj@ident

  data.to.plot <- data.to.plot %>% gather(key = genes.plot,
                                          value = expression, -c(cell, id))

  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>%
    summarize(avg.exp = mean(expm1(x = expression)),
              pct.exp = PercentAbove(x = expression, threshold = 0),
              n = n())

  if(!is.null(pct.expr.thresh)){
    data.to.plot <- data.to.plot %>% group_by(genes.plot) %>% filter(max(pct.exp) > pct.expr.thresh)
  }

  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>%
    mutate(avg.exp.scale = normalize(x = avg.exp))


  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot,
                                    levels = rev(x = sub(pattern = "-",
                                                         replacement = ".",
                                                         x = genes.plot)))

  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA

  data.to.plot <- as.data.frame(data.to.plot)

  cut.avg.exp.scale.data <- factor(round(data.to.plot$avg.exp.scale*100,0)+1)

  data.to.plot[,"feature"] <- cut.avg.exp.scale.data
  data.to.plot[,"size"] <- data.to.plot$pct.exp

  if(!is.null(add.group)){
    for(i in add.group){
      for(j in genes.plot){
        data.to.plot <- data.to.plot %>%
          do(add_row(.,id = i,
                     genes.plot = j,
                     avg.exp = 0,
                     pct.exp = 0,
                     n = 0,
                     avg.exp.scale=0,
                     size=0,
                     feature=factor(1))) %>%
          distinct()
      }
    }
  }

  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (colors.use %in% rownames(brewer.pal.info)){
    pal <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use,]$maxcolors,colors.use))(100)
  } else if (colors.use %in% viridis_palettes){
    pal <- viridis(n = 100, option = colors.use)
  } else {
    pal <- colors.use
  }

  data.to.plot <- data.to.plot %>% arrange(id)

  genes.plot <- data.to.plot$genes.plot

  if(flip){
    ax = ~id
    ay = ~genes.plot
  } else {
    ax = ~genes.plot
    ay = ~id
  }

  p <- plot_ly(data = data.to.plot,
               x = ax,
               y = ay,
               hoverinfo = 'text',
               text = ~paste('Group:', id,
                             '<br>Gene:', genes.plot,
                             '<br>Avg normalized expression:', round(avg.exp.scale,2),
                             '<br>% expressing:', round(pct.exp*100,2)),
               type = 'scatter',
               mode = 'markers',
               color = ~feature,
               colors = pal,
               cmin = pal[1],
               cmax = pal[100],
               size = ~size,
               sizes = c(0,10)*dot.scale,
               marker = list(opacity = opacity,
                             symbol = 'circle',
                             sizemode = 'diameter'),
               showlegend = FALSE,
               width = plot.width,
               height = plot.height) %>%
    layout(title = plot.title,
           titlefont = list(size = title.font.size),
           #font = list(size = font.size),
           xaxis = list(title = "Gene", tickangle = x.lab.rot, tickfont = list(size = x.font.size)),
           yaxis = list(title = "Cluster", tickfont = list(size = y.font.size), type = "category", categoryorder = "array", categoryarray = y.label.order),
           margin = list(b = 100, l = 250)
    )

  if (do.return) {
    return(p)
  } else if(isTRUE(export.df))
  {
    return(data.to.plot)
  } else {
    p
  }
  # print(length(unique(data.to.plot$id)))
  # print(length(unique(data.to.plot$genes.plot)))
}

#' PercentAbove
#'
#' Return the percentage of a list that is above a threshold
#' @param x A list of numeric values.
#' @param threshold A numeric threshold value.
#'
#' @return double
#' @export
#'
#' @examples
PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

#' Feature2Plotly
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of two different features.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' seuratObj@@dr slot
#'
#' @param seuratObj Seurat object
#' @param feature.1.use First variable to display. Currently only works with gene names
#' @param feature.2.use Second variable to display. Currently only works with gene names
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param pt.scale Factor by which to multiply the size of the points (default: 5)
#' @param pt.shape Shape to use for the points (default = circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors.use.1 Color palette to use for feature 1.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param colors.use.2 Color palette to use for feature 2.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param bins Number of bins to use in dividing expression levels. (default: 10)
#' @param plot.height Plot height in pixels (default: 900)
#' @param plot.width Plot width in pixels (default: 900)
#' @param plot.title  Display title with the name of the feature? (default TRUE)
#' @param pt.info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param legend Display legend? (default: TRUE)
#' @param legend.font.size Legend font size (default: 12)
#' @param do.return Return the plot object instead of displaying it (default: FALSE)
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom Seurat GetDimReduction
#' @importFrom Seurat FetchData
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom grDevices colorRampPalette
#'
#' @return
#' @export
#'
#' @examples
Feature2Plotly <- function(seuratObj,
                           feature.1.use = NULL,
                           feature.2.use = NULL,
                           do.return = FALSE,
                           pt.scale = 1.5,
                           pt.shape = "circle",
                           opacity = 0.5,
                           reduction.use = "tsne",
                           dim.1 = 1,
                           dim.2 = 2,
                           colors.use.1 = "Reds",
                           colors.use.2 = "Blues",
                           bins = 10,
                           plot.height = "750",
                           plot.width = "750",
                           pt.info = NULL,
                           legend = TRUE,
                           legend.font.size = 12){

  df <- as.data.frame(GetDimReduction(object = seuratObj,
                                      reduction.type = reduction.use,
                                      slot = "cell.embeddings"))
  dim.code <- GetDimReduction(
    object = seuratObj,
    reduction.type = reduction.use,
    slot = "key"
  )

  dim.axes <- colnames(
    GetDimReduction(
      object = seuratObj,
      reduction.type = reduction.use,
      slot = "cell.embeddings"
    )
  )

  dim.code <- c(dim.axes[[dim.1]], dim.axes[[dim.2]])
  df <- df[,dim.code]
  cell_names <- rownames(df)

  df$x <- df[,1]
  df$y <- df[,2]

  #if (is.null(feature)){ stop("No gene or feature given") }
  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (colors.use.1 %in% rownames(brewer.pal.info)){
    pal.1 <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use.1,]$maxcolors,colors.use.1))(bins)
  } else if (colors.use.1 %in% viridis_palettes){
    pal.1 <- viridis(n = bins, option = colors.use.1)
  } else {
    pal.1 <- colors.use.1
  }

  if (colors.use.2 %in% rownames(brewer.pal.info)){
    pal.2 <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use.2,]$maxcolors,colors.use.2))(bins)
  } else if (colors.use.2 %in% viridis_palettes){
    pal.2 <- viridis(n = bins, option = colors.use.2)
  } else {
    pal.2 <- colors.use.2
  }

  feature.1.data <- FetchData(object = seuratObj,
                              vars.all = feature.1.use,
                              use.scaled = TRUE)
  size.1.data <- FetchData(object = seuratObj,
                           vars.all = feature.1.use,
                           use.scaled = FALSE)
  feature.1.data[,1][feature.1.data[,1] == 0] <- NA
  feature.1.data <- as.matrix(feature.1.data)
  cut.feature.1.data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature.1.data), breaks = bins)))
  df[,"feature.1"] <- cut.feature.1.data
  df[,"size.1"] <- size.1.data[,1] * pt.scale

  feature.2.data <- FetchData(object = seuratObj,
                              vars.all = feature.2.use,
                              use.scaled = TRUE)
  size.2.data <- FetchData(object = seuratObj,
                           vars.all = feature.2.use,
                           use.scaled = FALSE)
  feature.2.data[,1][feature.2.data[,1] == 0] <- NA
  feature.2.data <- as.matrix(feature.2.data)
  cut.feature.2.data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature.2.data), breaks = bins)))
  df[,"feature.2"] <- cut.feature.2.data
  df[,"size.2"] <- size.2.data[,1] * pt.scale

  if(!is.null(pt.info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt.info
      rowinfo = ""
      for(j in 1:length(pt.info)){
        rowinfo <- paste0(rowinfo, " </br> ", pt.info[j], ": ", seuratObj@meta.data[i, pt.info[j]])
      }
      rowinfo <- paste0(rowinfo, ' </br> Expr ', feature.1.use, ': ', df[i,"feature.1"])
      rowinfo <- paste0(rowinfo, ' </br> Expr ', feature.2.use, ': ', df[i,"feature.2"])
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- meta.info
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               color = ~feature.1,
               mode = 'markers',
               colors = c(pal.1,pal.2),
               size = ~size.1,
               sizes = c(0,max(df$size.1)),
               marker = list(symbol = pt.shape,
                             opacity = opacity,
                             sizemode = "diameter"
               ),
               width = plot.width,
               height = plot.height,
               type = 'scattergl',
               legendgroup = 1,
               showlegend = legend,
               name = feature.1.use) %>%
    add_trace(df,
              x = ~x,
              y = ~y,
              color = ~feature.2,
              mode = 'markers',
              type = 'scattergl',
              size = ~size.2,
              sizes = c(0,max(df$size.2)),
              marker = list(symbol = pt.shape,
                            opacity = opacity,
                            sizemode = "diameter"),
              legendgroup = 1,
              showlegend = legend,
              name = feature.2.use) %>%
    layout(title = paste(feature.1.use, feature.2.use, sep = " x "),
           xaxis = list(title = dim.axes[as.numeric(dim.1)]),
           yaxis = list(title = dim.axes[as.numeric(dim.2)]),
           margin = c(100,NA,NA,NA)
    )

  if(!is.null(pt.info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = ~meta.info,
                           showlegend = FALSE
    )
  }

  p <- p %>% layout(legend = list(
    font = list(
      size = legend.font.size)))

  if (isTRUE(do.return)){
    return(p)
  } else {
    p
  }
}

#' Feature2Plotly3D
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' seuratObj@@dr slot
#'
#' @param seuratObj Seurat object
#' @param feature.1.use First variable to display. Currently only works with gene names
#' @param feature.2.use Second variable to display. Currently only works with gene names
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param dim.3 Dimension to display on the z-axis (default: 3)
#' @param pt.scale Factor by which to multiply the size of the points (default: 5)
#' @param pt.shape Shape to use for the points (default = circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors.use.1 Color palette to use for feature 1.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param colors.use.2 Color palette to use for feature 2.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param bins Number of bins to use in dividing expression levels. (default: 10)
#' @param plot.height Plot height in pixels (default: 900)
#' @param plot.width Plot width in pixels (default: 900)
#' @param plot.title  Display title with the name of the feature? (default TRUE)
#' @param plot.axes Display the major x, y, and z axes? (default: FALSE)
#' @param plot.grid Display the major unit tick marks? (default: FALSE)
#' @param pt.info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param legend Display legend? (default: TRUE)
#' @param legend.font.size Legend font size (default: 12)
#' @param do.return Return the plot object instead of displaying it (default: FALSE)
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom Seurat GetDimReduction
#' @importFrom Seurat FetchData
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom grDevices colorRampPalette
#'
#' @return If do.return is TRUE, a plotly object.
#' @export
#'
#' @examples
Feature2Plotly3D <- function(seuratObj,
                             feature.1.use = NULL,
                             feature.2.use = NULL,
                             do.return = FALSE,
                             pt.scale = 0.5,
                             pt.shape = "circle",
                             opacity = 1,
                             reduction.use = "tsne",
                             dim.1 = 1,
                             dim.2 = 2,
                             dim.3 = 3,
                             colors.use.1 = "Reds",
                             colors.use.2 = "Blues",
                             bins = 10,
                             plot.height = "750",
                             plot.width = "750",
                             pt.info = NULL,
                             legend = TRUE,
                             legend.font.size = 12,
                             plot.grid = FALSE,
                             plot.axes = FALSE){

  df <- as.data.frame(GetDimReduction(object = seuratObj,
                                      reduction.type = reduction.use,
                                      slot = "cell.embeddings"))
  dim.code <- GetDimReduction(
    object = seuratObj,
    reduction.type = reduction.use,
    slot = "key"
  )

  dim.axes <- colnames(
    GetDimReduction(
      object = seuratObj,
      reduction.type = reduction.use,
      slot = "cell.embeddings"
    )
  )

  dim.code <- c(dim.axes[[dim.1]], dim.axes[[dim.2]], dim.axes[[dim.3]])
  df <- df[,dim.code]
  cell_names <- rownames(df)

  df$x <- df[,1]
  df$y <- df[,2]
  df$z <- df[,3]
  rownames(df) <- cell_names

  #if (is.null(feature)){ stop("No gene or feature given") }

  feature.1.data <- FetchData(object = seuratObj,
                              vars.all = feature.1.use,
                              use.scaled = TRUE)
  size.1.data <- FetchData(object = seuratObj,
                           vars.all = feature.1.use,
                           use.scaled = FALSE)
  feature.1.data <- as.matrix(feature.1.data)
  cut.feature.1.data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature.1.data), breaks = bins)))
  cut.feature.1.data <- as.factor(as.numeric(as.factor(x = cut(x = as.numeric(x = feature.1.data), breaks = bins))))
  df[,"feature.1"] <- cut.feature.1.data
  df[,"size.1"] <- size.1.data[,1] * pt.scale

  feature.2.data <- FetchData(object = seuratObj,
                              vars.all = feature.2.use,
                              use.scaled = TRUE)
  size.2.data <- FetchData(object = seuratObj,
                           vars.all = feature.2.use,
                           use.scaled = FALSE)
  feature.2.data[,1][feature.2.data[,1] == 0] <- NA
  feature.2.data <- as.matrix(feature.2.data)
  cut.feature.2.data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature.2.data), breaks = bins)))
  df[,"feature.2"] <- cut.feature.2.data
  df[,"size.2"] <- size.2.data[,1] * pt.scale

  if(!is.null(pt.info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt.info
      rowinfo = ""
      for(j in 1:length(pt.info)){
        rowinfo <- paste0(rowinfo, " </br> ", pt.info[j], ": ", seuratObj@meta.data[i, pt.info[j]])
      }
      rowinfo <- paste0(rowinfo, ' </br> Expr ', feature.1.use, ': ', df[i,"feature.1"])
      rowinfo <- paste0(rowinfo, ' </br> Expr ', feature.2.use, ': ', df[i,"feature.2"])
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- meta.info
  }

  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (colors.use.1 %in% rownames(brewer.pal.info)){
    pal.1 <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use.1,]$maxcolors,colors.use.1))(bins)
  } else if (colors.use.1 %in% viridis_palettes){
    pal.1 <- viridis(n = bins, option = colors.use.1)
  } else {
    pal.1 <- colors.use.1
  }

  if (colors.use.2 %in% rownames(brewer.pal.info)){
    pal.2 <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use.2,]$maxcolors,colors.use.2))(bins)
  } else if (colors.use.2 %in% viridis_palettes){
    pal.2 <- viridis(n = bins, option = colors.use.2)
  } else {
    pal.2 <- colors.use.2
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               z = ~z,
               color = ~feature.1,
               mode = 'markers',
               colors = c(pal.1,pal.2),
               size = ~size.1,
               sizes = c(0,max(df$size.1)),
               marker = list(symbol = pt.shape,
                             opacity = opacity,
                             sizemode = "diameter"),
               width = plot.width,
               height = plot.height,
               type = 'scatter3d',
               showlegend = legend) %>%
    add_trace(df,
              x = ~x,
              y = ~y,
              z = ~z,
              color = ~feature.2,
              mode = 'markers',
              type = 'scatter3d',
              size = ~size.2,
              sizes = c(0,max(df$size.2)),
              marker = list(symbol = pt.shape,
                            opacity = opacity,
                            sizemode = "diameter"),
              showlegend = legend) %>%
    layout(title = paste(feature.1.use, feature.2.use, sep = " x "),
           scene = list(
             aspectratio = list(x = 1,y = 1,z = 1),
             camera = list(
               center = list(x = 0,y = 0,z = 0),
               eye = list(x = 2,y = -1,z = 0.5),
               up = list(x = 0,y = 0,z = 1)),
             dragmode = "turnable",
             xaxis = list(title = dim.axes[as.numeric(dim.1)], type = "double", showgrid = plot.grid, visible = plot.axes),
             yaxis = list(title = dim.axes[as.numeric(dim.2)], type = "double", showgrid = plot.grid, visible = plot.axes),
             zaxis = list(title = dim.axes[as.numeric(dim.3)], type = "double", showgrid = plot.grid, visible = plot.axes),
             margin = c(100,NA,NA,NA)))

  if(!is.null(pt.info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = ~meta.info,
                           showlegend = FALSE
    )
  }
  p <- p %>% layout(legend = list(
    font = list(
      size = legend.font.size)))

  if (isTRUE(do.return)){
    return(p)
  } else {
    p
  }
}

#' retrieveGO
#'
#' Retrieve the HUGO Gene Nomenclature Committee names associated with a GeneOntology term.
#'
#' @param term GO ID to retreive names for.
#' @param mart Biomart database to use (default: 'ensembl')
#' @param dataset Biomart dataset to use (default: 'hsapiens_gene_ensembl')
#'
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @return A list of HCGN names.
#' @export
#'
#' @examples genes_to_plot <- retrieveGO('GO:0046774')
#'
retrieveGO <- function(term, mart = 'ensembl', dataset='hsapiens_gene_ensembl'){
  gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                     filters = 'go',
                     values = term,
                     mart = useMart(mart,dataset=dataset)
                     )
  return(gene.data)
}

#' GOBubblePlotly
#'
#' Produces a Bubble Plot for the genes of a given GO term.
#'
#' @param seuratObj Seurat object
#' @param go_term Gene Ontology term identifier (i.e. GO:0046774)
#' @param group.by Factor by which to group cells.  (default: ident)
#' @param plot.height Plot height in pixels. (default: 900)
#' @param plot.width Plot width in pixels. (default: 900)
#' @param filter A list of gene names to filter the GO term members against. (default: seuratObj@var.genes)
#' @param ...options to pass to BubblePlotly
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#'
#' @return
#' @export
#'
#' @examples
GOBubblePlotly <- function(seuratObj,
                           go_term,
                           group.by = "ident",
                           plot.height = 900,
                           plot.width = 900,
                           filter = seuratObj@var.genes,
                           ...){

  go_genes_to_plot <- retrieveGO(go_term) %>%
                      select(hgnc_symbol) %>%
                      distinct() %>%
                      filter(hgnc_symbol %in% filter)

  if(length(go_genes_to_plot) > 0){
    BubblePlotly(seuratObj,
                 genes.plot = go_genes_to_plot$hgnc_symbol,
                 group.by = group.by,
                 plot.height = plot.height,
                 plot.width = plot.width,
                 ...)
  } else {
    print("No genes for that term are expressed in the dataset.")
  }

}
