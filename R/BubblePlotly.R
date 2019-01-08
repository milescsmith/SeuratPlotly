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
