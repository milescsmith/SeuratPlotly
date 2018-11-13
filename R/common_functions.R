#' @title PrepDf
#'
#' @description Helper function to extract dimensional reduction data from a Seurat object
#' and format it into a data frame for use in plotting.
#'
#' @param seuratObj Seurat object
#' @param reduction.use dimensional reduction to extract
#' @param dim.1 Dimension one.  X-axis data.
#' @param dim.2 Dimension two.  Y-axis data.
#' @param dim.3 Dimension three.  Z-axis data.
#' @param group.by Identity information to use for grouping.
#'   Permissible values include meta.data column names.
#'
#' @import magrittr
#' @importFrom Seurat GetDimReduction FetchData
#'
#' @return data.frame
#' @export
#'
#' @examples
PrepDf <- function(seuratObj,
                   reduction.use,
                   dim.1 = 1,
                   dim.2 = 2,
                   dim.3 = NULL,
                   group.by) {
  df <- GetDimReduction(
    object = seuratObj,
    reduction.type = reduction.use,
    slot = "cell.embeddings"
  ) %>%
    as.data.frame()

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

  if (is.null(dim.3)) {
    dim.code <- c(
      dim.axes[[dim.1]],
      dim.axes[[dim.2]]
    )
  } else {
    dim.code <- c(
      dim.axes[[dim.1]],
      dim.axes[[dim.2]],
      dim.axes[[dim.3]]
    )
  }

  df <- df[, dim.code]
  cell_names <- rownames(df)
  ident <- seuratObj@ident %>% as.factor()
  if (group.by != "ident") {
    ident <- FetchData(
      object = seuratObj,
      vars.all = group.by
    )[, 1] %>%
      as.factor()
  }

  column.titles <- letters[24:26]
  colnames(df) <- column.titles[c(
    dim.1,
    dim.2,
    dim.3
  )]
  rownames(df) <- cell_names
  df$ident <- ident

  return(df)
}

#' @title PrepInfo
#'
#' @description Using the members of pt.info, build a string for
#' each row that compiles the indicated pt.info from the meta.data
#' slot along with the necessary html tags and column names
#'
#' @param seuratObj Seurat object to extract data from
#' @param pt.info A list of meta.data columns to add to the hoverinfo popup.
#' @param df Plotting data frame to which to add the popup data.
#'
#' @import magrittr
#' @importFrom glue glue
#'
#' @return data.frame
#' @export
#'
#' @examples
PrepInfo <- function(seuratObj, pt.info, df) {
  if (!is.null(pt.info)) {
    meta.info <- list()
    # for each row
    for (i in seq(dim(df)[1])) {
      # for each member of pt.info
      rowinfo <- ""
      for (j in 1:length(pt.info)) {
        rowinfo <- glue("{rowinfo} </br> {pt.info[j]}: {seuratObj@meta.data[i, pt.info[j]]}")
      }
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- meta.info
  } else {
    df$meta.info <- df$ident
  }
  return(df)
}

#' @title PrepPalette
#'
#' @description Given a palette name and data frame with clustering or grouping
#' information in an 'ident' column, return a palette containing a color for each
#' unique cluster identity.
#'
#' @param df A graphing data.frame (from PrepDf) containing a column named 'ident'
#'   with group identities.
#' @param palette.use The name of a palette to use.  Must be a palette available in
#'   the Paletteer package.  If there are more unique identities than colors in the
#'   palette, additional values will be created by interpolation.
#'
#' @return A list containing color values.
#' @export
#'
#' @import paletteer
#' @import magrittr
#' @importFrom grDevices colorRampPalette
#'
#' @examples
PrepPalette <- function(df, palette.use) {
  bins <- length(unique(df[, "ident"]))

  if (palette.use %in% palettes_d_names$palette) {
    color.package <- palettes_d_names$package[which(palette.use == palettes_d_names$palette)]
    pal <- paletteer_d(
      package = !!color.package,
      palette = !!palette.use
    )
  } else if (palette.use %in% palettes_c_names$palette) {
    color.package <- palettes_c_names$package[which(palette.use == palettes_c_names$palette)]
    pal <- paletteer_c(
      package = !!color.package,
      palette = !!palette.use
    )
  } else {
    pal <- palette.use
  }
  pal <- colorRampPalette(pal)(bins)
  return(pal)
}
