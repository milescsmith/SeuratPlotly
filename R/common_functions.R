#' @title PrepDf
#'
#' @param seuratObj
#' @param reduction.use
#' @param dim.1
#' @param dim.2
#' @param dim.3
#' @param ident
#'
#' @import magrittr
#' @importFrom Seurat GetDimReduction FetchData
#'
#' @return
#' @export
#'
#' @examples
PrepDf <- function(seuratObj,
                   reduction.use,
                   dim.1 = 1,
                   dim.2 = 2,
                   dim.3 = NULL,
                   group.by) {
  titles <- letters[24:26]

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
        vars.all = group.by)[, 1] %>%
    as.factor()
  }

  colnames(df) <- letters[c(dim.1,
                            dim.2,
                            dim.3)]
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
#' @param seuratObj
#' @param pt.info
#' @param df
#'
#' @import magrittr
#' @importFrom glue glue
#'
#' @return
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
#' @param df
#' @param palette.use
#'
#' @return
#' @export
#'
#' @import paletteer
#' @import magrittr
#'
#' @examples
PrepPalette <- function(df, palette.use){

  bins = length(unique(df[,'ident']))

  if (palette.use %in% palettes_d_names$palette){
    color.package <- palettes_d_names$package[which(palette.use == palettes_d_names$palette)]
    pal <- paletteer_d(package = !!color.package,
                       palette = !!palette.use,
                       n = bins)
  } else if (palette.use %in% palettes_c_names$palette){
    color.package <- palettes_c_names$package[which(palette.use == palettes_c_names$palette)]
    pal <- paletteer_c(package = !!color.package,
                       palette = !!palette.use,
                       n = bins)
  } else {
    pal <- palette.use
  }

  return(pal)
}
