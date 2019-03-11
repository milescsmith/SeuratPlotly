#' @title PrepDf
#'
#' @description Helper function to extract dimensional reduction data from a Seurat object
#' and format it into a data frame for use in plotting.
#'
#' @param object scRNA-seq data object
#' @param reduction dimensional reduction to extract
#' @param dim.1 Dimension one.  X-axis data.
#' @param dim.2 Dimension two.  Y-axis data.
#' @param dim.3 Dimension three.  Z-axis data.
#' @param group.by Identity information to use for grouping.
#'   Permissible values include meta.data column names.
#'
#' @importFrom tibble rownames_to_column
#' @export
#'
#' @examples
PrepDf <- function(object, ...) {
  UseMethod("PrepDf")
}

# Seurat 2 objects
#' @rdname PrepDf
#' @import Seurat
#' @method PrepDf seurat
#' @return data.frame
PrepDf.seurat <- function(object,
                   reduction,
                   dim.1 = 1,
                   dim.2 = 2,
                   dim.3 = NULL,
                   group.by = NULL) {
  df <- GetDimReduction(
    object = object,
    reduction.type = reduction,
    slot = "cell.embeddings"
  ) %>%
    as.data.frame()

  dim.code <- GetDimReduction(
    object = object,
    reduction.type = reduction,
    slot = "key"
  )

  dim.axes <- colnames(
    GetDimReduction(
      object = object,
      reduction.type = reduction,
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

  column.titles <- letters[24:26]
  colnames(df) <- column.titles[c(
    dim.1,
    dim.2,
    dim.3
  )]
  rownames(df) <- cell_names

  if (!is.null(group.by)){
    if (group.by != "ident") {
      df$ident <- object$group.by %>%
        as.factor()
    } else if (group.by == "ident"){
      df$ident <- Idents(object) %>% as.factor()
    }
  }

  df %<>% rownames_to_column("cell")

  return(df)
}

# Seurat 3
#' @rdname PrepDf
#' @import Seurat
#' @method PrepDf Seurat
#' @return data.frame
PrepDf.Seurat <- function(object,
                          reduction,
                          dim.1 = 1,
                          dim.2 = 2,
                          dim.3 = NULL,
                          group.by = NULL) {
  df <- Embeddings(
    object = object,
    reduction = reduction) %>%
    as.data.frame()

  dim.axes <- colnames(df)

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
  ident <- Idents(object) %>% as.factor()

  column.titles <- letters[24:26]
  colnames(df) <- column.titles[c(
    dim.1,
    dim.2,
    dim.3
  )]
  rownames(df) <- cell_names

  if (!is.null(group.by)){
    if (group.by != "ident") {
      df$ident <- object$group.by %>%
        as.factor()
    } else if (group.by == "ident"){
      df$ident <- Idents(object) %>% as.factor()
    }
  }
  df %<>% rownames_to_column("cell")

  return(df)
}

#' @title PrepInfo
#'
#' @description Using the members of pt.info, build a string for
#' each row that compiles the indicated pt.info from the meta.data
#' slot along with the necessary html tags and column names
#'
#' @param object Seurat object to extract data from
#' @param pt.info A list of meta.data columns to add to the hoverinfo popup.
#' @param df Plotting data frame to which to add the popup data.
#'
#' @importFrom glue glue
#'
#' @return data.frame
#' @export
#'
#' @examples

PrepInfo <- function(object, ...) {
  UseMethod("PrepInfo")
}

# Seurat 2
#' @rdname PrepInfo
#' @method PrepInfo seurat
#' @return data.frame
PrepInfo.seurat <- function(object, pt.info, df) {
  if (!is.null(pt.info)) {
    meta.info <- list()
    # for each row
    for (i in seq(dim(df)[1])) {
      # for each member of pt.info
      rowinfo <- ""
      for (j in 1:length(pt.info)) {
        rowinfo <- glue("{rowinfo} </br> {pt.info[j]}: {object@meta.data[i, pt.info[j]]}")
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

# Seurat 3
#' @rdname PrepInfo
#' @method PrepInfo Seurat
#' @return data.frame
PrepInfo.Seurat <- function(object, pt.info, df) {
  if (!is.null(pt.info)) {
    feature_info <- FetchData(object = object, vars = pt.info) %>% rownames_to_column("cell")
    # for each row
    if ("cell" %nin% colnames(df)){
      df %<>% as.data.frame() %>% rownames_to_column("cell")
    }
    meta.info <- new("data.frame")
    for (i in seq(nrow(feature_info))) {
      # for each member of pt.info
      rowinfo <- ""
      for (j in 1:length(pt.info)) {
        if (is.numeric(feature_info[i, pt.info[j]])){
          rowinfo <- glue("{rowinfo} </br> {pt.info[j]}: {round(feature_info[i, pt.info[j]],digits = 3)}")
        } else {
          rowinfo <- glue("{rowinfo} </br> {pt.info[j]}: {feature_info[i, pt.info[j]]}")
        }

      }
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- meta.info
  } else {
    df$meta.info <- Idents(object)
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

#' @title GetFeatureValues
#'
#' @description Get expression levels of a given feature for each cell
#'
#' @param object Seurat object to get feature data from
#' @param df Plotting data frame to which to add the popup data.
#' @param features Features (genes expression, metadata, etc...) to retrieve data for
#' @param assay.use Assay to pull feature data from.  Default: "RNA"
#' @param slot.use Slot to pull feature data from. Default: "data"
#'
#' @importFrom glue glue
#'
#' @return data.frame
#' @export
#'
#' @examples

GetFeatureValues <- function(object, ...) {
  UseMethod("GetFeatureValues")
}

# Seurat 2
#' @rdname PrepInfo
#' @method PrepInfo seurat
#' @importFrom Seurat
#' @return data.frame
GetFeatureValues.seurat <- function(object,
                                    df,
                                    features,
                                    assay.use = "RNA",
                                    slot.use = "data") {
  if (slot.use == "scaled.data"){
    use.scaled <- TRUE
    use.raw <- FALSE
  } else if (slot.use == "counts" | slot.use == "raw.data"){
    use.raw <- FALSE
    use.scaled <- FALSE
  } else {
    use.raw <- FALSE
    use.scaled <- FALSE
  }
  feature_data <- FetchData(object, vars.all = features, use.scaled = use.scaled, use.raw = use.raw) %>% rownames_to_column("cell")
  df %<>% inner_join(feature_data, by = "cell")
  return(df)
}

# Seurat 3
#' @rdname PrepInfo
#' @method PrepInfo Seurat
#' @importFrom Seurat
#' @return data.frame
GetFeatureValues.Seurat <- function(object,
                                    df,
                                    features,
                                    assay.use = "RNA",
                                    slot.use = "data") {
  DefaultAssay(object) <- assay.use
  feature_data <- FetchData(object, vars = features, slot = slot.use) %>% rownames_to_column("cell")
  df %<>% inner_join(feature_data, by = "cell")
  return(df)
}
