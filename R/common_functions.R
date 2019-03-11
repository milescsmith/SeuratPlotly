#' @title PrepDf
#'
#' @description Helper function to extract dimensional reduction data from a Seurat object
#' and format it into a data frame for use in plotting.
#'
#' @param object scRNA-seq data object
#' @param reduction dimensional reduction to extract
#' @param dim_1 Dimension one.  X-axis data.
#' @param dim_2 Dimension two.  Y-axis data.
#' @param dim_3 Dimension three.  Z-axis data.
#' @param grouping Identity information to use for grouping.
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
                          dim_1 = 1,
                          dim_2 = 2,
                          dim_3 = NULL,
                          grouping = NULL) {
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

  if (is.null(dim_3)) {
    dim.code <- c(
      dim.axes[[dim_1]],
      dim.axes[[dim_2]]
    )
  } else {
    dim.code <- c(
      dim.axes[[dim_1]],
      dim.axes[[dim_2]],
      dim.axes[[dim_3]]
    )
  }

  df <- df[, dim.code]
  cell_names <- rownames(df)

  column.titles <- letters[24:26]
  colnames(df) <- column.titles[c(
    dim_1,
    dim_2,
    dim_3
  )]
  rownames(df) <- cell_names

  if (!is.null(grouping)){
    if (grouping != "ident") {
      df$ident <- object$grouping %>%
        as.factor()
    } else if (grouping == "ident"){
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
                          dim_1 = 1,
                          dim_2 = 2,
                          dim_3 = NULL,
                          grouping = NULL) {
  df <- Embeddings(
    object = object,
    reduction = reduction) %>%
    as.data.frame()

  dim.axes <- colnames(df)

  if (is.null(dim_3)) {
    dim.code <- c(
      dim.axes[[dim_1]],
      dim.axes[[dim_2]]
    )
  } else {
    dim.code <- c(
      dim.axes[[dim_1]],
      dim.axes[[dim_2]],
      dim.axes[[dim_3]]
    )
  }

  df <- df[, dim.code]
  cell_names <- rownames(df)
  ident <- Idents(object) %>% as.factor()

  column.titles <- letters[24:26]
  colnames(df) <- column.titles[c(
    dim_1,
    dim_2,
    dim_3
  )]
  rownames(df) <- cell_names

  if (!is.null(grouping)){
    if (grouping != "ident") {
      df$ident <- object$grouping %>%
        as.factor()
    } else if (grouping == "ident"){
      df$ident <- Idents(object) %>% as.factor()
    }
  }
  df %<>% rownames_to_column("cell")

  return(df)
}

#' @title PrepInfo
#'
#' @description Using the members of pt_info, build a string for
#' each row that compiles the indicated pt_info from the meta.data
#' slot along with the necessary html tags and column names
#'
#' @param object Seurat object to extract data from
#' @param pt_info A list of meta.data columns to add to the hoverinfo popup.
#' @param df Plotting data frame to which to add the popup data.
#'
#' @importFrom stringr str_glue
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
PrepInfo.seurat <- function(object, pt_info, df) {
  if (!is.null(pt_info)) {
    meta.info <- list()
    # for each row
    for (i in seq(dim(df)[1])) {
      # for each member of pt_info
      rowinfo <- ""
      for (j in 1:length(pt_info)) {
        rowinfo <- str_glue("{rowinfo} </br> {pt_info[j]}: {object@meta.data[i, pt_info[j]]}")
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
PrepInfo.Seurat <- function(object, pt_info, df) {
  if (!is.null(pt_info)) {
    feature_info <- FetchData(object = object, vars = pt_info) %>% rownames_to_column("cell")
    # for each row
    if ("cell" %nin% colnames(df)){
      df %<>% as.data.frame() %>% rownames_to_column("cell")
    }
    meta.info <- new("data.frame")
    for (i in seq(nrow(feature_info))) {
      # for each member of pt_info
      rowinfo <- ""
      for (j in 1:length(pt_info)) {
        if (is.numeric(feature_info[i, pt_info[j]])){
          rowinfo <- str_glue("{rowinfo} </br> {pt_info[j]}: {round(feature_info[i, pt_info[j]],digits = 3)}")
        } else {
          rowinfo <- str_glue("{rowinfo} </br> {pt_info[j]}: {feature_info[i, pt_info[j]]}")
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
#' @param palette_use The name of a palette to use.  Must be a palette available in
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
PrepPalette <- function(df,
                        palette_use) {
  bins <- length(unique(df[, "ident"]))

  if (palette_use %in% palettes_d_names$palette) {
    color.package <- palettes_d_names$package[which(palette_use == palettes_d_names$palette)]
    pal <- paletteer_d(
      package = !!color.package,
      palette = !!palette_use
    )
  } else if (palette_use %in% palettes_c_names$palette) {
    color.package <- palettes_c_names$package[which(palette_use == palettes_c_names$palette)]
    pal <- paletteer_c(
      package = !!color.package,
      palette = !!palette_use
    )
  } else {
    pal <- palette_use
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
#' @param feature feature (genes expression, metadata, etc...) to retrieve data for
#' @param assay_use Assay to pull feature data from.  Default: "RNA"
#' @param slot_use Slot to pull feature data from. Default: "data"
#'
#' @importFrom stringr str_glue
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
#' @import Seurat
#' @return data.frame
GetFeatureValues.seurat <- function(object,
                                    df,
                                    feature,
                                    bins = NULL,
                                    ...) {
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
  feature_data <- FetchData(object,
                            vars.all = feature,
                            use.scaled = use.scaled,
                            use.raw = use.raw) %<>% rownames_to_column("cell")

  colnames(feature_data) %<>% str_remove(pattern = key)

  if (!is.null(bins)){
    feature_data %<>% mutate_if(is.numeric,
                                list(~cut(.,
                                          breaks = bins))) %>%
      mutate_if(is.factor,
                list(as.numeric))
  }
  feature_data
  df %<>% inner_join(feature_data, by = "cell")
  return(df)
}

# Seurat 3
#' @rdname PrepInfo
#' @method PrepInfo Seurat
#' @import Seurat
#' @importFrom stringr str_glue
#' @importFrom purrr map_chr
#' @return data.frame
GetFeatureValues.Seurat <- function(object,
                                    df,
                                    feature,
                                    assay_use = "RNA",
                                    slot_use = "data",
                                    bins = NULL,
                                    suffix = NULL,
                                    ...) {
  key <- GetAssayData(object = object, assay = assay_use, slot = "key")
  feature_data <- FetchData(object,
                            vars = str_glue("{key}{feature}") %>% as.character(),
                            slot = slot_use) %>%
    rownames_to_column("cell")

  if (!is.null(suffix)){
    colnames(feature_data)[2] %<>% str_remove(pattern = key) %>% str_glue("_{suffix}")
  } else {
    colnames(feature_data)[2] %<>% str_remove(pattern = key)
  }


  if (!is.null(bins)){
    feature_data %<>% mutate_if(is.numeric,
                                list(~cut(.,
                                         breaks = bins))) %>%
      mutate_if(is.factor,
                list(as.numeric))
  }

  df %<>% inner_join(feature_data, by = "cell")
  return(df)
}

#' inverse match
#'
#' See \code{Hmisc::\link[Hmisc]{\%nin\%}} for details.
#'
#' @name %nin%
#' @rdname nin
#' @keywords internal
#' @export
#' @importFrom Hmisc %nin%
#' @usage lhs \%nin\% rhs
NULL
