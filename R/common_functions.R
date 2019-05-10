#' @title PrepDr
#'
#' @description Helper function to extract dimensional reduction data from a Seurat object
#' and format it into a data frame for use in plotting.
#'
#' @param object scRNA-seq data object
#' @param reduction dimensional reduction to extract
#' @param dim_1 Dimension one.  X-axis data.
#' @param dim_2 Dimension two.  Y-axis data.
#' @param dim_3 Dimension three.  Z-axis data.
#' @param grouping_var Identity information to use for grouping_var.
#' Permissible values include meta.data column names.
#' @param ... Not used
#'
#' @importFrom tibble as_tibble
#' @export
#'
PrepDr <- function(object, ...) {
  UseMethod("PrepDr")
}

# Seurat 3
#' @rdname PrepDr
#' @import Seurat
#' @method PrepDr Seurat
#' @return data.frame
#'
PrepDr.Seurat <- function(object,
                          reduction,
                          dim_1 = 1,
                          dim_2 = 2,
                          dim_3 = NULL,
                          grouping_var = NULL,
                          ...) {

  df <- Embeddings(object = object,
                   reduction = reduction) %>%
    as_tibble(rownames = "cell")

  column.titles <- letters[24:(22+ncol(df))]
  colnames(df)[2:ncol(df)] <- column.titles

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

  if (!is.null(grouping_var)){
    if (grouping_var != "ident") {
      df[[ident]] <- object[[grouping_var]] %>%
        as.factor()
    } else if (grouping_var == "ident"){
      df[[ident]] <- Idents(object) %>% as.factor()
    }
  }

  return(df)
}

#' @title PrepPalette
#'
#' @description Given a palette name and data frame with clustering or grouping_var
#' information in an 'ident' column, return a palette containing a color for each
#' unique cluster identity.
#'
#' @param df A graphing data.frame (from PrepDr) containing a column named 'ident'
#' with group identities.
#' @param palette The name of a palette to use.  Must be a palette available in
#' the Paletteer package.  If there are more unique identities than colors in the
#' palette, additional values will be created by interpolation.
#'
#' @return A list containing color values.
#' @export
#'
#' @import paletteer
#' @importFrom grDevices colorRampPalette
#'
PrepPalette <- function(df,
                        palette) {
  bins <- length(unique(df[["ident"]]))

  if (palette %in% palettes_d_names[[palette]]) {
    color.package <- palettes_d_names$package[which(palette == palettes_d_names$palette)]
    pal <- paletteer_d(
      package = !!color.package,
      palette = !!palette
    )
  } else if (palette %in% palettes_c_names[[palette]]) {
    color.package <- palettes_c_names$package[which(palette == palettes_c_names$palette)]
    pal <- paletteer_c(
      package = !!color.package,
      palette = !!palette
    )
  } else {
    pal <- palette
  }
  pal <- colorRampPalette(pal)(bins)
  return(pal)
}

#' @title PrepQuantitativePalette
#'
#' @description Given a palette name and data frame with clustering or grouping_var
#' information in an 'ident' column, return a palette containing a color for each
#' unique cluster identity.
#'
#' @param bins Number of discrete values that need colors.
#' @param palette The name of a palette to use.  Must be a palette available in
#' the Paletteer package.  If there are more unique identities than colors in the
#' palette, additional values will be created by interpolation.
#'
#' @return A list containing color values.
#' @export
#'
#' @importFrom paletteer paletteer_d paletteer_c
#' @importFrom grDevices colorRampPalette
#'
PrepQuantitativePalette <- function(bins,
                                   palette) {

  if (palette %in% palettes_d_names[[palette]]) {
    color.package <- palettes_d_names$package[which(palette == palettes_d_names$palette)]
    pal <- paletteer_d(
      package = !!color.package,
      palette = !!palette
    )
  } else if (palette %in% palettes_c_names[[palette]]) {
    color.package <- palettes_c_names$package[which(palette == palettes_c_names$palette)]
    pal <- paletteer_c(
      package = !!color.package,
      palette = !!palette
    )
  } else if (palette %in% palettes_dynamic_names[[palette]]) {
    color.package <- palettes_dynamic_names$package[which(palette == palettes_dynamic_names$palette)]
    pal <- paletteer_dynamic(
      package = !!color.package,
      palette = !!palette
    )

  } else {
    pal <- palette
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
#' @param features Feature (genes expression, metadata, etc...) to retrieve data for
#' @param assay Assay to pull feature data from.  Default: NULL
#' @param slot Slot to pull feature data from. Default: NULL
#' @param bins If provided, divide the values in the given number of bins.
#' @param suffix If given, add this suffix to the end of cells returned.
#' @param ... not used
#'
#' @importFrom stringr str_glue str_remove str_detect
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map_chr
#' @importFrom dplyr inner_join mutate_if
#'
#' @return data.frame
#' @export
#'
GetFeatureValues <- function(object, ...) {
  UseMethod("GetFeatureValues")
}


# Seurat 3
#' @rdname GetFeatureValues
#' @method GetFeatureValues Seurat
#' @importFrom Seurat FetchData
#' @return data.frame
#'
GetFeatureValues.Seurat <- function(object,
                                    features,
                                    df = NULL,
                                    assay = NULL,
                                    slot = NULL,
                                    bins = NULL,
                                    suffix = NULL,
                                    ...) {
  if (!is.null(assay)){
    feature_data <- FetchData(object,
                            vars = map_chr(features, function(x) str_glue("{assay}_{x}")),
                            slot = slot) %>%
    rownames_to_column("cell")

    colnames(feature_data) %<>% str_remove(pattern = str_glue("{assay}_"))
    if (!is.null(suffix)){
      colnames(feature_data) %<>% map_chr(function(i){
        if (!str_detect(string = i, pattern = "cell")){
          i <- str_glue("{i}_{suffix}")
        }
      })
    }
  } else {
    feature_data <- FetchData(object,
                              vars = features) %>%
    rownames_to_column("cell")

    if (!is.null(suffix)){
      colnames(feature_data) %<>% map_chr(function(i){
        if (!str_detect(string = i, pattern = "cell")){
          i <- str_glue("{i}_{suffix}")
        }
      })
    }
  }

  if (!is.null(bins)){
    feature_data %<>% mutate_if(is.numeric,
                                list(~cut(.,
                                          breaks = bins,
                                          labels = FALSE)-1)) %>%
      mutate_if(is.factor,
                list(as.numeric))
  }

  if (!is.null(df)){
    df %<>% inner_join(feature_data, by = "cell")
  } else {
    df <- feature_data
  }
  return(df)
}


#' @title PrepInfo
#'
#' @description Get meta data for data points
#'
#' @param object Seurat object to get feature data from.
#' @param pt_info Meta data columns to retreive data from.
#' @param df Plotting data frame to which to add the meta data.
#' @param ... not used
#'
#' @importFrom stringr str_glue str_remove str_detect
#' @importFrom tibble rownames_to_column
#' @importFrom purrr map_chr
#' @importFrom dplyr inner_join mutate_if
#'
#' @return
#' @export
#'
PrepInfo <- function(object, ...) {
  UseMethod("PrepInfo")
}

# Seurat 3
#' @rdname PrepInfo
#' @method PrepInfo Seurat
#' @import Seurat
#' @return data.frame
#'
PrepInfo.Seurat <- function(object,
                            pt_info,
                            df){

  intersecting_columns <- intersect(colnames(object[[]]),
                                    pt_info)
  if (length(intersecting_columns) != length(pt_info)){
    message(str_glue("The columns {setdiff(pt_info, intersecting_columns)} were not found"))
  }


}


















