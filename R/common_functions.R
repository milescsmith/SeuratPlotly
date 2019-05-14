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
#' @importFrom dplyr rename mutate inner_join select
#' @export
#'
PrepDr <- function(object, ...) {
  UseMethod("PrepDr")
}

# Seurat 3
#' @rdname PrepDr
#' @importFrom Seurat Embeddings Idents
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

  cell <- NULL
  value <- NULL

  df <- Embeddings(object = object,
                   reduction = reduction) %>%
    as_tibble(rownames = "cell") %>%
    select(cell,
           x = dim_1 + 1,
           y = dim_2 + 1,
           z = dim_3 + 1)

  if (is.null(grouping_var)){
    grouping_var <- "ident"
  }

  if (grouping_var != "ident") {
    idents <- object[[grouping_var]]
  } else if (grouping_var == "ident"){
    idents <- Idents(object)
  }

  idents %<>%
    as_tibble(rownames = "cell") %>%
    mutate(value = as.factor(value)) %>%
    rename(ident = value)

  df %<>% inner_join(idents)

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
#' @importFrom paletteer paletteer_d paletteer_c
#' @importFrom grDevices colorRampPalette
#'
PrepPalette <- function(df,
                        palette) {
  bins <- length(unique(df[["ident"]]))

  if (palette %in% palettes_d_names[["palette"]]) {
    color_package <- palettes_d_names$package[which(palette == palettes_d_names$palette)]
    pal <- paletteer_d(
      package = !!color_package,
      palette = !!palette
    )
  } else if (palette %in% palettes_c_names[["palette"]]) {
    color_package <- palettes_c_names$package[which(palette == palettes_c_names$palette)]
    pal <- paletteer_c(
      package = !!color_package,
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
#' @importFrom paletteer paletteer_d paletteer_c paletteer_dynamic
#' @importFrom grDevices colorRampPalette
#'
PrepQuantitativePalette <- function(bins,
                                   palette) {

  if (palette %in% palettes_d_names[["palette"]]) {
    color_package <- palettes_d_names$package[which(palette == palettes_d_names$palette)]
    pal <- paletteer_d(
      package = !!color_package,
      palette = !!palette
    )
  } else if (palette %in% palettes_c_names[["palette"]]) {
    color_package <- palettes_c_names$package[which(palette == palettes_c_names$palette)]
    pal <- paletteer_c(
      package = !!color_package,
      palette = !!palette,
      n = bins
    )
  } else if (palette %in% palettes_dynamic_names[["palette"]]) {
    color_package <- palettes_dynamic_names$package[which(palette == palettes_dynamic_names$palette)]
    pal <- paletteer_dynamic(
      package = !!color_package,
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
#' @importFrom tibble as_tibble
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
    as_tibble(rownames = "cell")

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
    as_tibble(rownames = "cell")
  }

  if (!is.null(suffix)){
    colnames(feature_data) <- map_chr(colnames(feature_data), function(i){
        if (!str_detect(string = i, pattern = "cell")){
          i <- str_glue("{i}_{suffix}")
        }
        i
      })
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
#' @description Using the members of pt_info, build a string for
#' each row that compiles the indicated pt_info from the meta.data
#' slot along with the necessary html tags and column names
#'
#' @param object Seurat object to extract data from
#' @param pt_info A list of meta.data columns to add to the hoverinfo popup.
#' @param df Plotting data frame to which to add the popup data.
#' @param ... Unused
#'
#' @importFrom tibble as_tibble
#' @importFrom stringr str_glue
#' @importFrom purrr map_chr
#'
#' @return data.frame
#' @export
#'
PrepInfo <- function(object, ...) {
  UseMethod("PrepInfo")
}

# Seurat 3
#' @rdname PrepInfo
#' @method PrepInfo Seurat
#' @importFrom Seurat FetchData
#' @return data.frame
#'
PrepInfo.Seurat <- function(object, pt_info, df) {
  if (!is.null(pt_info)) {
    feature_info <- FetchData(object = object, vars = pt_info) %>% as_tibble(rownames = "cell")
    # for each row
    if ("cell" %nin% colnames(df)){
      df %<>% as_tibble(rownames = "cell")
    }

    meta_info <- map_chr(seq(nrow(feature_info)), function(i) {
      # for each member of pt_info
      rowinfo <- ""
      for (j in 1:length(pt_info)) {
        if (is.numeric(feature_info[i, pt_info[j]])){
          rowinfo <- str_glue("{rowinfo} </br> {pt_info[j]}: {round(feature_info[i, pt_info[j]],digits = 3)}")
        } else {
          rowinfo <- str_glue("{rowinfo} </br> {pt_info[j]}: {feature_info[i, pt_info[j]]}")
        }
      }
      return(rowinfo)
    })

    df[["meta_info"]] <- meta_info
  } else {
    df[["meta_info"]] <- Idents(object)
  }
  return(df)
}
