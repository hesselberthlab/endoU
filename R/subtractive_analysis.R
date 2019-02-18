#' Normalize count data using the signal that occurs in the abscence of RNaseL
#'
#' @param data viral_cov_table
#'
#' @examples
#'
#' Spread the data by cell type or virus type and normalized count
#'
#' spread_data(viral_cov_table)
#'
#'
#' @export
#'

spread_data <- function(viral_cov_tbl, key) {
  key <- enquo(key)
  res <-  dplyr::select(viral_cov_tbl, -count, -dinuc) %>%
    spread(!!key, norm_count)
  res[is.na(res)] <- 0

  res

  }


#'
#' @param data viral_cov_table
#' @param key column for comparision, either "cell" or "virus"
#' @param m_col Name of mutant cell type or viral mutant of interest
#'
#' @examples
#'
#' iter_types(
#'    viral_cov_table,
#'    "cell",
#'    "RNaseL")
#'
#'
#' @export
#'

iter_types <- function(df, key , m_col) {


  cell_names <- unique(df[[key]])
  cell_names <- cell_names[cell_names != m_col]
  cell_names

}


#'
#' @param data count and position data spread by cell and normalized count
#' @param x_col list or name of WT cell types
#' @param m_col name of mutant cell type
#' @param x_col2 optional argument to provide all WT cell comparisions individually instead of in a list
#'
#' @examples
#'
#' subtract_norm(spread_tbl, "RNaseL")
#'
#'

subtract_norm <- function(df, x_col, m_col, x_col2 = NULL) {

  varname <- paste0(x_col, "_norm")
  res <- df

  res[,varname] <- res[,x_col] - res[,m_col]

  if(!is.null(x_col2)){
    varname2 <- paste0(x_col2, "_norm")
    res[,varname2] <- res[,x_col2] - res[,m_col]
  }

  res
}

#' @param list list of cell or virus types to interate through for normalizing to mutant
#' @param data count and position data spread by cell or virus and normalized count
#' @param m_col name of mutant cell or virus type
#'
#' @examples
#'
#' subtract_norm_all(
#'   cell_names,
#'   cell_spread_tbl,
#'   "RNaseL")
#'
#'
#' @export


subtract_norm_all <- function(names, df, m_col){

for (i in names) {

    df <- subtract_norm(df, i, m_col)

}
  df
}



