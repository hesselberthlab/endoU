#' Calculate fold change in cleavage activity in the abscence of RNaseL or EndoU actviity
#'
#' @param data viral_cov_tbl
#' @param m_col Name of mutant cell type or viral mutant of interest
#' @param key column for comparision, either "cell" or "virus"
#'
#' @examples
#'
#' foldChange_all(viral_cov_table,
#'    "RNaseL",
#'    "cell"
#' )
#'
#'
#' Function to calculate log2 fold change comparing mutant/WT

foldChange <- function(spread_tbl, x_col, m_col, x_col2 = NULL) {

  varname <- paste0(x_col, "_log2FC")
  res <- spread_tbl

  res[,varname] <- log2(res[,m_col] / res[,x_col])

  if(!is.null(x_col2)){
    varname2 <- paste0(x_col2, "_log2FC")
    res[,varname2] <- log2(res[,m_col2] / res[,x_col])
  }

  res
}

#' @export

foldChange_all <- function(df, m_col, key) {

  spread_tbl <- df %>%
      dplyr::select(-count, -dinuc) %>%
      spread(!!key, norm_count)
  spread_tbl[is.na(spread_tbl)] <- 0

  iter_list <- endoU::iter_types(df, key, m_col)
  iter_list

  for (i in iter_list) {

      spread_tbl <- foldChange(spread_tbl, i, m_col)

  }

  spread_tbl

}


#' Idenitfy the top sites of cleavage by fold change analysis
#'
#' @param data fold_change_tbl
#' @param cols_gather range of column numbers containing fold change values
#' @param top number of top value positions to identify
#'
#' @examples
#'
#' top_sites(fold_change_tbl, cols_gather = 7:8, top = 20)
#'
#'
#'
#' @export

top_sites <- function(df, cols_gather, top) {

  res <- df %>%
    gather(sample_comp, foldChange, cols_gather) %>%
    dplyr::filter(foldChange != "NaN" & foldChange != "-Inf" ) %>%
    arrange(foldChange)

  top_pos <- res %>%
    dplyr::select(pos) %>%
    unique() %>%
    head(top)

  res_final <- left_join(top_pos, res)
  res_final
}

