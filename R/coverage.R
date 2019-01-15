#' Read position and count data into table
#'
#' @param x path to file names
#'
#' @examples
#' "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project/viral_bg/neg"
#' a <- list_names(x)
#'
#' @export

list_names <- function(x) {
  location = x
  names = list.files(location, full.names = T)
}


#' Function to read in position and count data

read_file <- function(x) {
  df <- readr::read_tsv(x, col_names = c("chrom", "start", "end", "count", "normalized_count", "dinuc"))
  df$name <- basename(x)
  df
}

#' Load count and posiiton data into table by cell, virus, and time
#'
#' @param x list of file names
#' @return table with count data organized by cell, virus, and time
#'
#' @examples
#' coverage_table(a, read_file)
#'
#'
#' @export

coverage_table <- function(x) {
  purrr::map_df(x, read_file) %>%
  dplyr::mutate(name = str_replace(name, ".mhv.neg.dinuc.bg", "")) %>%
  separate(name, into = c('cell', 'virus', 'time'), sep = '_')
}


