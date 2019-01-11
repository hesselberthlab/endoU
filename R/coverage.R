#' Read position and count data into table
#'
#' @param x list of file names
#'
#' @examples
#' # "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project/viral_bg/neg"
#' list_names(x)
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

#' Load count and posiiton data into table by cell, virus, and time and merge with dinucleotide table
#'
#' @param x list of file names
#' @param y read_file function
#'
#' @examples
#' # path to file
#' coverage_table(x, y)
#'
#' @export

coverage_table <- function(x) {
  purrr::map_df(names, read_file) %>%
  mutate(name = str_replace(name, ".mhv.neg.dinuc.bg", "")) %>%
  separate(name, into = c('cell', 'virus', 'time'), sep = '_')
}



