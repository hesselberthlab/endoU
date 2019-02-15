#' Read position and count data into table
#'
#' @param x path to file names
#'
#' @examples
#' "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project/viral_bg/neg"
#' filenames <- list_names(x)
#'
#' @export

list_names <- function(path) {
  names <- list.files(path, full.names = TRUE)
}


#' Function to read in position and count data

read_file <- function(x) {
  tab <- suppressMessages(readr::read_tsv(
    x,
    col_names = c(
      "chrom", "pos", "count", "normalized_count"
    )
  ))
  tab$filename <- basename(x)
  tab
}

#' Load count and posiiton data into table by cell, virus, and time
#' @param x list of file names
#' @return table with count data organized by cell, virus, and time
#'
#'
#' @examples
#' coverage_table(filenames)
#'
#'
#' @export

coverage_table <- function(filenames) {
  purrr::map_df(filenames, read_file) %>%
    group_by(filename) %>%
    mutate(ref = str_split(filename, "\\.")[[1]][[1]]) %>%
    ungroup() %>%
    select(-filename) %>%
    separate(ref, into = c('cell', 'virus', 'time'), sep = '_')
}


#' Complete coverage table by filling in missing data and adding in dinucleotide information by position
#' @param x table with position and count data
#' @param y path to a table with dinucleotide data by position
#' @return table with position, count, and dinucleotide data
#'
#'
#' @examples
#' complete_table(table, dinuc_path)
#'
#'
#' @export

complete_table <- function(table, dinuc_path) {
  x <- complete(table, cell, virus, time, pos, fill = list(count = 1))
  y <- suppressMessages(readr::read_tsv(
    dinuc_path,
    col_names = c(
      "pos", "dinuc"
    )
  ))

  final_table <- inner_join(x, y, by = "pos")

  final_table
}

#' Normalize coverage table using the sum of aligned reads per library
#' @param x path to depth files containing all aligned reads in each libary
#' @param y complete table with position and count data for specific RNA of interest
#' @return table with position, count, normalized count, and dinucleotide data
#'
#'
#' @examples
#' complete_table(path, table)
#'
#'
#' @export

normalize_table <- function(path, table){

  filenames <- list_names(path)

  all_reads <- coverage_table(filenames) %>%
    dplyr::select(count, cell, virus, time) %>%
    dplyr::group_by(cell, virus, time) %>%
    dplyr::mutate(total_umi_reads = sum(count)) %>%
    dplyr::select(-count) %>%
    unique()

  final_table <- inner_join(table, all_reads, by = c("cell", "virus", "time")) %>%
    dplyr::mutate(norm_count = (count/total_umi_reads)*100) %>%
    dplyr::select(-chrom, -normalized_count, -total_umi_reads)

  final_table
}




