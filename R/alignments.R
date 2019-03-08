#' Calcuate total reads per library and total reads per RNA type
#'
#' @param x list object of filenames containing depth files
#'
#' @examples
#' count_table <- sum_all_reads(filenames)
#'
#'
#' functions needed to read in files and assign RNA type


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


library_table <- function(filenames) {
  purrr::map_df(filenames, read_file) %>%
    group_by(filename) %>%
    mutate(ref = str_split(filename, "\\.")[[1]][[1]]) %>%
    ungroup() %>%
    select(-filename) %>%
    separate(ref, into = c('cell', 'virus', 'time'), sep = '_') %>%
    separate(chrom, into = c('RNA_type', 'extra'), sep = '\\|') %>%
    dplyr::select(-extra)
}


#' @export
#'

sum_all_reads <- function(filenames) {

  all_reads <- library_table(filenames)

  total_table <- all_reads %>%
    dplyr::group_by(cell, virus, time) %>%
    dplyr::mutate(total_umi_reads = sum(count)) %>%
    dplyr::select(cell, virus, time, total_umi_reads) %>%
    unique()

  rna_type_reads <- function(x) {

    reads <- x %>%
      group_by(cell, virus, time) %>%
      mutate(total_count = sum(count)) %>%
      select(cell, virus, time, total_count) %>%
      unique()

  }

  res <- all_reads %>%
    split(.$RNA_type) %>%
    map_dfr(rna_type_reads, .id = "RNA_type") %>%
    mutate(RNA_type = paste(RNA_type, "total", sep="_")) %>%
    spread(key=RNA_type, value=total_count) %>%
    left_join(total_table)

  res

}


#' Calcuate total reads per library and total reads per RNA type
#'
#' @param x list object of filenames containing depth files
#' @param cold_ids column numbers of count data for RNA types to compare to total reads in library
#'
#' @examples
#'
#' final_table <- percent_total(count_table, 4:11)
#'
#' @export
#'

percent_total <- function(x, col_ids){
  percent_cdna <- x %>%
    mutate_at(.fun = funs(percent = ./total_umi_reads * 100), .vars = col_ids)
  return(percent_cdna)
}














