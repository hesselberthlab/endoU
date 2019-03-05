#' Calculate dinucleotide frequencies observed from experimntal data
#'
#' @param x viral_cov_tbl
#'
#'
#' @examples
#'
#' exp_freqs <- calc_freqs(viral_cov_tbl)
#'
#'
#' function to generate counts of dinucleotides per grouped library

frequencies <- function(x) {

  ## make df of dinucleotides
  dinuc <- c("AA","GA", "CA", "UA","AC", "CC", "GC", "UC", "AG", "CG", "GG", "UG", "AU", "CU", "GU", "UU")
  df_dinuc <- data_frame(dinuc)


  ## summarize the data by dinuc
  df <- aggregate(cbind(count) ~ dinuc, data = x, sum)
  df <- left_join(df_dinuc, df) %>%
    mutate(frequency = count/sum(count)) %>%
    mutate(count_sum = sum(count))
}

#'
#'
#'
#' @export

calc_freqs <- function(x) {

  ## add group ID for each library
  df1 <- x %>%
    mutate(Group = group_indices(., cell, virus, time))

  ## map frequencies function to the df by group
  res <- df1 %>%
    split(.$Group) %>%
    map_dfr(frequencies, .id = "Group")

  df2 <- df1 %>%
    select(cell, virus, time, Group)

  ## join with original table to regain library type information
  res <- res %>%
    mutate(Group = as.character(Group)) %>%
    left_join(
      df2 %>%
        mutate(Group = as.character(Group)) %>%
        unique(.), by = "Group"
    )

  res
}


#' Calculate enrichment and significance of dinucleotides captured in the
#' experimental data compared to prevalence of dinucelotide in RNA sequence of interest
#'
#' @param obs_freq table contianing dincucleotide frequencies for RNA sequence of interest
#' @param exp_freq table containing dinucleotide frequencies calculated from experimental data
#'
#' @examples
#' combined_freqs <- calc_freq_enrichment(mhv_freqs, exp_freqs)
#'
#' @export
#'

calc_freq_enrichment <- function(observed_freq, experimental_freq) {

  #remove count column from observed frequencies table
  obs_freqs <- observed_freq

  #calculate frequency enrichment comparing experimental to observed
  combined_freqs <- experimental_freq %>%
    left_join(obs_freqs, by = "dinuc") %>%
    mutate(log2_enrichment = log2(frequency.x/frequency.y))

  combined_freqs

}


#' Calculate percent of dinucloetide capture per total umi corrected reads in the library
#'
#' @param filenames list of filenames and locations contianing the library information
#' @param table table summarized counts for each dinucleotide for each library
#'
#' @examples
#'
#' calc_dinuc_percent(filenames, combined_freqs)
#'
#' @export
#'


calc_dinuc_percent <- function(filenames, table){

  #read data and sum all counts per library
  all_reads <- coverage_table(filenames) %>%
    dplyr::select(count, cell, virus, time) %>%
    dplyr::group_by(cell, virus, time) %>%
    dplyr::mutate(total_umi_reads = sum(count)) %>%
    dplyr::select(-count) %>%
    unique()

  #calculate percent of reads captured at each dinucleotide
  res <- left_join(combined_freqs, all_reads, by = c("cell", "virus", "time")) %>%
    mutate(percent_reads = (count/total_umi_reads)*100)

  res

}

