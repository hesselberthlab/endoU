#' Calculate dinucleotide frequencies observed from experimntal data
#'
#' @param x viral_cov_tbl
#'
#'
#' @examples
#'
#' obs_freqs <- calc_freqs(viral_cov_tbl)
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


#' Calculate enrichment of dinucleotides captured in the
#' experimental data compared to prevalence of dinucelotide in RNA sequence of interest
#'
#' @param obs_freq table contianing dincucleotide frequencies for RNA sequence of interest
#' @param exp_freq table containing dinucleotide frequencies calculated from experimental data
#'
#' @examples
#' combined_freqs <- calc_freq_enrichment(mhv_freq, obs_freqs)
#'
#' @export
#'

calc_freq_enrichment <- function(expected_freq, observed_freq) {

  #calculate frequency enrichment comparing observed to expected

  combined_freqs <- observed_freq %>%
    left_join(expected_freq, by = "dinuc") %>%
    mutate(log2_enrichment = log2(frequency.x/frequency.y))

  combined_freqs

}

#' Calculate significance of enrichment values comparing dinucleotides captured in the
#' experimental data compared to prevalence of dinucelotide in RNA sequence of interest
#'
#' @param x table generated from calculating dinucelotide enrichment
#'
#'
#' @examples
#' cal_enrichment_sigs(combined_freqs)
#'
#'
#'
#' function to execute fisher exact test by row of dataframe comapring the odds ratio
#' of obtaining a speciic dinucleotide in the expected data to the observed data

fisher.byrow <- function(x) {
  mx <- matrix(
    c(
      x$count.y,
      (x$count_sum.y - x$count.y),
      x$count.x,
      (x$count_sum.x - x$count.x)),
      2, 2
  )

  fisher.test(mx, alternative="less")$p.value
}

#'
#' @importFrom qvalue qvalue
#'
#' @export
#'
#'

cal_enrichment_sigs <- function(x) {

  ## turn counts data into lists

  tbl <- x %>%
    select(Group, dinuc, starts_with("count")) %>%
    group_by(Group, dinuc) %>%
    nest()

  ## apply fisher test by row and join back with original data frame

  tbl2 <- mutate(tbl, pval = map_dbl(tbl$data, fisher.byrow)) %>%
          select(-data) %>%
          left_join(x)

  ## calculate qvals

  res <- mutate(tbl2, qval = qvalue::qvalue(tbl2$pval, pi0 = 1)$qvalues)

  res
}

usethis::use_package("qvalue")

#' Calculate percent of dinucloetide capture per total umi corrected reads in the library
#'
#' @param filenames list of filenames and locations contianing the library information
#' @param table table summarized counts for each dinucleotide for each library
#'
#' @examples
#'
#' calc_dinuc_percent(filenames, obs_freqs)
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
  res <- left_join(table, all_reads, by = c("cell", "virus", "time")) %>%
    mutate(dinuc_percent_reads = (count/total_umi_reads)*100)

  res

}



