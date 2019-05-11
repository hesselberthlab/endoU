library(endoU)
library(tidyverse)
library(fs)
library(usethis)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_1"))
viral_tabs <- fs::path_join(c(basedir, "viral_bg/neg"))
dinuc_positions <- fs::path_join(c(basedir, "dinuc_byposition/mhv.dinuc.antisense.txt"))
dinuc_positions_1 <- fs::path_join(c(basedir, "dinuc_byposition/mhv.dinuc.antisense_+1.txt"))
dinuc_positions_2 <- fs::path_join(c(basedir, "dinuc_byposition/mhv.dinuc.antisense_+2.txt"))
dinuc_positions_3 <- fs::path_join(c(basedir, "dinuc_byposition/mhv.dinuc.antisense_+3.txt"))

filenames <- list_names(norm_tabs)

library_sizes <- coverage_table(filenames) %>%
  dplyr::select(count, cell, virus, time) %>%
  dplyr::group_by(cell, virus, time) %>%
  dplyr::mutate(total_umi_reads = sum(count)) %>%
  dplyr::select(-count) %>%
  unique()

viral_cov_tbl <- coverage_table(list_names(viral_tabs))
viral_cov_tbl <- complete_table(viral_cov_tbl, dinuc_positions)
viral_cov_tbl <- normalize_table(library_sizes, viral_cov_tbl)

viral_cov_tbl1 <- coverage_table(list_names(viral_tabs))
viral_cov_tbl1 <- complete_table(viral_cov_tbl1, dinuc_positions_1)
viral_cov_tbl1 <- normalize_table(library_sizes, viral_cov_tbl1)

viral_cov_tbl2 <- coverage_table(list_names(viral_tabs))
viral_cov_tbl2 <- complete_table(viral_cov_tbl2, dinuc_positions_2)
viral_cov_tbl2 <- normalize_table(library_sizes, viral_cov_tbl2)

viral_cov_tbl3 <- coverage_table(list_names(viral_tabs))
viral_cov_tbl3 <- complete_table(viral_cov_tbl3, dinuc_positions_3)
viral_cov_tbl3 <- normalize_table(library_sizes, viral_cov_tbl3)


usethis::use_data(viral_cov_tbl, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_tbl1, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_tbl2, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_tbl3, compress = "xz", overwrite = TRUE)


set.seed(42)
viral_cov_small <- dplyr::sample_n(viral_cov_tbl, size = 10000) %>%
  dplyr::arrange(cell, virus, time, pos)

usethis::use_data(viral_cov_small, compress = "xz", overwrite = TRUE)

viral_cov_freq <- coverage_table(list_names(viral_tabs))
viral_cov_freq <- endoU::only_dinuc_complete(viral_cov_freq, dinuc_positions)
viral_cov_freq <- normalize_table(library_sizes, viral_cov_freq)

viral_cov_freq1 <- coverage_table(list_names(viral_tabs))
viral_cov_freq1 <- endoU::only_dinuc_complete(viral_cov_freq1, dinuc_positions_1)
viral_cov_freq1 <- normalize_table(library_sizes, viral_cov_freq1)

viral_cov_freq2 <- coverage_table(list_names(viral_tabs))
viral_cov_freq2 <- endoU::only_dinuc_complete(viral_cov_freq2, dinuc_positions_2)
viral_cov_freq2 <- normalize_table(library_sizes, viral_cov_freq2)

viral_cov_freq3 <- coverage_table(list_names(viral_tabs))
viral_cov_freq3 <- endoU::only_dinuc_complete(viral_cov_freq3, dinuc_positions_3)
viral_cov_freq3 <- normalize_table(library_sizes, viral_cov_freq3)


usethis::use_data(viral_cov_freq, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_freq1, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_freq2, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_freq3, compress = "xz", overwrite = TRUE)
