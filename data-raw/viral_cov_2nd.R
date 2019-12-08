library(endoU)
library(tidyverse)
library(fs)
library(usethis)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_2"))
viral_tabs <- fs::path_join(c(basedir, "viral_bg/neg/exp2"))
dinuc_positions <- fs::path_join(c(basedir, "dinuc_byposition/mhv.dinuc.antisense.txt"))
dinuc_positions_1 <- fs::path_join(c(basedir, "dinuc_byposition/mhv.dinuc.antisense_+1.txt"))
dinuc_positions_2 <- fs::path_join(c(basedir, "dinuc_byposition/mhv.dinuc.antisense_+2.txt"))
dinuc_positions_3 <- fs::path_join(c(basedir, "dinuc_byposition/mhv.dinuc.antisense_+3.txt"))
single_nuc_positions <- fs::path_join(c(basedir, "dinuc_byposition/mhv.nuc.antisense.txt"))


filenames <- list_names(norm_tabs)

library_sizes_exp2 <- coverage_table(filenames) %>%
  dplyr::select(count, cell, virus, time) %>%
  dplyr::group_by(cell, virus, time) %>%
  dplyr::mutate(total_umi_reads = sum(count)) %>%
  dplyr::select(-count) %>%
  unique()

viral_cov_tbl_exp2 <- coverage_table(list_names(viral_tabs))
viral_cov_tbl_exp2 <- complete_table(viral_cov_tbl_exp2, dinuc_positions)
viral_cov_tbl_exp2 <- normalize_table(library_sizes_exp2, viral_cov_tbl_exp2)

viral_cov_tbl1_exp2 <- coverage_table(list_names(viral_tabs))
viral_cov_tbl1_exp2 <- complete_table(viral_cov_tbl1_exp2, dinuc_positions_1)
viral_cov_tbl1_exp2 <- normalize_table(library_sizes_exp2, viral_cov_tbl1_exp2)

viral_cov_tbl2_exp2 <- coverage_table(list_names(viral_tabs))
viral_cov_tbl2_exp2 <- complete_table(viral_cov_tbl2_exp2, dinuc_positions_2)
viral_cov_tbl2_exp2 <- normalize_table(library_sizes_exp2, viral_cov_tbl2_exp2)

viral_cov_tbl3_exp2 <- coverage_table(list_names(viral_tabs))
viral_cov_tbl3_exp2 <- complete_table(viral_cov_tbl3_exp2, dinuc_positions_3)
viral_cov_tbl3_exp2 <- normalize_table(library_sizes_exp2, viral_cov_tbl3_exp2)

viral_cov_nuc_exp2 <- coverage_table(list_names(viral_tabs))
viral_cov_nuc_exp2 <- complete_table(viral_cov_nuc_exp2, single_nuc_positions)
viral_cov_nuc_exp2 <- normalize_table(library_sizes_exp2, viral_cov_nuc_exp2)

usethis::use_data(viral_cov_tbl_exp2, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_tbl1_exp2, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_tbl2_exp2, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_tbl3_exp2, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_nuc_exp2, compress = "xz", overwrite = TRUE)

viral_cov_freq_exp2 <- coverage_table(list_names(viral_tabs))
viral_cov_freq_exp2 <- endoU::only_dinuc_complete(viral_cov_freq_exp2, dinuc_positions)
viral_cov_freq_exp2 <- normalize_table(library_sizes_exp2, viral_cov_freq_exp2)

viral_cov_freq1_exp2 <- coverage_table(list_names(viral_tabs))
viral_cov_freq1_exp2 <- endoU::only_dinuc_complete(viral_cov_freq1_exp2, dinuc_positions_1)
viral_cov_freq1_exp2 <- normalize_table(library_sizes_exp2, viral_cov_freq1_exp2)

viral_cov_freq2_exp2 <- coverage_table(list_names(viral_tabs))
viral_cov_freq2_exp2 <- endoU::only_dinuc_complete(viral_cov_freq2_exp2, dinuc_positions_2)
viral_cov_freq2_exp2 <- normalize_table(library_sizes_exp2, viral_cov_freq2_exp2)

viral_cov_freq3_exp2 <- coverage_table(list_names(viral_tabs))
viral_cov_freq3_exp2 <- endoU::only_dinuc_complete(viral_cov_freq3_exp2, dinuc_positions_3)
viral_cov_freq3_exp2 <- normalize_table(library_sizes_exp2, viral_cov_freq3_exp2)


usethis::use_data(viral_cov_freq_exp2, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_freq1_exp2, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_freq2_exp2, compress = "xz", overwrite = TRUE)
usethis::use_data(viral_cov_freq3_exp2, compress = "xz", overwrite = TRUE)
