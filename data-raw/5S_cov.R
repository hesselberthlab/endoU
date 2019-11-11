library(endoU)
library(tidyverse)
library(fs)
library(usethis)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"

r5S_tabs <- fs::path_join(c(basedir, "5S/new"))
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_1"))
dinuc_counts <- fs::path_join(c(basedir, "dinuc_byposition/r5s.dinuc.txt"))
dinuc_counts1 <- fs::path_join(c(basedir, "dinuc_byposition/r5s.dinuc_+1.txt"))

filenames <- list_names(norm_tabs)

library_sizes <- coverage_table(filenames) %>%
  dplyr::select(count, cell, virus, time) %>%
  dplyr::group_by(cell, virus, time) %>%
  dplyr::mutate(total_umi_reads = sum(count)) %>%
  dplyr::select(-count) %>%
  unique()

r5S_cov_tbl <- coverage_table(list_names(r5S_tabs))
r5S_cov_tbl <- complete_table(r5S_cov_tbl, dinuc_counts)
r5S_cov_tbl <- normalize_table(library_sizes, r5S_cov_tbl)

r5S_cov_tbl1 <- coverage_table(list_names(r5S_tabs))
r5S_cov_tbl1 <- complete_table(r5S_cov_tbl1, dinuc_counts1)
r5S_cov_tbl1 <- normalize_table(library_sizes, r5S_cov_tbl1)

usethis::use_data(r5S_cov_tbl, compress = "xz", overwrite = TRUE)
usethis::use_data(r5S_cov_tbl1, compress = "xz", overwrite = TRUE)
