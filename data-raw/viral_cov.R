library(endoU)
library(tidyverse)
library(fs)
library(usethis)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"

viral_tabs <- fs::path_join(c(basedir, "viral_bg/neg"))
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_1"))
dinuc_counts <- fs::path_join(c(basedir, "dinuc_byposition/mhv.dinuc.neg.table"))

viral_cov_tbl <- coverage_table(list_names(viral_tabs))
viral_cov_tbl <- complete_table(viral_cov_tbl, dinuc_counts)
viral_cov_tbl <- normalize_table(norm_tabs, viral_cov_tbl)

usethis::use_data(viral_cov_tbl, compress = "xz", overwrite = TRUE)

set.seed(42)
viral_cov_small <- sample_n(viral_cov_tbl, size = 10000) %>%
  arrange(cell, virus, time, pos)

usethis::use_data(viral_cov_small, compress = "xz", overwrite = TRUE)

viral_cov_freq <- coverage_table(list_names(viral_tabs))
viral_cov_freq <- only_dinuc_complete(viral_cov_freq, dinuc_counts)
viral_cov_freq <- normalize_table(norm_tabs, viral_cov_freq)

usethis::use_data(viral_cov_freq, compress = "xz", overwrite = TRUE)
