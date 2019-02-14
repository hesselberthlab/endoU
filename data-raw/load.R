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
