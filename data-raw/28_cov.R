library(endoU)
library(tidyverse)
library(fs)
library(usethis)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"

r28S_tabs <- fs::path_join(c(basedir, "28S/exp1/new"))
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_1"))
dinuc_counts <- fs::path_join(c(basedir, "dinuc_byposition/r28s.dinuc.txt"))

r28S_cov_tbl <- coverage_table(list_names(r28S_tabs))
r28S_cov_tbl <- complete_table(r28S_cov_tbl, dinuc_counts)
r28S_cov_tbl <- normalize_table(norm_tabs, r28S_cov_tbl)

usethis::use_data(r28S_cov_tbl, compress = "xz", overwrite = TRUE)
