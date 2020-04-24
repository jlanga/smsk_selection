#!/usr/bin/env Rscript
suppressMessages(library("optparse", quietly = T))
suppressMessages(library("tidyverse", quietly = T))


# library('optparse')

option_list <- list(
  make_option(c("-e", "--ete3"), help="Results from ete3", dest="ete3"),
  make_option(c("-f", "--fastcodeml"), help = "Results from fastcodeml", dest="fastcodeml"),
  make_option(c("-o", "--output"), help = "TSV with the adjusted p-values results", dest="output"),
  make_option(c("-p", "--pvalue"), help = "p-value threshold", default = 0.05, dest="p_threshold")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(OptionParser(option_list=option_list))


if(!all(c("ete3", "fastcodeml", "output") %in% names(arguments))) {
  stop(sprintf(parser@usage))
}


ete3 <- read_tsv(arguments$ete3, col_types = cols()) %>%
  select(orthogroup, omega_zero, pvalue) %>%
  mutate(program = "ete3")

fastcodeml <- read_tsv(arguments$fastcodeml, col_types = cols()) %>%
  select(orthogroup, omega_zero, pvalue) %>%
  mutate(program = "fastcodeml")

pvalues <- bind_rows(ete3, fastcodeml) %>%
  group_by(orthogroup) %>%
  mutate(
    pvalue_adjusted = p.adjust(p = pvalue, method = "BH"),
    method = "BH"
  ) %>%
  mutate(
    selected = all(pvalue_adjusted <= arguments$p_threshold) && n() == 6) %>%
  write_tsv(path = arguments$output)
