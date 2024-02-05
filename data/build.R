#!/usr/bin/env Rscript
library(targets)

Sys.setenv(TAR_PROJECT = "expt1")
tar_make()

Sys.setenv(TAR_PROJECT = "expt2")
tar_make()

Sys.setenv(TAR_PROJECT = "expt3")
tar_make()
