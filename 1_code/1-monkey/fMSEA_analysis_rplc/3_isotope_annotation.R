no_source()
####data preparation

library(tidymass)
library(tidyverse)

###RPLC pos
rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/fMSEA/annotate_isotopes4fmsea.R")
source("1-code/fMSEA/data_checking.R")

setwd(r4projects::get_project_wd())

load(
  "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/object_list"
)

dir.create("3-data_analysis/monkey_fMSEA_analysis_rplc/3_isotope_annotation")
setwd("3-data_analysis/monkey_fMSEA_analysis_rplc/3_isotope_annotation")

object_list <-
  annotate_isotopes4fmsea(
    object_list = object_list,
    threads = 3,
    rt_tol = 5,
    mz_tol = 15,
    intensity_tol = 0.5,
    max_isotope = 2
  )

save(object_list, file = "object_list")
