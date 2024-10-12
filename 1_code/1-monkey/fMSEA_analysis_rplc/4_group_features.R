no_source()
####data preparation

library(tidymass)
library(tidyverse)

rm(list = ls())
gc()

setwd(r4projects::get_project_wd())
source("1-code/fMSEA/group_feature4fmsea.R")
source("1-code/fMSEA/data_checking.R")

#####load database
setwd(r4projects::get_project_wd())

load(
  "3-data_analysis/monkey_fMSEA_analysis_rplc/3_isotope_annotation/object_list"
)

dir.create("3-data_analysis/monkey_fMSEA_analysis_rplc/4_group_features")
setwd("3-data_analysis/monkey_fMSEA_analysis_rplc/4_group_features")

annotation_table <-
  group_features4fmsea(object_list = object_list,
                       rt.match.tol = 10,
                       threads = 3)

save(annotation_table, file = "annotation_table", compress = "xz")


