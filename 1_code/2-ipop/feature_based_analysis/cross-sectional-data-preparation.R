no_function()

###
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1-code/tools.R")
library(tidyverse)
library(tidymass)

load("3-data_analysis/ipop_metabolomics_data_preparation/peak/object")

setwd("3-data_analysis/ipop_metabolomics_data_preparation/peak")

object

dim(object)

object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(adjusted_age))

####only remain the health vist
object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(CL4 == "Healthy")


###for each participant, calculate the median value
object_corss_section <-
  massdataset::summarise_samples(object = object,
                                 what = "mean_intensity",
                                 group_by = "subject_id")

save(object_corss_section, file = "object_corss_section")
