no_source()
####data preparation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/fMSEA/remove_redundancy.R")
source("1-code/fMSEA/calculate_redundancy.R")
source("1-code/tools.R")

load("3-data_analysis/monkey_fMSEA_analysis/5_redundant_removal/annotation_table")

load(
  "3-data_analysis/monkey_fMSEA_analysis/1_data_preparation/RPLC/POS/object_rplc_pos"
)
load(
  "3-data_analysis/monkey_fMSEA_analysis/1_data_preparation/RPLC/NEG/object_rplc_neg"
)
load(
  "3-data_analysis/monkey_fMSEA_analysis/1_data_preparation/HILIC/POS/object_hilic_pos"
)
load(
  "3-data_analysis/monkey_fMSEA_analysis/1_data_preparation/HILIC/NEG/object_hilic_neg"
)

dir.create("3-data_analysis/monkey_fMSEA_analysis/7_metabolite_set_enrichment_analsyis")
setwd("3-data_analysis/monkey_fMSEA_analysis/7_metabolite_set_enrichment_analsyis")


######calculate the correlation between feature and ages
##RPLC pos
temp_object <-
  object_rplc_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject" & !is.na(age)) %>%
  `+`(1) %>%
  log(2) %>%
  scale()

cor_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    temp_data$sex[temp_data$sex == "M"] <- 1
    temp_data$sex[temp_data$sex == "F"] <- 0
    temp_data$sex <- as.numeric(temp_data$sex)
    result <-
      lm(formula = value ~ sex, data = temp_data)
    temp_data$value <- result$residuals
    
    cor_result <-
      cor.test(temp_data$value, temp_data$age, method = "spearman")
    data.frame(cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

object_rplc_pos <-
  object_rplc_pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(cor_p = cor_data$cor_p,
                spearman_cor = cor_data$spearman_cor)



##RPLC neg
temp_object <-
  object_rplc_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject" & !is.na(age)) %>%
  `+`(1) %>%
  log(2) %>%
  scale()

cor_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    temp_data$sex[temp_data$sex == "M"] <- 1
    temp_data$sex[temp_data$sex == "F"] <- 0
    temp_data$sex <- as.numeric(temp_data$sex)
    result <-
      lm(formula = value ~ sex, data = temp_data)
    temp_data$value <- result$residuals
    
    cor_result <-
      cor.test(temp_data$value, temp_data$age, method = "spearman")
    data.frame(cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

object_rplc_neg <-
  object_rplc_neg %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(cor_p = cor_data$cor_p,
                spearman_cor = cor_data$spearman_cor)

##HILIC pos
temp_object <-
  object_hilic_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject" & !is.na(age)) %>%
  `+`(1) %>%
  log(2) %>%
  scale()

cor_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    temp_data$sex[temp_data$sex == "M"] <- 1
    temp_data$sex[temp_data$sex == "F"] <- 0
    temp_data$sex <- as.numeric(temp_data$sex)
    result <-
      lm(formula = value ~ sex, data = temp_data)
    temp_data$value <- result$residuals
    
    cor_result <-
      cor.test(temp_data$value, temp_data$age, method = "spearman")
    data.frame(cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

object_hilic_pos <-
  object_hilic_pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(cor_p = cor_data$cor_p,
                spearman_cor = cor_data$spearman_cor)


##HILIC neg
temp_object <-
  object_hilic_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject" & !is.na(age)) %>%
  `+`(1) %>%
  log(2) %>%
  scale()

cor_data <-
  seq_len(nrow(temp_object)) %>%
  purrr::map(function(i) {
    cat(i, " ")
    value <-
      as.numeric(unlist(temp_object[i, , drop = TRUE]))
    sample_info <-
      temp_object@sample_info
    temp_data <-
      data.frame(sample_info, value)
    temp_data$sex[temp_data$sex == "M"] <- 1
    temp_data$sex[temp_data$sex == "F"] <- 0
    temp_data$sex <- as.numeric(temp_data$sex)
    result <-
      lm(formula = value ~ sex, data = temp_data)
    temp_data$value <- result$residuals
    
    cor_result <-
      cor.test(temp_data$value, temp_data$age, method = "spearman")
    data.frame(cor_p = cor_result$p.value,
               spearman_cor = cor_result$estimate)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

object_hilic_neg <-
  object_hilic_neg %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(cor_p = cor_data$cor_p,
                spearman_cor = cor_data$spearman_cor)

extract_variable_info(object_rplc_neg) %>%
  dplyr::select(variable_id, spearman_cor)
