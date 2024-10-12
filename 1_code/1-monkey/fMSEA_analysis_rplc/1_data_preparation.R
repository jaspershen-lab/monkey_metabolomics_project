no_source()

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())

#####RPLC positive
load("3-data_analysis/monkey_metabolomics_data_cleaning/RPLC/POS/object_pos2")
object_rplc_pos <- object_pos2
load("3-data_analysis/monkey_metabolomics_data_cleaning/RPLC/NEG/object_neg2")
object_rplc_neg <- object_neg2

load("3-data_analysis/monkey_metabolomics_data_cleaning/HILIC/POS/object_pos2")
object_hilic_pos <- object_pos2
load("3-data_analysis/monkey_metabolomics_data_cleaning/HILIC/NEG/object_neg2")
object_hilic_neg <- object_neg2

dir.create("3-data_analysis/monkey_fMSEA_analysis_rplc/1_data_preparation/", recursive = TRUE)
setwd("3-data_analysis/monkey_fMSEA_analysis_rplc/1_data_preparation/")

##remove redundant samples
###RPLC positive
colnames(object_rplc_pos)

object_rplc_pos <-
  object_rplc_pos %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(-contains("blkA")) %>%
  dplyr::select(-contains("blkB"))

# ##QC8
# plot(log(object_rplc_pos$QC8, 2),
#      log(object_rplc_pos$QC8_180818105226, 2))
#
# ##17O
# plot(log(object_rplc_pos$`17O`, 2),
#      log(object_rplc_pos$`17O_180818114837`, 2))
#
# ##47Z
# plot(log(object_rplc_pos$`47Z`, 2),
#      log(object_rplc_pos$`47Z_180818111119`, 2))
#
# ##47Z
# plot(log(object_rplc_pos$`65B`, 2),
#      log(object_rplc_pos$`65B_180818112958`, 2))

object_rplc_pos <-
  object_rplc_pos %>%
  dplyr::select(-c(
    QC8_180818105226,
    `17O_180818114837`,
    `47Z_180818111119`,
    `65B_180818112958`
  )) %>%
  dplyr::rename(`44Y` = `44Y_180818120716`)

object_rplc_pos <-
  object_rplc_pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(mode = "RPLC_positive") %>%
  dplyr::mutate(column = "RPLC",
                polarity = "positive")

###RPLC neg
sort(colnames(object_rplc_neg))

object_rplc_neg <-
  object_rplc_neg %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(-blk1_180816162622) %>%
  dplyr::rename(QC6 = `QC6_180818103341`)

# ##20F
# plot(log(object_rplc_neg$`20F_180819095042`, 2),
#      log(object_rplc_neg$`20F`, 2))

object_rplc_neg <-
  object_rplc_neg %>%
  dplyr::select(-c(`20F_180819095042`))

object_rplc_neg <-
  object_rplc_neg %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(mode = "RPLC_negative") %>%
  dplyr::mutate(column = "RPLC",
                polarity = "negative")

# ###HILIC positive
# sort(colnames(object_hilic_pos))
# 
# object_hilic_pos <-
#   object_hilic_pos %>%
#   activate_mass_dataset(what = "expression_data") %>%
#   dplyr::select(-contains("blkA")) %>%
#   dplyr::select(-contains("blkB"))
# 
# 
# object_hilic_pos <-
#   object_hilic_pos %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::mutate(mode = "HILIC_positive") %>%
#   dplyr::mutate(column = "HILIC",
#                 polarity = "positive")
# 
# ###HILIC neg
# sort(colnames(object_hilic_neg))
# 
# # ##13Q
# # plot(log(object_hilic_neg$`13Q`, 2),
# #      log(object_hilic_neg$`13Q_180804135213`, 2))
# #
# # ##54V
# # plot(log(object_hilic_neg$`54V`, 2),
# #      log(object_hilic_neg$`54V_180804095144`, 2))
# #
# # ##54V
# # plot(log(object_hilic_neg$`66F`, 2),
# #      log(object_hilic_neg$`66F_180804092500`, 2))
# #
# # ##QC2
# # plot(
# #   log(object_hilic_neg$`QC2_180804125843`, 2),
# #   log(object_hilic_neg$`QC2_180804132529`, 2)
# # )
# 
# object_hilic_neg <-
#   object_hilic_neg %>%
#   activate_mass_dataset(what = "expression_data") %>%
#   dplyr::select(
#     -c(
#       `13Q_180804135213`,
#       `54V_180804095144`,
#       `66F_180804092500`,
#       `QC2_180804132529`
#     )
#   ) %>%
#   dplyr::rename(`QC2` = `QC2_180804125843`)
# 
# object_hilic_neg <-
#   object_hilic_neg %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::mutate(mode = "HILIC_negative") %>%
#   dplyr::mutate(column = "HILIC",
#                 polarity = "negative")
# 
# ###rename variable id
# object_rplc_pos <-
#   object_rplc_pos %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::mutate(variable_id = paste(variable_id, "RPLC_pos", sep = "_"))
# 
# object_rplc_neg <-
#   object_rplc_neg %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::mutate(variable_id = paste(variable_id, "RPLC_neg", sep = "_"))
# 
# object_hilic_pos <-
#   object_hilic_pos %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::mutate(variable_id = paste(variable_id, "HILIC_pos", sep = "_"))
# 
# object_hilic_neg <-
#   object_hilic_neg %>%
#   activate_mass_dataset(what = "variable_info") %>%
#   dplyr::mutate(variable_id = paste(variable_id, "HILIC_neg", sep = "_"))
# 
# object_hilic_pos <-
#   object_hilic_pos %>%
#   activate_mass_dataset(what = "expression_data") %>%
#   dplyr::select(-blk4)
# 
# object_hilic_neg <-
#   object_hilic_neg %>%
#   activate_mass_dataset(what = "expression_data") %>%
#   dplyr::select(-blk4)

intersect(colnames(object_rplc_pos), colnames(object_rplc_neg))
# intersect(colnames(object_hilic_pos), colnames(object_hilic_neg))
# intersect(colnames(object_rplc_pos), colnames(object_hilic_pos))
# setdiff(colnames(object_hilic_pos), colnames(object_rplc_pos))

object_rplc_pos <-
  object_rplc_pos[, colnames(object_rplc_pos)]

object_rplc_neg <-
  object_rplc_neg[, colnames(object_rplc_pos)]

# object_hilic_pos <-
#   object_hilic_pos[, colnames(object_rplc_pos)]
# 
# object_hilic_neg <-
#   object_hilic_neg[, colnames(object_rplc_pos)]

metadata <-
  readxl::read_xlsx(
    here::here(
      "3-data_analysis/monkey_metabolomics_data_preparation/metadata.xlsx"
    ),
    col_types = c("text", "text", "date", "text", "text", "text", "text")
  ) %>%
  dplyr::select(`ANIMAL ID`:SIRE)

metadata <-
  metadata %>%
  dplyr::rename(
    subject_id = `ANIMAL ID`,
    sample_id = `ANIMAL ID`,
    dob = DOB,
    sex = SEX,
    mother = DAME,
    Father = SIRE
  ) %>%
  dplyr::mutate(age = as.numeric(as.Date("2018-7-1") - as.Date(dob)) / 365)

object_rplc_pos <-
  object_rplc_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::left_join(metadata, by = "sample_id")

object_rplc_neg <-
  object_rplc_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::left_join(metadata, by = "sample_id")

# object_hilic_pos <-
#   object_hilic_pos %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::left_join(metadata, by = "sample_id")
# 
# object_hilic_neg <-
#   object_hilic_neg %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::left_join(metadata, by = "sample_id")

sample_list <-
  readxl::read_xlsx(
    here::here(
      "3-data_analysis/monkey_metabolomics_data_preparation/sample_list.xlsx"
    )
  )

sample_list <-
  sample_list %>%
  dplyr::rename(sample_id = `File Name`) %>%
  dplyr::distinct(sample_id)

sample_list <-
  sample_list %>%
  dplyr::filter(sample_id %in% colnames(object_rplc_pos)) %>%
  dplyr::select(sample_id) %>%
  dplyr::mutate(injection.order = 1:nrow(.))

object_rplc_pos <-
  object_rplc_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::select(-injection.order) %>%
  dplyr::left_join(sample_list, by = "sample_id")

object_rplc_neg <-
  object_rplc_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::select(-injection.order) %>%
  dplyr::left_join(sample_list, by = "sample_id")

# object_hilic_pos <-
#   object_hilic_pos %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::select(-injection.order) %>%
#   dplyr::left_join(sample_list, by = "sample_id")
# 
# object_hilic_neg <-
#   object_hilic_neg %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::select(-injection.order) %>%
#   dplyr::left_join(sample_list, by = "sample_id")

####add MS2
object_rplc_pos <-
  object_rplc_pos %>%
  massdataset::mutate_ms2(
    column = "rp",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "RPLC/POS/"
  )

####add MS2
object_rplc_neg <-
  object_rplc_neg %>%
  massdataset::mutate_ms2(
    column = "rp",
    polarity = "negative",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "RPLC/NEG/"
  )

# ####add MS2
# object_hilic_pos <-
#   object_hilic_pos %>%
#   massdataset::mutate_ms2(
#     column = "hilic",
#     polarity = "positive",
#     ms1.ms2.match.mz.tol = 15,
#     ms1.ms2.match.rt.tol = 30,
#     path = "HILIC/POS/"
#   )
# 
# 
# ####add MS2
# object_hilic_neg <-
#   object_hilic_neg %>%
#   massdataset::mutate_ms2(
#     column = "hilic",
#     polarity = "negative",
#     ms1.ms2.match.mz.tol = 15,
#     ms1.ms2.match.rt.tol = 30,
#     path = "HILIC/NEG/"
#   )


####calculate correlation

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

sample_id_qc <-
  object_rplc_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

object_rplc_pos <-
  object_rplc_pos %>%
  mutate_mean_intensity(according_to_samples = sample_id_qc)

sample_id_qc <-
  object_rplc_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

object_rplc_neg <-
  object_rplc_neg %>%
  mutate_mean_intensity(according_to_samples = sample_id_qc)

save(object_rplc_pos, file = "RPLC/POS/object_rplc_pos")
save(object_rplc_neg, file = "RPLC/NEG/object_rplc_neg")

# save(object_hilic_pos, file = "HILIC/POS/object_hilic_pos")
# save(object_hilic_neg, file = "HILIC/NEG/object_hilic_neg")
