####This should be run in workstation
library(tidymass)
library(tidyverse)

# setwd(r4projects::get_project_wd())

rm(list = ls())

###RPLC pos
load(
  "D:/shenxt/project/metabolomics_aging/monkey_data/metabolomics/R_Macaque/mzxml/RPLC/pos/Result/object"
)

object_pos <- object

##remove some QC samples
object_pos <-
  object_pos %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(-contains("QC_"))

setwd(r4projects::get_project_wd())
setwd("3-data_analysis/metabolomics_data_cleaning/RPLC/POS/")

###data quality assessment
object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class != "Blank") %>% 
  massqc::massqc_report(path = "before_data_cleaning")

qc_id <-
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

subject_id <-
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Subject") %>%
  pull(sample_id)

object_pos <-
  object_pos %>%
  mutate_variable_na_freq(according_to_samples = qc_id) %>%
  mutate_variable_na_freq(according_to_samples = subject_id)

##remove the feature that have more than 20% NA in qc or more than 50% NA in subject
object_pos <-
  object_pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(na_freq < 0.2 & na_freq.1 < 0.5)

dim(object_pos)

outlier_samples <-
  object_pos %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

outlier_table <-
  extract_outlier_table(outlier_samples)

outlier_table %>%
  dplyr::filter(according_to_na)


##for blank samples, use the 0
blank_id <- 
  object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::filter(class == 'Blank') %>% 
  pull(sample_id)

object_pos <-
  impute_mv(object = object_pos,
            sample_id = blank_id,
            method = "zero")

object_pos <-
  impute_mv(object = object_pos,
            sample_id = c(qc_id, subject_id),
            method = "knn")

object_pos <-
  normalize_data(object_pos, method = "median")

object_pos2 <-
  integrate_data(object_pos, method = "subject_median")

object_pos2

###data quality
object_pos2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class != "Blank") %>%
  massqc::massqc_report(path = "after_data_cleaning")


save(object_pos2, file = "object_pos2")




###RPLC neg
load(
  "D:/shenxt/project/metabolomics_aging/monkey_data/metabolomics/R_Macaque/mzxml/RPLC/neg/Result/object"
)

object_neg <- object

object_neg <- 
object_neg %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::mutate(group = case_when(
    stringr::str_detect(sample_id, "QC") ~ "QC",
    TRUE ~ group
  )) %>% 
  dplyr::mutate(class  = case_when(
    stringr::str_detect(sample_id, "QC") ~ "QC",
    TRUE ~ class 
  ))

object_neg

##remove some QC samples
object_neg <-
  object_neg %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(-contains("QC_"))

setwd(r4projects::get_project_wd())
setwd("3-data_analysis/metabolomics_data_cleaning/RPLC/NEG/")

###data quality assessment
object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class != "Blank") %>% 
  massqc::massqc_report(path = "before_data_cleaning")

qc_id <-
  object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

subject_id <-
  object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Subject") %>%
  pull(sample_id)

object_neg <-
  object_neg %>%
  mutate_variable_na_freq(according_to_samples = qc_id) %>%
  mutate_variable_na_freq(according_to_samples = subject_id)

##remove the feature that have more than 20% NA in qc or more than 50% NA in subject
object_neg <-
  object_neg %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(na_freq < 0.2 & na_freq.1 < 0.5)

dim(object_neg)

outlier_samples <-
  object_neg %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

outlier_table <-
  extract_outlier_table(outlier_samples)

outlier_table %>%
  dplyr::filter(according_to_na)

object_neg <-
  object_neg %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::filter(sample_id != "QC6")

##for blank samples, use the 0
blank_id <- 
  object_neg %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::filter(class == 'Blank') %>% 
  pull(sample_id)

blank_id

qc_id <-
  object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

subject_id <-
  object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Subject") %>%
  pull(sample_id)

object_neg <-
  impute_mv(object = object_neg,
            sample_id = blank_id,
            method = "zero")

object_neg <-
  impute_mv(object = object_neg,
            sample_id = c(qc_id, subject_id),
            method = "knn")

sum(is.na(object_neg))

object_neg <-
  normalize_data(object_neg, method = "median")

object_neg2 <-
  integrate_data(object_neg, method = "subject_median")

object_neg2

###data quality
object_neg2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class != "Blank") %>%
  massqc::massqc_report(path = "after_data_cleaning")


save(object_neg2, file = "object_neg2")









####This should be run in workstation
library(tidymass)
library(tidyverse)
rm(list = ls())

###HILIC pos
load(
  "D:/shenxt/project/metabolomics_aging/monkey_data/metabolomics/R_Macaque/mzxml/HILIC/pos/Result/object"
)

object_pos <- object
dim(object_pos)

##remove some QC samples
object_pos <-
  object_pos %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(-contains("QC_"))

setwd(r4projects::get_project_wd())
setwd("3-data_analysis/metabolomics_data_cleaning/HILIC/POS/")

# ###data quality assessment
# object_pos %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class != "Blank") %>% 
#   massqc::massqc_report(path = "before_data_cleaning")

qc_id <-
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

subject_id <-
  object_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Subject") %>%
  pull(sample_id)

object_pos <-
  object_pos %>%
  mutate_variable_na_freq(according_to_samples = qc_id) %>%
  mutate_variable_na_freq(according_to_samples = subject_id)

##remove the feature that have more than 20% NA in qc or more than 50% NA in subject
object_pos <-
  object_pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(na_freq < 0.2 & na_freq.1 < 0.5)

dim(object_pos)

outlier_samples <-
  object_pos %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

outlier_table <-
  extract_outlier_table(outlier_samples)

outlier_table %>%
  dplyr::filter(according_to_na)

##for blank samples, use the 0
blank_id <- 
  object_pos %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::filter(class == 'Blank') %>% 
  pull(sample_id)

sum(is.na(object_pos))

object_pos <-
  impute_mv(object = object_pos,
            sample_id = blank_id,
            method = "zero")
sum(is.na(object_pos))

object_pos <-
  impute_mv(object = object_pos,
            sample_id = c(qc_id, subject_id),
            method = "knn")
sum(is.na(object_pos))

object_pos <-
  normalize_data(object_pos, method = "median")

object_pos2 <-
  integrate_data(object_pos, method = "subject_median")

object_pos2

###data quality
object_pos2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class != "Blank") %>%
  massqc::massqc_report(path = "after_data_cleaning")


save(object_pos2, file = "object_pos2")




###HILIC neg
load(
  "D:/shenxt/project/metabolomics_aging/monkey_data/metabolomics/R_Macaque/mzxml/HILIC/neg/Result/object"
)

object_neg <- object

object_neg <- 
  object_neg %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::mutate(group = case_when(
    stringr::str_detect(sample_id, "QC") ~ "QC",
    TRUE ~ group
  )) %>% 
  dplyr::mutate(class  = case_when(
    stringr::str_detect(sample_id, "QC") ~ "QC",
    TRUE ~ class 
  ))

object_neg

##remove some QC samples
object_neg <-
  object_neg %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(-contains("QC_"))

setwd(r4projects::get_project_wd())
setwd("3-data_analysis/metabolomics_data_cleaning/HILIC/NEG/")

###data quality assessment
object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class != "Blank") %>% 
  massqc::massqc_report(path = "before_data_cleaning")

qc_id <-
  object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

subject_id <-
  object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Subject") %>%
  pull(sample_id)

object_neg <-
  object_neg %>%
  mutate_variable_na_freq(according_to_samples = qc_id) %>%
  mutate_variable_na_freq(according_to_samples = subject_id)

##remove the feature that have more than 20% NA in qc or more than 50% NA in subject
object_neg <-
  object_neg %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(na_freq < 0.2 & na_freq.1 < 0.5)

dim(object_neg)

outlier_samples <-
  object_neg %>%
  `+`(1) %>%
  log(2) %>%
  scale() %>%
  detect_outlier()

outlier_samples

outlier_table <-
  extract_outlier_table(outlier_samples)

outlier_table %>%
  dplyr::filter(according_to_na)

object_neg <-
  object_neg %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::filter(!sample_id %in% c("QC2", "QC2_180804083128"))

##for blank samples, use the 0
blank_id <- 
  object_neg %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::filter(class == 'Blank') %>% 
  pull(sample_id)

blank_id

qc_id <-
  object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

subject_id <-
  object_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(group == "Subject") %>%
  pull(sample_id)

object_neg <-
  impute_mv(object = object_neg,
            sample_id = blank_id,
            method = "zero")

object_neg <-
  impute_mv(object = object_neg,
            sample_id = c(qc_id, subject_id),
            method = "knn")

sum(is.na(object_neg))

object_neg <-
  normalize_data(object_neg, method = "median")

object_neg2 <-
  integrate_data(object_neg, method = "subject_median")

object_neg2

###data quality
object_neg2 %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class != "Blank") %>%
  massqc::massqc_report(path = "after_data_cleaning")


save(object_neg2, file = "object_neg2")



