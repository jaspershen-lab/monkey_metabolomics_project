no_source()
####data preparation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/fMSEA/group_features.R")

#####load database
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/hmdb_ms1.rda"
)

setwd(r4projects::get_project_wd())

####RPLC pos
load(
  "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/RPLC/POS/annotation_table_rplc_pos"
)
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/RPLC/POS/expression_data_rplc_pos"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/RPLC/POS/sample_info_rplc_pos"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/RPLC/POS/variable_info_rplc_pos"
# )

####RPLC neg
load(
  "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/RPLC/NEG/annotation_table_rplc_neg"
)
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/RPLC/NEG/expression_data_rplc_neg"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/RPLC/NEG/sample_info_rplc_neg"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/RPLC/NEG/variable_info_rplc_neg"
# )

####HILIC pos
load(
  "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/HILIC/POS/annotation_table_hilic_pos"
)
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/HILIC/POS/expression_data_hilic_pos"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/HILIC/POS/sample_info_hilic_pos"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/HILIC/POS/variable_info_hilic_pos"
# )

####HILIC neg
load(
  "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/HILIC/NEG/annotation_table_hilic_neg"
)
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/HILIC/NEG/expression_data_hilic_neg"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/HILIC/NEG/sample_info_hilic_neg"
# )
# load(
#   "3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation/HILIC/NEG/variable_info_hilic_neg"
# )

dir.create("3-data_analysis/monkey_fMSEA_analysis/4_group_features")
setwd("3-data_analysis/monkey_fMSEA_analysis/4_group_features")

annotation_table_rplc_pos <-
  annotation_table_rplc_pos %>%
  dplyr::mutate(mode = "RPLC",
                polarity = "positive")

annotation_table_rplc_neg <-
  annotation_table_rplc_neg %>%
  dplyr::mutate(mode = "RPLC",
                polarity = "negative")

annotation_table_hilic_pos <-
  annotation_table_hilic_pos %>%
  dplyr::mutate(mode = "HILIC",
                polarity = "positive")

annotation_table_hilic_neg <-
  annotation_table_hilic_neg %>%
  dplyr::mutate(mode = "HILIC",
                polarity = "negative")

annotation_table_rplc_pos <-
  annotation_table_rplc_pos %>%
  dplyr::mutate(variable_id_old = variable_id) %>%
  dplyr::mutate(variable_id = paste(variable_id, mode, polarity, sep = "_"))

annotation_table_rplc_neg <-
  annotation_table_rplc_neg %>%
  dplyr::mutate(variable_id_old = variable_id) %>%
  dplyr::mutate(variable_id = paste(variable_id, mode, polarity, sep = "_"))

annotation_table_hilic_pos <-
  annotation_table_hilic_pos %>%
  dplyr::mutate(variable_id_old = variable_id) %>%
  dplyr::mutate(variable_id = paste(variable_id, mode, polarity, sep = "_"))

annotation_table_hilic_neg <-
  annotation_table_hilic_neg %>%
  dplyr::mutate(variable_id_old = variable_id) %>%
  dplyr::mutate(variable_id = paste(variable_id, mode, polarity, sep = "_"))

annotation_table <-
  rbind(
    annotation_table_rplc_pos,
    annotation_table_rplc_neg,
    annotation_table_hilic_pos,
    annotation_table_hilic_neg
  )

dim(annotation_table)

annotation_table <-
  group_features(annotation_table = annotation_table,
                 rt.match.tol = 5)

save(annotation_table, file = "annotation_table", compress = "xz")

# save(expression_data_rplc_pos, file = "expression_data_rplc_pos")
# save(sample_info_rplc_pos, file = "sample_info_rplc_pos")
# save(variable_info_rplc_pos, file = "variable_info_rplc_pos")
# 
# save(expression_data_rplc_neg, file = "expression_data_rplc_neg")
# save(sample_info_rplc_neg, file = "sample_info_rplc_neg")
# save(variable_info_rplc_neg, file = "variable_info_rplc_neg")
# 
# save(expression_data_hilic_pos, file = "expression_data_hilic_pos")
# save(sample_info_hilic_pos, file = "sample_info_hilic_pos")
# save(variable_info_hilic_pos, file = "variable_info_hilic_pos")
# 
# save(expression_data_hilic_neg, file = "expression_data_hilic_neg")
# save(sample_info_hilic_neg, file = "sample_info_hilic_neg")
# save(variable_info_hilic_neg, file = "variable_info_hilic_neg")

