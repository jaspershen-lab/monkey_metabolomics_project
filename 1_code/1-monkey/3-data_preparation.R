no_source()

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())

load(
  "3-data_analysis/monkey_metabolomics_metabolite_annotation/RPLC/POS/object_rplc_pos"
)

load(
  "3-data_analysis/monkey_metabolomics_metabolite_annotation/RPLC/NEG/object_rplc_neg"
)

load(
  "3-data_analysis/monkey_metabolomics_metabolite_annotation/HILIC/POS/object_hilic_pos"
)

load(
  "3-data_analysis/monkey_metabolomics_metabolite_annotation/HILIC/NEG/object_hilic_neg"
)

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
  dplyr::mutate(mode = "RPLC_positive")

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
  dplyr::mutate(mode = "RPLC_negative")

###HILIC positive
sort(colnames(object_hilic_pos))

object_hilic_pos <-
  object_hilic_pos %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(-contains("blkA")) %>%
  dplyr::select(-contains("blkB"))


object_hilic_pos <-
  object_hilic_pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(mode = "HILIC_positive")


###HILIC neg
sort(colnames(object_hilic_neg))

# ##13Q
# plot(log(object_hilic_neg$`13Q`, 2),
#      log(object_hilic_neg$`13Q_180804135213`, 2))
#
# ##54V
# plot(log(object_hilic_neg$`54V`, 2),
#      log(object_hilic_neg$`54V_180804095144`, 2))
#
# ##54V
# plot(log(object_hilic_neg$`66F`, 2),
#      log(object_hilic_neg$`66F_180804092500`, 2))
#
# ##QC2
# plot(
#   log(object_hilic_neg$`QC2_180804125843`, 2),
#   log(object_hilic_neg$`QC2_180804132529`, 2)
# )

object_hilic_neg <-
  object_hilic_neg %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(
    -c(
      `13Q_180804135213`,
      `54V_180804095144`,
      `66F_180804092500`,
      `QC2_180804132529`
    )
  ) %>%
  dplyr::rename(`QC2` = `QC2_180804125843`)

object_hilic_neg <-
  object_hilic_neg %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(mode = "HILIC_negative")

###rename variable id
object_rplc_pos <-
  object_rplc_pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste(variable_id, "RPLC_pos", sep = "_"))

object_rplc_neg <-
  object_rplc_neg %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste(variable_id, "RPLC_neg", sep = "_"))

object_hilic_pos <-
  object_hilic_pos %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste(variable_id, "HILIC_pos", sep = "_"))

object_hilic_neg <-
  object_hilic_neg %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(variable_id = paste(variable_id, "HILIC_neg", sep = "_"))

object_hilic_pos <-
  object_hilic_pos %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(-blk4)

object_hilic_neg <-
  object_hilic_neg %>%
  activate_mass_dataset(what = "expression_data") %>%
  dplyr::select(-blk4)

intersect(colnames(object_rplc_pos), colnames(object_rplc_neg))
intersect(colnames(object_hilic_pos), colnames(object_hilic_neg))
intersect(colnames(object_rplc_pos), colnames(object_hilic_pos))
setdiff(colnames(object_hilic_pos), colnames(object_rplc_pos))

object_rplc_pos <-
  object_rplc_pos[, colnames(object_rplc_pos)]

object_rplc_neg <-
  object_rplc_neg[, colnames(object_rplc_pos)]

object_hilic_pos <-
  object_hilic_pos[, colnames(object_rplc_pos)]

object_hilic_neg <-
  object_hilic_neg[, colnames(object_rplc_pos)]

setwd(r4projects::get_project_wd())

setwd("3-data_analysis/monkey_metabolomics_data_preparation/")

metadata <-
  readxl::read_xlsx("metadata.xlsx",
                    col_types = c("text", "text", "date", "text", "text", "text", "text")) %>%
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

object_hilic_pos <-
  object_hilic_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::left_join(metadata, by = "sample_id")

object_hilic_neg <-
  object_hilic_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::left_join(metadata, by = "sample_id")

sample_list <-
  readxl::read_xlsx("sample_list.xlsx")

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

object_hilic_pos <-
  object_hilic_pos %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::select(-injection.order) %>%
  dplyr::left_join(sample_list, by = "sample_id")

object_hilic_neg <-
  object_hilic_neg %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::select(-injection.order) %>%
  dplyr::left_join(sample_list, by = "sample_id")

# save(object_rplc_pos, file = "object_rplc_pos")
# save(object_rplc_neg, file = "object_rplc_neg")
# save(object_hilic_pos, file = "object_hilic_pos")
# save(object_hilic_neg, file = "object_hilic_neg")

####peaks
object <-
  rbind(object_rplc_pos,
        object_rplc_neg,
        object_hilic_pos,
        object_hilic_neg)

colnames(object@sample_info)

colnames(object)

object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::select(-contains("_2")) %>%
  dplyr::select(-contains("_3")) %>%
  dplyr::select(-contains("_4"))

# save(object, file = "peak/object")

export_mass_dataset(object, file_type = "xlsx", path = "peak/")

###metabolites
object <-
  object %>%
  activate_mass_dataset(what = "annotation_table") %>%
  dplyr::filter(!is.na(Level))

dim(object)

###remove redundant metabolites
variable_info <-
  extract_variable_info(object)

library(plyr)

variable_info <-
  variable_info %>%
  plyr::dlply(.variables = .(Compound.name))

variable_info <-
  variable_info %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    } else{
      x %>%
        dplyr::filter(Level == min(Level)) %>%
        dplyr::filter(SS == max(SS)) %>%
        dplyr::filter(Total.score == max(Total.score)) %>%
        head(1)
    }
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

object <-
  object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% variable_info$variable_id)

library(plyr)

variable_info <-
  extract_variable_info(object) %>%
  plyr::dlply(.variables = .(Compound.name)) %>%
  purrr::map(function(x) {
    if (nrow(x) == 1) {
      return(x)
    } else{
      x %>%
        dplyr::filter(Level == min(Level)) %>%
        dplyr::filter(SS == max(SS, na.rm = TRUE)) %>%
        dplyr::filter(Total.score == max(Total.score, na.rm = TRUE)) %>%
        head(1)
    }
  }) %>%
  dplyr::bind_rows()

object <-
  object %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% variable_info$variable_id)

save(object, file = "metabolite/object")

export_mass_dataset(object = object,
                    file_type = "xlsx",
                    path = "metabolite/")
