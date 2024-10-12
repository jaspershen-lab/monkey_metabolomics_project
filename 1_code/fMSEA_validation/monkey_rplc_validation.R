no_source()
####data preparation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())

source("1-code/fMSEA/annotate_isotopes4fmsea.R")
source("1-code/fMSEA/annotate_metabolite4fmsea.R")
source("1-code/fMSEA/calculate_redundancy.R")
source("1-code/fMSEA/data_checking.R")
source("1-code/fMSEA/do_fmsea.R")
source("1-code/fMSEA/group_feature4fmsea.R")
source("1-code/fMSEA/remove_redundancy4fmsea.R")

load(
  "3-data_analysis/monkey_fMSEA_analysis_rplc/1_data_preparation/RPLC/POS/object_rplc_pos"
)

object_rplc_pos

load(
  "3-data_analysis/monkey_fMSEA_analysis_rplc/1_data_preparation/RPLC/NEG/object_rplc_neg"
)

object_rplc_neg

dir.create("3-data_analysis/fMSEA_validation/monkey_rplc", recursive = TRUE)
setwd("3-data_analysis/fMSEA_validation/monkey_rplc")

object_list <-
  list(object_rplc_pos,
       object_rplc_neg)

# rm(list = c("object_rplc_pos", "object_rplc_neg"))
# gc()

###load databases
load(
  here::here(
    "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/hmdb_ms2.rda"
  )
)

load(
  here::here(
    "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/massbank_ms2.rda"
  )
)

load(
  here::here(
    "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/metlin_ms2.rda"
  )
)

load(
  here::here(
    "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/mona_ms2.rda"
  )
)

load(
  here::here(
    "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/mpsnyder_rplc_ms2.rda"
  )
)

load(
  here::here(
    "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/nist_ms2.rda"
  )
)

load(
  here::here(
    "3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/database/hmdb_ms1.rda"
  )
)

###only remain the metabolites have HMDB.ID, CAS.ID or KEGG.ID
hmdb_ms2 <-
  hmdb_ms2 %>%
  dplyr::filter(!is.na(HMDB.ID))

massbank_ms2 <-
  massbank_ms2 %>%
  dplyr::filter(!is.na(HMDB.ID))

metlin_ms2 <-
  metlin_ms2 %>%
  dplyr::filter(!is.na(HMDB.ID))

mpsnyder_rplc_ms2 <-
  mpsnyder_rplc_ms2 %>%
  dplyr::filter(!is.na(HMDB.ID))

mona_ms2 <-
  mona_ms2 %>%
  dplyr::filter(!is.na(HMDB.ID))

nist_ms2 <-
  nist_ms2 %>%
  dplyr::filter(!is.na(HMDB.ID))

hmdb_ms1 <-
  hmdb_ms1 %>%
  dplyr::filter(!is.na(HMDB.ID)) %>%
  dplyr::filter(!is.na(Formula)) %>%
  dplyr::filter(nchar(Formula) > 4) %>%
  dplyr::filter(!stringr::str_detect(Formula, "Cu|Mg|Al|Si|Ca"))


# ####Do fMSEA using the public MS2 library
# system.time(
#   do_fmsea(
#     object_list = object_list,
#     databases_list = list(
#       hmdb_ms2,
#       massbank_ms2,
#       metlin_ms2,
#       mona_ms2,
#       nist_ms2,
#       hmdb_ms1
#     ),
#     use_default_database = FALSE,
#     rt.match.tol = 30,
#     threads = 3,
#     remove_fragment_intensity_cutoff = 0.01,
#     isotope_match_rt_tol = 10,
#     isotope_match_mz_tol = 15,
#     isotope_match_intensity_tol = 0.5,
#     max_isotope = 2,
#     grouping_rt_tol = 10,
#     score_cutoff = 20,
#     candidate.num = 10000,
#     path = "level2_and_3"
#   )
# )

# ####Do fMSEA only using the MS1 HMDB database
# system.time(
#   do_fmsea(
#     object_list = object_list,
#     databases_list = list(
#       hmdb_ms1
#     ),
#     use_default_database = FALSE,
#     rt.match.tol = 30,
#     threads = 3,
#     remove_fragment_intensity_cutoff = 0.01,
#     isotope_match_rt_tol = 10,
#     isotope_match_mz_tol = 15,
#     isotope_match_intensity_tol = 0.5,
#     max_isotope = 2,
#     grouping_rt_tol = 10,
#     score_cutoff = 20,
#     path = "only_levele3"
#   )
# )

# # #####in-house library
# object_list <-
#   annotate_metabolite4fmsea(
#     object_list = object_list,
#     databases_list = list(mpsnyder_rplc_ms2),
#     use_default_database = FALSE,
#     same_lc_with_inhouse_library = TRUE,
#     rt.match.tol = 30,
#     threads = 3,
#     remove_fragment_intensity_cutoff = 0.01,
#     path = "Level1_annotation"
#   )
#
# annotation_table_pos <-
# object_list[[1]]@annotation_table %>%
#   dplyr::mutate(variable_id = paste0(variable_id, "_RPLC_positive"))
#
# annotation_table_neg <-
#   object_list[[2]]@annotation_table %>%
#   dplyr::mutate(variable_id = paste0(variable_id, "_RPLC_negative"))
#
# annotation_table_level1 <-
#   rbind(annotation_table_pos,
#         annotation_table_neg)
#
# save(annotation_table_level1, file = 'annotation_table_level1')
load("annotation_table_level1")

annotation_table_level1 <-
  annotation_table_level1 %>%
  plyr::dlply(.variables = plyr::.(variable_id)) %>%
  purrr::map(function(x) {
    x %>%
      dplyr::filter(SS == max(SS)) %>%
      head(1)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

dim(annotation_table_level1)

######annotation result using the HMDB MS1 database
load("only_levele3/redundancy_removing/annotation_table")
calculate_redundancy(annotation_table = annotation_table)

# temp <-
# annotation_table %>%
#   dplyr::filter(variable_id %in% no_matched_variable)

validation_result <-
  annotation_table_level1[, c(
    "variable_id",
    "Compound.name",
    "CAS.ID",
    "HMDB.ID",
    "KEGG.ID",
    "Lab.ID",
    "Adduct",
    "mz.error",
    "RT.error",
    "SS",
    "Total.score",
    "Database",
    "Level"
  )] %>%
  dplyr::left_join(annotation_table[, c(
    "variable_id",
    "Compound.name",
    "CAS.ID",
    "HMDB.ID",
    "KEGG.ID",
    "Lab.ID",
    "Adduct",
    "mz.error",
    "RT.error",
    "SS",
    "Total.score",
    "Database",
    "Level",
    "compound_class",
    "score",
    "Formula"
  )],
  by = "variable_id")

validation_result <-
  validation_result %>%
  dplyr::left_join(
    rbind(
      object_list[[1]]@variable_info[, c("variable_id", "mz")] %>%
        dplyr::mutate(variable_id = paste(variable_id, "RPLC_positive", sep = "_")),
      object_list[[2]]@variable_info[, c("variable_id", "mz")] %>%
        dplyr::mutate(variable_id = paste(variable_id, "RPLC_negative", sep = "_"))
    ),
    by = "variable_id"
  ) %>%
  plyr::dlply(.variables = plyr::.(variable_id)) %>%
  purrr::map(function(x) {
    x[, c(
      "variable_id",
      "mz",
      "Adduct.x",
      "Adduct.y",
      "Compound.name.x",
      "Compound.name.y",
      "HMDB.ID.x",
      "HMDB.ID.y",
      "score"
    )] %>%
      dplyr::arrange(desc(score))
  })

validation_result[[1]]

result <-
  validation_result %>%
  purrr::map(function(x) {
    idx <-
      which(x$HMDB.ID.x == x$HMDB.ID.y)
    if (length(idx) == 0) {
      temp <- x[1,]
      temp$Compound.name.y <- NA
      temp$HMDB.ID.y <- NA
      temp$score <- NA
      temp$rank <- NA
      return(temp)
    }
    x$rank <- idx
    x[idx,]
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

result %>%
  dplyr::mutate(idx = 1:nrow(result)) %>%
  dplyr::arrange(rank) %>%
  dplyr::mutate(idx = factor(idx, levels = idx)) %>%
  ggplot(aes(idx, rank)) +
  geom_point(aes(size = score)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

sum(is.na(result$rank))

result %>%
  dplyr::filter(is.na(rank))

# annotation_table_level1 %>%
#   dplyr::filter(variable_id == "M162T34_POS_RPLC_positive" &
#                   HMDB.ID == "HMDB0000062")

# x <-
# metid::annotate_single_peak_mass_dataset(object = object_rplc_neg,
#                                          variable_id = "M564T771_NEG",
#                                          based_on_rt = FALSE,
#                                          based_on_ms2 = FALSE,
#                                          add_to_annotation_table = FALSE,
#                                          polarity = "negative",
#                                          column = "rp",
#                                          candidate.num = 100,
#                                          database = hmdb_ms1)

which(x$HMDB.ID == "HMDB0000207")

# load("only_levele3/metabolite_annotation/object_list")
#
# temp_pos <-
#   temp %>%
#   dplyr::filter(stringr::str_detect(variable_id, "positive"))
#
# temp_neg <-
#   temp %>%
#   dplyr::filter(stringr::str_detect(variable_id, "negative"))
#
# stringr::str_replace(temp_pos$variable_id, "_RPLC_positive", "") %in% object_list[[1]]@annotation_table$variable_id
# stringr::str_replace(temp_neg$variable_id, "_RPLC_negative", "") %in% object_list[[2]]@annotation_table$variable_id
#
# result_pos <-
# object_list[[1]]@annotation_table %>%
#   dplyr::filter(variable_id %in% stringr::str_replace(temp_pos$variable_id, "_RPLC_positive", ""))
#
# result_neg <-
#   object_list[[2]]@annotation_table %>%
#   dplyr::filter(variable_id %in% stringr::str_replace(temp_neg$variable_id, "_RPLC_negative", ""))
#
# result_pos %>%
#   dplyr::filter(HMDB.ID %in% temp_pos$HMDB.ID.x)
#

# ######annotation result using the HMDB MS1 database and public databases
# load("level2_and_3/redundancy_removing/annotation_table")
# calculate_redundancy(annotation_table = annotation_table)
#
# validation_result <-
#   annotation_table_level1[, c(
#     "variable_id",
#     "Compound.name",
#     "CAS.ID",
#     "HMDB.ID",
#     "KEGG.ID",
#     "Lab.ID",
#     "Adduct",
#     "mz.error",
#     "RT.error",
#     "SS",
#     "Total.score",
#     "Database",
#     "Level"
#   )] %>%
#   dplyr::left_join(annotation_table[, c(
#     "variable_id",
#     "Compound.name",
#     "CAS.ID",
#     "HMDB.ID",
#     "KEGG.ID",
#     "Lab.ID",
#     "Adduct",
#     "mz.error",
#     "RT.error",
#     "SS",
#     "Total.score",
#     "Database",
#     "Level",
#     "compound_class",
#     "score",
#     "Formula"
#   )],
#   by = "variable_id")
#
# validation_result <-
#   validation_result %>%
#   plyr::dlply(.variables = plyr::.(variable_id)) %>%
#   purrr::map(function(x) {
#     x[, c("Compound.name.x",
#           "Compound.name.y",
#           "HMDB.ID.x",
#           "HMDB.ID.y",
#           "score")] %>%
#       dplyr::arrange(desc(score))
#   })
#
# validation_result[[1]]
#
# result <-
#   validation_result %>%
#   purrr::map(function(x) {
#     idx <-
#       which(x$HMDB.ID.x == x$HMDB.ID.y)
#     if (length(idx) == 0) {
#       temp <- x[1, ]
#       temp$Compound.name.y <- NA
#       temp$HMDB.ID.y <- NA
#       temp$score <- NA
#       temp$rank <- NA
#       return(temp)
#     }
#     x$rank <- idx
#     x[idx, ]
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
#
# result %>%
#   dplyr::mutate(idx = 1:nrow(result)) %>%
#   dplyr::arrange(rank) %>%
#   dplyr::mutate(idx = factor(idx, levels = idx)) %>%
#   ggplot(aes(idx, rank)) +
#   geom_point(aes(size = score)) +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank())
#
# sum(is.na(result$rank))
#
# result %>%
#   dplyr::filter(is.na(rank))
