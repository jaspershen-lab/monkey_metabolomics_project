no_source()
####data preparation

library(tidymass)
library(tidyverse)

###RPLC pos
rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/fMSEA/annotate_isotope.R")

#####load database
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/hmdb_ms1.rda"
)

setwd(r4projects::get_project_wd())

load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/RPLC/POS/object_rplc_pos"
)

load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/RPLC/NEG/object_rplc_neg"
)

load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/HILIC/POS/object_hilic_pos"
)

load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/HILIC/NEG/object_hilic_neg"
)

dir.create("3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation")
setwd("3-data_analysis/monkey_fMSEA_analysis/3_isotope_annotation")

###object_rplc_pos
###for one feature, if it have Level 1 or Level 2 annotation, remove Level 3 annotation
library(plyr)

object_rplc_pos@annotation_table <-
  object_rplc_pos@annotation_table %>%
  plyr::dlply(.variables = .(variable_id)) %>%
  purrr::map(function(x) {
    if (any(x$Level == 1) | any(x$Level == 2)) {
      x <-
        x %>%
        dplyr::filter(Level != 3)
    }
    
    if (any(x$Level == 1)) {
      x <-
        x %>%
        dplyr::filter(Level == 1)
    }
    return(x)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

annotation_table_rplc_pos <-
  object_rplc_pos@annotation_table %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Formula")],
                   by = "HMDB.ID")

variable_info_rplc_pos <-
  object_rplc_pos@variable_info

annotation_table_rplc_pos <-
  annotation_table_rplc_pos %>%
  dplyr::left_join(variable_info_rplc_pos[, c("variable_id", "mz", "rt")],
                   by = "variable_id")

annotation_table_rplc_pos <-
  annotation_table_rplc_pos %>%
  dplyr::filter(!is.na(Formula)) %>%
  dplyr::filter(nchar(Formula) > 3)

library(progress)

new_annotation_table_rplc_pos <-
  annotation_table_rplc_pos %>%
  dplyr::distinct(variable_id, Formula, Adduct, .keep_all = TRUE)

pb <-
  progress::progress_bar$new(total = nrow(new_annotation_table_rplc_pos))

isotope_table <-
  seq_len(nrow(new_annotation_table_rplc_pos)) %>%
  purrr::map(function(i) {
    pb$tick()
    temp <-
      annotate_isotope(
        formula = new_annotation_table_rplc_pos$Formula[i],
        adduct = new_annotation_table_rplc_pos$Adduct[i],
        feature_mz = new_annotation_table_rplc_pos$mz[i],
        feature_rt = new_annotation_table_rplc_pos$rt[i],
        variable_info = variable_info_rplc_pos,
        rt_tol = 5,
        mz_tol = 15,
        max_isotope = 3
      )
    if (!is.null(temp)) {
      temp <-
        data.frame(
          original_variable_id = new_annotation_table_rplc_pos$variable_id[i],
          Adduct = new_annotation_table_rplc_pos$Adduct[i],
          Formula = new_annotation_table_rplc_pos$Formula[i],
          temp
        )
    }
    return(temp)
    
  })

isotope_table <-
  isotope_table %>%
  dplyr::bind_rows()

isotope_table <-
  annotation_table_rplc_pos[, c(
    "variable_id",
    "Formula",
    "Adduct",
    "Compound.name",
    "CAS.ID",
    "HMDB.ID",
    "KEGG.ID",
    "Lab.ID",
    "Database"
  )] %>%
  dplyr::left_join(isotope_table,
                   by = c("variable_id" = "original_variable_id",
                          "Formula", "Adduct")) %>%
  dplyr::select(-c(feature_mz, feature_rt, variable_id)) %>%
  dplyr::rename(RT.error = rt.error,
                variable_id = variable_id.y) %>%
  dplyr::filter(!is.na(isotopes)) %>%
  dplyr::mutate(
    ms2_files_id = NA,
    ms2_spectrum_id = NA,
    mz.match.score = NA,
    RT.match.score = NA,
    CE = NA,
    SS = NA,
    Total.score = NA,
    Level = 3
  ) %>%
  dplyr::select(
    variable_id,
    ms2_files_id,
    ms2_spectrum_id,
    Compound.name,
    isotopes,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Lab.ID,
    Adduct,
    mz.error,
    mz.match.score,
    RT.error,
    RT.match.score,
    CE,
    SS,
    Total.score,
    Database,
    Level,
    Formula,
    mz,
    rt
  )

annotation_table_rplc_pos <-
  annotation_table_rplc_pos %>%
  dplyr::mutate(isotopes = "[M]") %>%
  dplyr::select(
    variable_id,
    ms2_files_id,
    ms2_spectrum_id,
    Compound.name,
    isotopes,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Lab.ID,
    Adduct,
    mz.error,
    mz.match.score,
    RT.error,
    RT.match.score,
    CE,
    SS,
    Total.score,
    Database,
    Level,
    Formula,
    mz,
    rt
  )

dim(annotation_table_rplc_pos)
dim(isotope_table)

annotation_table_rplc_pos <-
  rbind(annotation_table_rplc_pos,
        isotope_table)

dir.create("RPLC/POS/", showWarnings = FALSE, recursive = TRUE)
save(annotation_table_rplc_pos, file = "RPLC/POS/annotation_table_rplc_pos")

variable_info_rplc_pos <-
  object_rplc_pos@variable_info
save(variable_info_rplc_pos, file = "RPLC/POS/variable_info_rplc_pos")

sample_info_rplc_pos <-
  object_rplc_pos@sample_info
save(sample_info_rplc_pos, file = "RPLC/POS/sample_info_rplc_pos")

expression_data_rplc_pos <-
  object_rplc_pos@expression_data
save(expression_data_rplc_pos, file = "RPLC/POS/expression_data_rplc_pos")






###object_rplc_neg
###for one feature, if it have Level 1 or Level 2 annotation, remove Level 3 annotation
library(plyr)
object_rplc_neg@annotation_table <-
  object_rplc_neg@annotation_table %>%
  plyr::dlply(.variables = .(variable_id)) %>%
  purrr::map(function(x) {
    if (any(x$Level == 1) | any(x$Level == 2)) {
      x <-
        x %>%
        dplyr::filter(Level != 3)
    }
    
    if (any(x$Level == 1)) {
      x <-
        x %>%
        dplyr::filter(Level == 1)
    }
    return(x)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

annotation_table_rplc_neg <-
  object_rplc_neg@annotation_table %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Formula")],
                   by = "HMDB.ID")

variable_info_rplc_neg <-
  object_rplc_neg@variable_info

annotation_table_rplc_neg <-
  annotation_table_rplc_neg %>%
  dplyr::left_join(variable_info_rplc_neg[, c("variable_id", "mz", "rt")],
                   by = "variable_id")

annotation_table_rplc_neg <-
  annotation_table_rplc_neg %>%
  dplyr::filter(!is.na(Formula)) %>%
  dplyr::filter(nchar(Formula) > 3)

library(progress)

new_annotation_table_rplc_neg <-
  annotation_table_rplc_neg %>%
  dplyr::distinct(variable_id, Formula, Adduct, .keep_all = TRUE)

pb <-
  progress::progress_bar$new(total = nrow(new_annotation_table_rplc_neg))

isotope_table <-
  seq_len(nrow(new_annotation_table_rplc_neg)) %>%
  purrr::map(function(i) {
    pb$tick()
    temp <-
      annotate_isotope(
        formula = new_annotation_table_rplc_neg$Formula[i],
        adduct = new_annotation_table_rplc_neg$Adduct[i],
        feature_mz = new_annotation_table_rplc_neg$mz[i],
        feature_rt = new_annotation_table_rplc_neg$rt[i],
        variable_info = variable_info_rplc_neg,
        rt_tol = 5,
        mz_tol = 15,
        max_isotope = 3
      )
    if (!is.null(temp)) {
      temp <-
        data.frame(
          original_variable_id = new_annotation_table_rplc_neg$variable_id[i],
          Adduct = new_annotation_table_rplc_neg$Adduct[i],
          Formula = new_annotation_table_rplc_neg$Formula[i],
          temp
        )
    }
    return(temp)
    
  })

isotope_table <-
  isotope_table %>%
  dplyr::bind_rows()

isotope_table <-
  annotation_table_rplc_neg[, c(
    "variable_id",
    "Formula",
    "Adduct",
    "Compound.name",
    "CAS.ID",
    "HMDB.ID",
    "KEGG.ID",
    "Lab.ID",
    "Database"
  )] %>%
  dplyr::left_join(isotope_table,
                   by = c("variable_id" = "original_variable_id",
                          "Formula", "Adduct")) %>%
  dplyr::select(-c(feature_mz, feature_rt, variable_id)) %>%
  dplyr::rename(RT.error = rt.error,
                variable_id = variable_id.y) %>%
  dplyr::filter(!is.na(isotopes)) %>%
  dplyr::mutate(
    ms2_files_id = NA,
    ms2_spectrum_id = NA,
    mz.match.score = NA,
    RT.match.score = NA,
    CE = NA,
    SS = NA,
    Total.score = NA,
    Level = 3
  ) %>%
  dplyr::select(
    variable_id,
    ms2_files_id,
    ms2_spectrum_id,
    Compound.name,
    isotopes,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Lab.ID,
    Adduct,
    mz.error,
    mz.match.score,
    RT.error,
    RT.match.score,
    CE,
    SS,
    Total.score,
    Database,
    Level,
    Formula,
    mz,
    rt
  )

annotation_table_rplc_neg <-
  annotation_table_rplc_neg %>%
  dplyr::mutate(isotopes = "[M]") %>%
  dplyr::select(
    variable_id,
    ms2_files_id,
    ms2_spectrum_id,
    Compound.name,
    isotopes,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Lab.ID,
    Adduct,
    mz.error,
    mz.match.score,
    RT.error,
    RT.match.score,
    CE,
    SS,
    Total.score,
    Database,
    Level,
    Formula,
    mz,
    rt
  )

dim(annotation_table_rplc_neg)
dim(isotope_table)

annotation_table_rplc_neg <-
  rbind(annotation_table_rplc_neg,
        isotope_table)

dir.create("RPLC/NEG/", showWarnings = FALSE, recursive = TRUE)
save(annotation_table_rplc_neg, file = "RPLC/NEG/annotation_table_rplc_neg")

variable_info_rplc_neg <-
  object_rplc_neg@variable_info
save(variable_info_rplc_neg, file = "RPLC/NEG/variable_info_rplc_neg")

sample_info_rplc_neg <-
  object_rplc_neg@sample_info
save(sample_info_rplc_neg, file = "RPLC/NEG/sample_info_rplc_neg")

expression_data_rplc_neg <-
  object_rplc_neg@expression_data
save(expression_data_rplc_neg, file = "RPLC/NEG/expression_data_rplc_neg")


###object_hilic_pos
###for one feature, if it have Level 1 or Level 2 annotation, remove Level 3 annotation
library(plyr)
object_hilic_pos@annotation_table <-
  object_hilic_pos@annotation_table %>%
  plyr::dlply(.variables = .(variable_id)) %>%
  purrr::map(function(x) {
    if (any(x$Level == 1) | any(x$Level == 2)) {
      x <-
        x %>%
        dplyr::filter(Level != 3)
    }
    
    if (any(x$Level == 1)) {
      x <-
        x %>%
        dplyr::filter(Level == 1)
    }
    return(x)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

annotation_table_hilic_pos <-
  object_hilic_pos@annotation_table %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Formula")],
                   by = "HMDB.ID")

variable_info_hilic_pos <-
  object_hilic_pos@variable_info

annotation_table_hilic_pos <-
  annotation_table_hilic_pos %>%
  dplyr::left_join(variable_info_hilic_pos[, c("variable_id", "mz", "rt")],
                   by = "variable_id")

annotation_table_hilic_pos <-
  annotation_table_hilic_pos %>%
  dplyr::filter(!is.na(Formula)) %>%
  dplyr::filter(nchar(Formula) > 3)

library(progress)

new_annotation_table_hilic_pos <-
  annotation_table_hilic_pos %>%
  dplyr::distinct(variable_id, Formula, Adduct, .keep_all = TRUE)

pb <-
  progress::progress_bar$new(total = nrow(new_annotation_table_hilic_pos))

isotope_table <-
  seq_len(nrow(new_annotation_table_hilic_pos)) %>%
  purrr::map(function(i) {
    pb$tick()
    temp <-
      annotate_isotope(
        formula = new_annotation_table_hilic_pos$Formula[i],
        adduct = new_annotation_table_hilic_pos$Adduct[i],
        feature_mz = new_annotation_table_hilic_pos$mz[i],
        feature_rt = new_annotation_table_hilic_pos$rt[i],
        variable_info = variable_info_hilic_pos,
        rt_tol = 5,
        mz_tol = 15,
        max_isotope = 3
      )
    if (!is.null(temp)) {
      temp <-
        data.frame(
          original_variable_id = new_annotation_table_hilic_pos$variable_id[i],
          Adduct = new_annotation_table_hilic_pos$Adduct[i],
          Formula = new_annotation_table_hilic_pos$Formula[i],
          temp
        )
    }
    return(temp)
    
  })

isotope_table <-
  isotope_table %>%
  dplyr::bind_rows()

isotope_table <-
  annotation_table_hilic_pos[, c(
    "variable_id",
    "Formula",
    "Adduct",
    "Compound.name",
    "CAS.ID",
    "HMDB.ID",
    "KEGG.ID",
    "Lab.ID",
    "Database"
  )] %>%
  dplyr::left_join(isotope_table,
                   by = c("variable_id" = "original_variable_id",
                          "Formula", "Adduct")) %>%
  dplyr::select(-c(feature_mz, feature_rt, variable_id)) %>%
  dplyr::rename(RT.error = rt.error,
                variable_id = variable_id.y) %>%
  dplyr::filter(!is.na(isotopes)) %>%
  dplyr::mutate(
    ms2_files_id = NA,
    ms2_spectrum_id = NA,
    mz.match.score = NA,
    RT.match.score = NA,
    CE = NA,
    SS = NA,
    Total.score = NA,
    Level = 3
  ) %>%
  dplyr::select(
    variable_id,
    ms2_files_id,
    ms2_spectrum_id,
    Compound.name,
    isotopes,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Lab.ID,
    Adduct,
    mz.error,
    mz.match.score,
    RT.error,
    RT.match.score,
    CE,
    SS,
    Total.score,
    Database,
    Level,
    Formula,
    mz,
    rt
  )

annotation_table_hilic_pos <-
  annotation_table_hilic_pos %>%
  dplyr::mutate(isotopes = "[M]") %>%
  dplyr::select(
    variable_id,
    ms2_files_id,
    ms2_spectrum_id,
    Compound.name,
    isotopes,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Lab.ID,
    Adduct,
    mz.error,
    mz.match.score,
    RT.error,
    RT.match.score,
    CE,
    SS,
    Total.score,
    Database,
    Level,
    Formula,
    mz,
    rt
  )

dim(annotation_table_hilic_pos)
dim(isotope_table)

annotation_table_hilic_pos <-
  rbind(annotation_table_hilic_pos,
        isotope_table)

dir.create("HILIC/POS/", showWarnings = FALSE, recursive = TRUE)
save(annotation_table_hilic_pos, file = "HILIC/POS/annotation_table_hilic_pos")

variable_info_hilic_pos <-
  object_hilic_pos@variable_info
save(variable_info_hilic_pos, file = "HILIC/POS/variable_info_hilic_pos")

sample_info_hilic_pos <-
  object_hilic_pos@sample_info
save(sample_info_hilic_pos, file = "HILIC/POS/sample_info_hilic_pos")

expression_data_hilic_pos <-
  object_hilic_pos@expression_data
save(expression_data_hilic_pos, file = "HILIC/POS/expression_data_hilic_pos")






###object_hilic_neg
###for one feature, if it have Level 1 or Level 2 annotation, remove Level 3 annotation
library(plyr)
object_hilic_neg@annotation_table <-
  object_hilic_neg@annotation_table %>%
  plyr::dlply(.variables = .(variable_id)) %>%
  purrr::map(function(x) {
    if (any(x$Level == 1) | any(x$Level == 2)) {
      x <-
        x %>%
        dplyr::filter(Level != 3)
    }
    
    if (any(x$Level == 1)) {
      x <-
        x %>%
        dplyr::filter(Level == 1)
    }
    return(x)
  }) %>%
  dplyr::bind_rows() %>%
  as.data.frame()

annotation_table_hilic_neg <-
  object_hilic_neg@annotation_table %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Formula")],
                   by = "HMDB.ID")

variable_info_hilic_neg <-
  object_hilic_neg@variable_info

annotation_table_hilic_neg <-
  annotation_table_hilic_neg %>%
  dplyr::left_join(variable_info_hilic_neg[, c("variable_id", "mz", "rt")],
                   by = "variable_id")

annotation_table_hilic_neg <-
  annotation_table_hilic_neg %>%
  dplyr::filter(!is.na(Formula)) %>%
  dplyr::filter(nchar(Formula) > 3)

library(progress)

new_annotation_table_hilic_neg <-
  annotation_table_hilic_neg %>%
  dplyr::distinct(variable_id, Formula, Adduct, .keep_all = TRUE)

pb <-
  progress::progress_bar$new(total = nrow(new_annotation_table_hilic_neg))

isotope_table <-
  seq_len(nrow(new_annotation_table_hilic_neg)) %>%
  purrr::map(function(i) {
    pb$tick()
    temp <-
      annotate_isotope(
        formula = new_annotation_table_hilic_neg$Formula[i],
        adduct = new_annotation_table_hilic_neg$Adduct[i],
        feature_mz = new_annotation_table_hilic_neg$mz[i],
        feature_rt = new_annotation_table_hilic_neg$rt[i],
        variable_info = variable_info_hilic_neg,
        rt_tol = 5,
        mz_tol = 15,
        max_isotope = 3
      )
    if (!is.null(temp)) {
      temp <-
        data.frame(
          original_variable_id = new_annotation_table_hilic_neg$variable_id[i],
          Adduct = new_annotation_table_hilic_neg$Adduct[i],
          Formula = new_annotation_table_hilic_neg$Formula[i],
          temp
        )
    }
    return(temp)
    
  })

isotope_table <-
  isotope_table %>%
  dplyr::bind_rows()

isotope_table <-
  annotation_table_hilic_neg[, c(
    "variable_id",
    "Formula",
    "Adduct",
    "Compound.name",
    "CAS.ID",
    "HMDB.ID",
    "KEGG.ID",
    "Lab.ID",
    "Database"
  )] %>%
  dplyr::left_join(isotope_table,
                   by = c("variable_id" = "original_variable_id",
                          "Formula", "Adduct")) %>%
  dplyr::select(-c(feature_mz, feature_rt, variable_id)) %>%
  dplyr::rename(RT.error = rt.error,
                variable_id = variable_id.y) %>%
  dplyr::filter(!is.na(isotopes)) %>%
  dplyr::mutate(
    ms2_files_id = NA,
    ms2_spectrum_id = NA,
    mz.match.score = NA,
    RT.match.score = NA,
    CE = NA,
    SS = NA,
    Total.score = NA,
    Level = 3
  ) %>%
  dplyr::select(
    variable_id,
    ms2_files_id,
    ms2_spectrum_id,
    Compound.name,
    isotopes,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Lab.ID,
    Adduct,
    mz.error,
    mz.match.score,
    RT.error,
    RT.match.score,
    CE,
    SS,
    Total.score,
    Database,
    Level,
    Formula,
    mz,
    rt
  )

annotation_table_hilic_neg <-
  annotation_table_hilic_neg %>%
  dplyr::mutate(isotopes = "[M]") %>%
  dplyr::select(
    variable_id,
    ms2_files_id,
    ms2_spectrum_id,
    Compound.name,
    isotopes,
    CAS.ID,
    HMDB.ID,
    KEGG.ID,
    Lab.ID,
    Adduct,
    mz.error,
    mz.match.score,
    RT.error,
    RT.match.score,
    CE,
    SS,
    Total.score,
    Database,
    Level,
    Formula,
    mz,
    rt
  )

dim(annotation_table_hilic_neg)
dim(isotope_table)

annotation_table_hilic_neg <-
  rbind(annotation_table_hilic_neg,
        isotope_table)

dir.create("HILIC/NEG/", showWarnings = FALSE, recursive = TRUE)
save(annotation_table_hilic_neg, file = "HILIC/NEG/annotation_table_hilic_neg")

variable_info_hilic_neg <-
  object_hilic_neg@variable_info
save(variable_info_hilic_neg, file = "HILIC/NEG/variable_info_hilic_neg")

sample_info_hilic_neg <-
  object_hilic_neg@sample_info
save(sample_info_hilic_neg, file = "HILIC/NEG/sample_info_hilic_neg")

expression_data_hilic_neg <-
  object_hilic_neg@expression_data
save(expression_data_hilic_neg, file = "HILIC/NEG/expression_data_hilic_neg")
