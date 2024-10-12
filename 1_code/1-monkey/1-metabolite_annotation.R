###
no_source

####This should be run in workstation

library(tidymass)
library(tidyverse)

###RPLC pos
rm(list = ls())

setwd(r4projects::get_project_wd())
load("3-data_analysis/monkey_metabolomics_data_cleaning/RPLC/POS/object_pos2")
object_pos2

load("3-data_analysis/monkey_metabolomics_metabolite_annotation/hmdb_ms2.rda")
load("3-data_analysis/monkey_metabolomics_metabolite_annotation/massbank_ms2.rda")
load("3-data_analysis/monkey_metabolomics_metabolite_annotation/metlin_ms2.rda")
load("3-data_analysis/monkey_metabolomics_metabolite_annotation/mona_ms2.rda")
load(
  "3-data_analysis/monkey_metabolomics_metabolite_annotation/mpsnyder_hilic_ms2.rda"
)
load(
  "3-data_analysis/monkey_metabolomics_metabolite_annotation/mpsnyder_rplc_ms2.rda"
)
load("3-data_analysis/monkey_metabolomics_metabolite_annotation/nist_ms2.rda")

###only remain the metabolites have HMDB.ID, CAS.ID or KEGG.ID
hmdb_ms2 <-
  hmdb_ms2 %>%
  filter(!is.na(HMDB.ID) | !is.na(KEGG.ID) | !is.na(CAS.ID))

massbank_ms2 <-
  massbank_ms2 %>%
  filter(!is.na(HMDB.ID) | !is.na(KEGG.ID) | !is.na(CAS.ID))

metlin_ms2 <-
  metlin_ms2 %>%
  filter(!is.na(HMDB.ID) | !is.na(KEGG.ID) | !is.na(CAS.ID))

mpsnyder_hilic_ms2 <-
  mpsnyder_hilic_ms2 %>%
  filter(!is.na(HMDB.ID) | !is.na(KEGG.ID) | !is.na(CAS.ID))

mpsnyder_rplc_ms2 <-
  mpsnyder_rplc_ms2 %>%
  filter(!is.na(HMDB.ID) | !is.na(KEGG.ID) | !is.na(CAS.ID))

mona_ms2 <-
  mona_ms2 %>%
  filter(!is.na(HMDB.ID) | !is.na(KEGG.ID) | !is.na(CAS.ID))

nist_ms2 <-
  nist_ms2 %>%
  filter(!is.na(HMDB.ID) | !is.na(KEGG.ID) | !is.na(CAS.ID))

setwd("3-data_analysis/monkey_metabolomics_metabolite_annotation/RPLC/POS/")



object_rplc_pos <-
  object_pos2 %>%
  massdataset::mutate_ms2(
    column = "rp",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "."
  )

###mpsnyder_rplc
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = mpsnyder_rplc_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

##MassBank
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = massbank_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###METLIN
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = metlin_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###MONA
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = mona_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###NIST
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = nist_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )


###HMDB
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = hmdb_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

save(object_rplc_pos, file = "object_rplc_pos")


####output MS2 matching plot
annotation_table <-
  extract_annotation_table(object = object_rplc_pos)

dir.create("ms2_plot")

seq_len(nrow(annotation_table)) %>%
  purrr::walk(function(i) {
    cat(i, " ")
    database_name <-
      annotation_table$Database[i]
    database <-
      switch(
        EXPR = database_name,
        "MoNA_2022-04-27" = mona_ms2,
        "MassBank_2022-04-27" = massbank_ms2,
        "HMDB_2022-04-11" = hmdb_ms2,
        "NIST_20220425" = nist_ms2,
        "Michael_Snyder_RPLC_20220424" = mpsnyder_rplc_ms2,
        "METLIN_20220425" = metlin_ms2
      )
    
    plot <-
      metid::ms2_plot_mass_dataset(
        object = object_rplc_pos,
        variable_id = annotation_table$variable_id[i],
        polarity = "positive",
        database = database
      )
    if (class(plot) == "list") {
      plot <- plot[[1]]
    }
    
    ggsave(
      plot,
      filename = file.path(
        "ms2_plot",
        paste0(
          annotation_table$variable_id[i],
          "_",
          annotation_table$Database[i],
          ".pdf"
        )
      ),
      width = 9,
      height = 7
    )
  })


###RPLC neg

setwd(r4projects::get_project_wd())
load("3-data_analysis/monkey_metabolomics_data_cleaning/RPLC/NEG/object_neg2")
object_neg2

setwd("3-data_analysis/monkey_metabolomics_metabolite_annotation/RPLC/NEG/")

object_rplc_neg <-
  object_neg2 %>%
  massdataset::mutate_ms2(
    column = "rp",
    polarity = "negative",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "."
  )

###mpsnyder_rplc
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = mpsnyder_rplc_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

##MassBank
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = massbank_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###METLIN
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = metlin_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###MONA
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = mona_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###NIST
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = nist_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )


###HMDB
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 3,
    database = hmdb_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

save(object_rplc_neg, file = "object_rplc_neg")






####output MS2 matching plot
annotation_table <-
  extract_annotation_table(object = object_rplc_neg)

dir.create("ms2_plot")

seq_len(nrow(annotation_table)) %>%
  purrr::walk(function(i) {
    cat(i, " ")
    database_name <-
      annotation_table$Database[i]
    database <-
      switch(
        EXPR = database_name,
        "MoNA_2022-04-27" = mona_ms2,
        "MassBank_2022-04-27" = massbank_ms2,
        "HMDB_2022-04-11" = hmdb_ms2,
        "NIST_20220425" = nist_ms2,
        "Michael_Snyder_RPLC_20220424" = mpsnyder_rplc_ms2,
        "METLIN_20220425" = metlin_ms2
      )
    
    plot <-
      metid::ms2_plot_mass_dataset(
        object = object_rplc_neg,
        variable_id = annotation_table$variable_id[i],
        polarity = "negative",
        database = database
      )
    if (class(plot) == "list") {
      plot <- plot[[1]]
    }
    
    ggsave(
      plot,
      filename = file.path(
        "ms2_plot",
        paste0(
          annotation_table$variable_id[i],
          "_",
          annotation_table$Database[i],
          ".pdf"
        )
      ),
      width = 9,
      height = 7
    )
  })









setwd(r4projects::get_project_wd())
load("3-data_analysis/monkey_metabolomics_data_cleaning/HILIC/POS/object_pos2")
setwd("3-data_analysis/monkey_metabolomics_metabolite_annotation/HILIC/POS/")

object_hilic_pos <-
  object_pos2 %>%
  massdataset::mutate_ms2(
    column = "hilic",
    polarity = "positive",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "."
  )

###mpsnyder_hilic
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = mpsnyder_hilic_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

##MassBank
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = massbank_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###METLIN
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = metlin_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###MONA
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = mona_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###NIST
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = nist_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###HMDB
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = hmdb_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

save(object_hilic_pos, file = "object_hilic_pos")


####output MS2 matching plot
annotation_table <-
  extract_annotation_table(object = object_hilic_pos)

dir.create("ms2_plot")

seq_len(nrow(annotation_table)) %>%
  purrr::walk(function(i) {
    cat(i, " ")
    database_name <-
      annotation_table$Database[i]
    database <-
      switch(
        EXPR = database_name,
        "MoNA_2022-04-27" = mona_ms2,
        "MassBank_2022-04-27" = massbank_ms2,
        "HMDB_2022-04-11" = hmdb_ms2,
        "NIST_20220425" = nist_ms2,
        "Michael_Snyder_HILIC_20220424" = mpsnyder_hilic_ms2,
        "METLIN_20220425" = metlin_ms2
      )
    
    plot <-
      metid::ms2_plot_mass_dataset(
        object = object_hilic_pos,
        variable_id = annotation_table$variable_id[i],
        polarity = "positive",
        database = database
      )
    if (class(plot) == "list") {
      plot <- plot[[1]]
    }
    
    ggsave(
      plot,
      filename = file.path(
        "ms2_plot",
        paste0(
          annotation_table$variable_id[i],
          "_",
          annotation_table$Database[i],
          ".pdf"
        )
      ),
      width = 9,
      height = 7
    )
  })


###HILIC neg

setwd(r4projects::get_project_wd())
load("3-data_analysis/monkey_metabolomics_data_cleaning/HILIC/NEG/object_neg2")
object_neg2

setwd("3-data_analysis/monkey_metabolomics_metabolite_annotation/HILIC/NEG/")

object_hilic_neg <-
  object_neg2 %>%
  massdataset::mutate_ms2(
    column = "hilic",
    polarity = "negative",
    ms1.ms2.match.mz.tol = 15,
    ms1.ms2.match.rt.tol = 30,
    path = "."
  )

###mpsnyder_hilic
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = mpsnyder_hilic_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

##MassBank
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = massbank_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###METLIN
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = metlin_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###MONA
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = mona_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

###NIST
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = nist_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )


###HMDB
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 3,
    database = hmdb_ms2,
    threads = 10,
    remove_fragment_intensity_cutoff = 0.01
  )

save(object_hilic_neg, file = "object_hilic_neg")






####output MS2 matching plot
annotation_table <-
  extract_annotation_table(object = object_hilic_neg)

dir.create("ms2_plot")

seq_len(nrow(annotation_table)) %>%
  purrr::walk(function(i) {
    cat(i, " ")
    database_name <-
      annotation_table$Database[i]
    database <-
      switch(
        EXPR = database_name,
        "MoNA_2022-04-27" = mona_ms2,
        "MassBank_2022-04-27" = massbank_ms2,
        "HMDB_2022-04-11" = hmdb_ms2,
        "NIST_20220425" = nist_ms2,
        "Michael_Snyder_HILIC_20220424" = mpsnyder_hilic_ms2,
        "METLIN_20220425" = metlin_ms2
      )
    
    plot <-
      metid::ms2_plot_mass_dataset(
        object = object_hilic_neg,
        variable_id = annotation_table$variable_id[i],
        polarity = "negative",
        database = database
      )
    if (class(plot) == "list") {
      plot <- plot[[1]]
    }
    
    ggsave(
      plot,
      filename = file.path(
        "ms2_plot",
        paste0(
          annotation_table$variable_id[i],
          "_",
          annotation_table$Database[i],
          ".pdf"
        )
      ),
      width = 9,
      height = 7
    )
  })
