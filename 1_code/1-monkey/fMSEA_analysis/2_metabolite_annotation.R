no_source()
####data preparation

library(tidymass)
library(tidyverse)

###RPLC pos
rm(list = ls())

setwd(r4projects::get_project_wd())

#####load database
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/hmdb_ms2.rda"
)
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/massbank_ms2.rda"
)
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/metlin_ms2.rda"
)
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/mona_ms2.rda"
)
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/mpsnyder_hilic_ms2.rda"
)
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/mpsnyder_rplc_ms2.rda"
)
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/nist_ms2.rda"
)
load(
  "3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/database/hmdb_ms1.rda"
)

###only remain the metabolites have HMDB.ID, CAS.ID or KEGG.ID
hmdb_ms2 <-
  hmdb_ms2 %>%
  filter(!is.na(HMDB.ID))

massbank_ms2 <-
  massbank_ms2 %>%
  filter(!is.na(HMDB.ID))

metlin_ms2 <-
  metlin_ms2 %>%
  filter(!is.na(HMDB.ID))

mpsnyder_hilic_ms2 <-
  mpsnyder_hilic_ms2 %>%
  filter(!is.na(HMDB.ID))

mpsnyder_rplc_ms2 <-
  mpsnyder_rplc_ms2 %>%
  filter(!is.na(HMDB.ID))

mona_ms2 <-
  mona_ms2 %>%
  filter(!is.na(HMDB.ID))

nist_ms2 <-
  nist_ms2 %>%
  filter(!is.na(HMDB.ID))

setwd(r4projects::get_project_wd())

load(
  "3-data_analysis/monkey_fMSEA_analysis/1_data_preparation/RPLC/POS/object_rplc_pos"
)

object_rplc_pos

setwd("3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/RPLC/POS/")

###mpsnyder_rplc
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = mpsnyder_rplc_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

##MassBank
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = massbank_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###METLIN
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = metlin_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###MONA
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = mona_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###NIST
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = nist_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###HMDB MS2
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = hmdb_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )


###HMDB MS1
object_rplc_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_pos,
    polarity = "positive",
    rt.match.tol = 100000,
    column = "rp",
    candidate.num = 100,
    database = hmdb_ms1,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

save(object_rplc_pos, file = "object_rplc_pos")


###RPLC neg
setwd(r4projects::get_project_wd())
load("3-data_analysis/monkey_metabolomics_data_cleaning/RPLC/NEG/object_neg2")
object_neg2

setwd("3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/RPLC/NEG/")

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
    candidate.num = 100,
    database = mpsnyder_rplc_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

##MassBank
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = massbank_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###METLIN
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = metlin_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###MONA
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = mona_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###NIST
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = nist_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )


###HMDB
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = hmdb_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###HMDB MS1
object_rplc_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_rplc_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "rp",
    candidate.num = 100,
    database = hmdb_ms1,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

save(object_rplc_neg, file = "object_rplc_neg")



setwd(r4projects::get_project_wd())
load("3-data_analysis/monkey_metabolomics_data_cleaning/HILIC/POS/object_pos2")

setwd("3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/HILIC/POS/")
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
    candidate.num = 100,
    database = mpsnyder_hilic_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

##MassBank
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = massbank_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###METLIN
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = metlin_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###MONA
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = mona_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###NIST
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = nist_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###HMDB
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = hmdb_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###HMDB MS1
object_hilic_pos <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_pos,
    polarity = "positive",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = hmdb_ms1,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

save(object_hilic_pos, file = "object_hilic_pos")


###HILIC neg
setwd(r4projects::get_project_wd())
load("3-data_analysis/monkey_metabolomics_data_cleaning/HILIC/NEG/object_neg2")
object_neg2

setwd("3-data_analysis/monkey_metabolomics_metabolite_annotation/HILIC/NEG/")
setwd("3-data_analysis/monkey_fMSEA_analysis/2_metabolite_annotation/HILIC/NEG/")

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
    candidate.num = 100,
    database = mpsnyder_hilic_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

##MassBank
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = massbank_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###METLIN
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = metlin_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###MONA
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = mona_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###NIST
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = nist_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )


###HMDB
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = hmdb_ms2,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

###HMDB MS1
object_hilic_neg <-
  metid::annotate_metabolites_mass_dataset(
    object = object_hilic_neg,
    polarity = "negative",
    rt.match.tol = 30,
    column = "hilic",
    candidate.num = 100,
    database = hmdb_ms1,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

save(object_hilic_neg, file = "object_hilic_neg")
