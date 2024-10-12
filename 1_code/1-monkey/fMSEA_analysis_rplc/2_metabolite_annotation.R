no_source()
####data preparation

library(tidymass)
library(tidyverse)

###RPLC pos
rm(list = ls())

setwd(r4projects::get_project_wd())

source("1-code/fMSEA/annotate_metabolite4fmsea.R")
source("1-code/fMSEA/data_checking.R")

load(
  "3-data_analysis/monkey_fMSEA_analysis_rplc/1_data_preparation/RPLC/POS/object_rplc_pos"
)

object_rplc_pos

load(
  "3-data_analysis/monkey_fMSEA_analysis_rplc/1_data_preparation/RPLC/NEG/object_rplc_neg"
)

object_rplc_neg

setwd("3-data_analysis/monkey_fMSEA_analysis_rplc/2_metabolite_annotation/")

object_list <-
  list(object_rplc_pos,
       object_rplc_neg)

rm(list = c("object_rplc_pos", "object_rplc_neg"))
gc()

###metabolite annotation
object_list <-
  annotate_metabolite4fmsea(
    object_list = object_list,
    use_default_database = TRUE,
    rt.match.tol = 30,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01
  )

save(object_list, file = "object_list")

# load("../5_redundant_removal/temp_data")
# 
# temp_data1 <- 
# temp_data %>% 
#   dplyr::filter(stringr::str_detect(Adduct, "\\(M\\+H-2H2O\\)\\+"))
# 
# seq_len(nrow(temp_data1)) %>% 
# purrr::map(function(i){
#   object_list[[1]]@annotation_table %>% 
#     dplyr::filter(HMDB.ID == temp_data1$HMDB.ID[i] & Adduct == temp_data1$Adduct[i])  
# }) %>% 
#   dplyr::bind_rows()




