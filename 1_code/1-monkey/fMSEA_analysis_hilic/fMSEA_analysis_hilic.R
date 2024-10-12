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
source("1-code/fMSEA/whole_plot.R")

load(
  "3-data_analysis/monkey_fMSEA_analysis/1_data_preparation/HILIC/POS/object_hilic_pos"
)

object_hilic_pos

load(
  "3-data_analysis/monkey_fMSEA_analysis/1_data_preparation/HILIC/NEG/object_hilic_neg"
)

object_hilic_neg

dir.create("3-data_analysis/monkey_fMSEA_analysis_hilic/whole_workflow",
           recursive = TRUE)
setwd("3-data_analysis/monkey_fMSEA_analysis_hilic/whole_workflow")

object_list <-
  list(object_hilic_pos,
       object_hilic_neg)

# rm(list = c("object_rplc_pos", "object_rplc_neg"))
# gc()

#####1. Metabolite annotation
object_list <-
  annotate_metabolite4fmsea(
    object_list = object_list,
    use_default_database = TRUE,
    rt.match.tol = 30,
    threads = 3,
    remove_fragment_intensity_cutoff = 0.01,
    path = "metabolite_annotation"
  )

###2. Isotope annotation
object_list <-
  annotate_isotopes4fmsea(
    object_list = object_list,
    threads = 3,
    rt_tol = 10,
    mz_tol = 15,
    intensity_tol = 0.5,
    max_isotope = 2
  )

dir.create("isotope_annotation")

save(object_list,
     file = file.path("isotope_annotation", "object_list"))

load("isotope_annotation/object_list")

# ####3. grouping features
# annotation_table <-
#   group_features4fmsea(object_list = object_list,
#                        rt.match.tol = 10,
#                        threads = 3)
# 
# dir.create("feature_grouping")
# 
# save(annotation_table,
#      file = file.path("feature_grouping", "annotation_table"))

load(file.path("feature_grouping", "annotation_table"))

# ## 4. redundancy removing
# annotation_table <-
#   remove_redundancy4fmsea(
#     annotation_table = annotation_table,
#     score_cutoff = 20,
#     threads = 3
#   )
# 
# dir.create("redundancy_removing")
# 
# save(annotation_table,
#      file = file.path("redundancy_removing", "annotation_table"))

load(file.path("redundancy_removing", "annotation_table"))

calculate_redundancy(annotation_table = annotation_table)

#### 5. optimization
dim(annotation_table)

length(unique(annotation_table$variable_id))
length(unique(annotation_table$compound_class))
length(unique(annotation_table$HMDB.ID))

library(metpath)
data("hmdb_pathway")

#get the class of pathways
pathway_class <-
  metpath::pathway_class(hmdb_pathway)

remain_idx <-
  which(unlist(pathway_class) == "Metabolic;primary_pathway")

remain_idx

hmdb_pathway <-
  hmdb_pathway[remain_idx]

calculate_redundancy(annotation_table = annotation_table)

unique(annotation_table$variable_id)
unique(annotation_table$HMDB.ID)

library(data.table)

# temp <-
#   data.table::data.table(annotation_table)[, .N, by = .(variable_id, HMDB.ID)]
#
# temp <-
#   temp[which(temp$N > 1), ]
#
# temp

temp <-
  data.table::data.table(annotation_table)[, (number = length(unique(compound_class))), by = .(HMDB.ID)]

temp <-
  temp[which(temp$N > 1),]

temp


annotation_table %>%
  dplyr::count(variable_id) %>%
  ggplot(aes(n)) +
  geom_histogram(color = "black", binwidth = 3)

annotation_table %>%
  dplyr::count(variable_id) %>%
  dplyr::pull(n) %>%
  median()

annotation_table %>%
  dplyr::count(variable_id) %>%
  dplyr::pull(n) %>%
  mean()

###use several pathways as example
pathway <- hmdb_pathway[1:2]

plot_whole(annotation_table = annotation_table,
           pathway = hmdb_pathway[1:3])



######
