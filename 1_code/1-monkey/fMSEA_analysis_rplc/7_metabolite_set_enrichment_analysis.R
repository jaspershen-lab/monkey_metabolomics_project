no_source()
####data preparation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/fMSEA/calculate_redundancy.R")
source("1-code/fMSEA/run_msea.R")
source("1-code/tools.R")
source("1-code/fMSEA/whole_plot.R")

setwd("3-data_analysis/monkey_fMSEA_analysis_rplc/whole_workflow")

load(file.path("redundancy_removing", "annotation_table"))
load(file.path("isotope_annotation", "object_list"))

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

temp <-
  data.table::data.table(annotation_table)[, (number = length(unique(compound_class))), by = .(HMDB.ID)]

temp <-
  temp[which(temp$N > 1),]

temp

###use several pathways as example
pathway <- hmdb_pathway[1:2]

plot_whole(annotation_table = annotation_table,
           pathway = hmdb_pathway[1:3])

annotation_table$variable_id

###Histidine Metabolism
metabolite_set <-
  hmdb_pathway

feature_table <- annotation_table

feature_table <-
  annotation_table %>%
  dplyr::rename(condition = spearman_cor) %>%
  dplyr::arrange(dplyr::desc(condition))

#########fMSEA
###this is the weight for runing enrichment score
exponent = 1
perm_num = 1000
min_size = 5
max_size = 1000
pvalue_cutoff = 0.2
p_adjust_method = "fdr"
seed = FALSE
verbose = TRUE

result <-
  run_msea(
    feature_table = feature_table,
    metabolite_set = metabolite_set,
    exponent = 1,
    perm_num = 1000,
    min_size = 5,
    max_size = 1000,
    pvalue_cutoff = 0.25,
    p_adjust_method = "fdr",
    seed = TRUE,
    verbose = TRUE
  )
#
# msea_plot(x = result, title = result@result$Description[1])
#
# msea_plot(x = head_result[[idx]],
#           title = paste(head_result[[idx]]@result$p_adjust))
#
#
# msea_plot(x = tail_result[[idx]],
#           title = paste(tail_result[[idx]]@result$p_adjust))
