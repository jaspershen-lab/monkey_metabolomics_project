no_source()
####data preparation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1-code/fMSEA/remove_redundancy.R")
source("1-code/fMSEA/calculate_redundancy.R")
source("1-code/tools.R")

load("3-data_analysis/monkey_fMSEA_analysis/5_redundant_removal/annotation_table")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/variable_info_hilic_pos")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/variable_info_hilic_neg")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/variable_info_rplc_pos")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/variable_info_rplc_neg")
# 
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/sample_info_hilic_pos")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/sample_info_hilic_neg")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/sample_info_rplc_pos")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/sample_info_rplc_neg")
# 
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/expression_data_rplc_pos")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/expression_data_rplc_neg")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/expression_data_hilic_pos")
# load("3-data_analysis/monkey_fMSEA_analysis/4_group_features/expression_data_hilic_neg")

dir.create("3-data_analysis/monkey_fMSEA_analysis/6_pmd_analysis")
setwd("3-data_analysis/monkey_fMSEA_analysis/6_pmd_analysis")

library(pmd)
# data("spmeinvivo")
# str(spmeinvivo)

####RPLC positive
library(massdataset)
rplc_pos_object <-
  massdataset::create_mass_dataset(expression_data = expression_data_rplc_pos, 
                                   sample_info = sample_info_rplc_pos, 
                                   variable_info = variable_info_rplc_pos)
pmd_data  <-
  rplc_pos_object %>% 
  activate_mass_dataset(what = "sample_info") %>% 
  dplyr::filter(class == "QC" | class == "Subject") %>% 
  convert_mass_dataset2pmd()

pmd <- getpaired(pmd_data, rtcutoff = 5)
plotrtg(pmd)
plotpaired(pmd)

diff <- unique(pmd$paired$diff2)[1]
index <- pmd$paired$diff2 == diff
plotpaired(pmd,index)

std <- getstd(pmd)
plotstd(std)

plotpca(std$data[std$stdmassindex,],
        lv = as.numeric(as.factor(std$group$sample_group)),
        main = paste(sum(std$stdmassindex),"independent peaks"))

x <-
  data.frame(
    variable_id = rownames(std$data),
    mz = std$mz,
    rt = std$rt,
    rtcluster = std$rtcluster,
    index = std$stdmassindex
  ) %>%
  dplyr::arrange(rt) %>%
  plyr::dlply(.variables = .(rtcluster))

y <-
annotation_table %>% 
  dplyr::filter(mode == "RPLC" & polarity == "positive") %>% 
  dplyr::filter(variable_id_old %in% x[[1]]$variable_id) %>% 
  dplyr::left_join(x[[1]] %>% dplyr::select(-c(mz, rt)), 
                   by = c("variable_id_old" = "variable_id")) %>% 
  plyr::dlply(.variables = .(compound_class))


