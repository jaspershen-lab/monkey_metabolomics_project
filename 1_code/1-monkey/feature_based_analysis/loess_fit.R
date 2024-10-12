no_function()

###
setwd(r4projects::get_project_wd())
rm(list = ls())
source("1-code/tools.R")
load("3-data_analysis/monkey_metabolomics_data_preparation/peak/object")

dir.create("3-data_analysis/monkey_feature_analysis/loess_fit",
           recursive = TRUE)
setwd("3-data_analysis/monkey_feature_analysis/loess_fit")

object

dim(object)

object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(!is.na(age))

####fore each molecule, we just use the loess method to
####impute the space between the first sample and the last sample

expression_data <-
  extract_expression_data(object)
sample_info <-
  extract_sample_info(object)
variable_info <-
  extract_variable_info(object)

dim(expression_data)
library(plyr)

# new_expression_data <-
#   vector(mode = "list", length = nrow(variable_info))
#
# for (i in seq_along(variable_info$variable_id)) {
#   temp_variable_id <- variable_info$variable_id[i]
#   cat(i, " ")
#   temp_data <-
#     data.frame(value = as.numeric(expression_data[temp_variable_id, ]),
#                sample_info)
#
#   optimize_span <-
#     optimize_loess_span(
#       x = temp_data$age,
#       y = temp_data$value,
#       span_range = c(0.3, 0.4, 0.5, 0.6)
#     )
#
#   span <-
#     optimize_span[[1]]$span[which.min(optimize_span[[1]]$rmse)]
#
#   value <- temp_data$value
#   age <- temp_data$age
#
#   ls_reg <-
#     loess(value ~ age,
#           span = span)
#
#   prediction_value =
#     predict(ls_reg,
#             newdata = data.frame(age = age))
#   new_expression_data[[i]] <- as.numeric(prediction_value)
# }
#
# new_expression_data <-
#   new_expression_data %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
#
#
# colnames(new_expression_data) <- colnames(expression_data)
# rownames(new_expression_data) <- rownames(expression_data)
#
# save(new_expression_data, file = "new_expression_data")

load("new_expression_data")

negative_idx <-
  new_expression_data %>%
  apply(1, function(x) {
    sum(x < 0)
  }) %>%
  `>`(0) %>%
  which()

plot(as.numeric(new_expression_data[negative_idx[1], ]))

new_expression_data[which(new_expression_data < 0, arr.ind = TRUE)] <-
  0

# cor_data <-
#   seq_len(nrow(new_expression_data)) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     temp <-
#       cor.test(as.numeric(expression_data[i,]),
#                as.numeric(new_expression_data[i,]),
#                method = "spearman")
#     data.frame(cor = as.numeric(temp$estimate),
#                p_value = as.numeric(temp$p.value))
#   }) %>%
#   dplyr::bind_rows() %>%
#   as.data.frame()
#
# cor_data$variable_id <-
#   variable_info$variable_id
#
# save(cor_data, file = "cor_data")
#
load("cor_data")

new_sample_info <-
  sample_info

new_variable_info <-
  variable_info

save(new_sample_info, file = "new_sample_info")
save(new_variable_info, file = "new_variable_info")

library(massdataset)

object <-
  create_mass_dataset(
    expression_data = new_expression_data,
    sample_info = new_sample_info,
    variable_info = new_variable_info
  )

save(object, file = "object")




sum(new_expression_data < 0)
