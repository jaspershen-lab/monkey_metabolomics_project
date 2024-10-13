no_source()
####This should be run in workstation

library(tidymass)
library(tidyverse)

rm(list = ls())

setwd(r4projects::get_project_wd())
source("1_code/100_tools.R")

load("3_data_analysis/monkey_metabolomics_data_preparation/metabolite/object")

dir.create("3_data_analysis/monkey_sex")
setwd("3_data_analysis/monkey_sex")

# ####match samples
# object <-
# object %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(class == "Subject")
#
# sample_info <-
#   extract_sample_info(object)
#
# sample_info_male <-
#   sample_info %>%
#   dplyr::filter(sex == "M")
#
# sample_info_female <-
#   sample_info %>%
#   dplyr::filter(sex == "F")
#
#
# match_result <-
# seq_len(nrow(sample_info_male)) %>%
#   purrr::map(function(i){
#     age_error <-
#       abs(sample_info_male$age[i] - sample_info_female$age)
#     data.frame(sample_id1 = sample_info_male$sample_id[i],
#                sample_id2 = sample_info_female$sample_id,
#                age_error) %>%
#       dplyr::filter(!is.na(age_error)) %>%
#       dplyr::arrange(age_error) %>%
#       dplyr::filter(age_error < 1) %>%
#       head(10)
#   })
#
# temp <- NULL
# final_match_result <- vector(mode = "list", length = length(match_result))
#   for(i in 1:length(match_result)){
#   final_match_result[[i]] <-
#     match_result[[i]] %>%
#       dplyr::filter(!sample_id2 %in% temp$sample_id2) %>%
#       head(1)
#   temp <- rbind(
#     temp,
#     final_match_result[[i]]
#   )
#   }
#
#
# final_match_result <-
#   final_match_result %>%
#   dplyr::bind_rows()
#
# unique(final_match_result$sample_id1)
# unique(final_match_result$sample_id2)
#
#
#
# object <-
#   object %>%
#   activate_mass_dataset(what = "sample_info") %>%
#   dplyr::filter(sample_id %in% c(
#     final_match_result$sample_id1,
#     final_match_result$sample_id2
#   ))
# save(final_match_result, file = "final_match_result)
# save(object, file = "object")

load("object")

####PCA analysis
pca_object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject") %>%
  `+`(1) %>%
  log(2) %>%
  scale %>%
  massstat::run_pca()

plot <-
  massstat::pca_score_plot(
    object = object,
    pca_object = pca_object,
    color_by = "sex",
    frame = FALSE
  ) +
  scale_fill_manual(values = sex_color)

plot

plot
# ggsave(plot, filename = "pca_sex_plot.pdf",
#        width = 9, height = 7)

####heatmap to show the different heatmap
library(ComplexHeatmap)

temp_object <-
  object %>%
  `+`(1) %>%
  log(2) %>%
  scale()

library(circlize)

# cor_fun <-
#   circlize::colorRamp2(
#     colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
#     breaks = seq(
#       min(extract_sample_info(temp_object)$age),
#       max(extract_sample_info(temp_object)$age),
#       length.out = 9
#     )
#   )

h1 <-
  HeatmapAnnotation(sex = extract_sample_info(temp_object)$sex,
                    col = list(sex = sex_color))

variable_info <-
  extract_variable_info(temp_object)

load(
  here::here(
    "3_data_analysis/monkey_metabolomics_metabolite_annotation/hmdb_ms1.rda"
  )
)

variable_info <-
  variable_info %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Super_class", "Class")],
                   by = "HMDB.ID")

h2 <-
  rowAnnotation(
    Super_class = variable_info$Super_class,
    col = list(Super_class = metabolite_super_class_color)
  )

plot <-
  ComplexHeatmap::Heatmap(
    matrix = temp_object@expression_data,
    name = "z-score",
    show_column_names = FALSE,
    show_row_names = FALSE,
    row_names_gp = gpar(cex = 0.2),
    column_names_gp = gpar(cex = 0.2),
    top_annotation = h1,
    right_annotation = h2,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D",
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D"
  )

plot <- ggplotify::as.ggplot(plot)
plot
# ggsave(plot, filename = "all_heatmap.pdf", width = 10, height = 7)

######wilcox test
temp_object <-
  object %>%
  `+`(1) %>%
  log(2) %>%
  scale()

male_sample_id <-
  temp_object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(sex == "M") %>%
  pull(sample_id)

female_sample_id <-
  temp_object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(sex == "F") %>%
  pull(sample_id)

temp_object <-
  temp_object %>%
  massstat::mutate_p_value(
    control_sample_id = male_sample_id,
    case_sample_id = female_sample_id,
    method = "t.test",
    p_adjust_methods = "BH"
  ) %>%
  massstat::mutate_fc(
    control_sample_id = male_sample_id,
    case_sample_id = female_sample_id,
    mean_median = "median"
  )

###volcano plot
plot <- 
massstat::volcano_plot(
  object = temp_object,
  fc_column_name = "fc",
  p_value_column_name = "p_value", 
  fc_up_cutoff = 1, 
  fc_down_cutoff = 1, 
  point_size_scale = "p_value", 
  add_text = TRUE,
  text_from = "Compound.name",
) +
  scale_color_manual(values = c("UP" = unname(marker_color["Up"]),
                                "DOWN" = unname(marker_color["Down"]),
                                "NO" = unname(marker_color["No"])))

plot

variable_info <-
  extract_variable_info(temp_object)


# ggsave(plot,
#        filename = "volcano_plot.pdf",
#        width = 9,
#        height = 7)

####differential expressional metabolites
dim(variable_info)

sex_marker <-
  variable_info %>%
  dplyr::arrange(p_value_adjust) %>% 
  head(20)

pca_object <-
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject") %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% sex_marker$variable_id) %>%
  `+`(1) %>%
  log(2) %>%
  scale %>%
  massstat::run_pca()

plot <-
  massstat::pca_score_plot(
    object = object %>%
      activate_mass_dataset(what = "sample_info") %>%
      dplyr::filter(class == "Subject"),
    pca_object = pca_object,
    color_by = "sex",
    frame = TRUE
  ) +
  scale_fill_manual(values = sex_color) +
  scale_color_manual(values = sex_color) 

plot

plot1 <-
  pca_object$x %>%
  as.data.frame() %>%
  dplyr::select(PC1) %>%
  tibble::rownames_to_column(var = "sample_id") %>%
  dplyr::left_join(object@sample_info[, c("sample_id", "age", "sex")],
                   by = "sample_id") %>%
  ggplot(aes(x = PC1, y = age)) +
  geom_bar(
    aes(fill = sex),
    stat = "identity",
    width = 0.1,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = sex_color) +
  theme_base +
  labs(x = "", y = "Age (years)") +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

library(patchwork)

plot <-
  plot1 + plot + plot_layout(ncol = 1, heights = c(1, 5))
plot
# ggsave(plot, filename = "pca_sex_plot_with_markers.pdf",
#        width = 9, height = 7)

####heatmap to show the different heatmap
library(ComplexHeatmap)

temp_object2 <-
  temp_object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(class == "Subject") %>%
  dplyr::filter(!is.na(age)) %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(variable_id %in% sex_marker$variable_id)

library(circlize)

cor_fun <-
  circlize::colorRamp2(
    colors = RColorBrewer::brewer.pal(n = 9, name = "Reds"),
    breaks = seq(
      min(extract_sample_info(temp_object2)$age),
      max(extract_sample_info(temp_object2)$age),
      length.out = 9
    )
  )

h1 <-
  HeatmapAnnotation(sex = extract_sample_info(temp_object2)$sex,
                    col = list(sex = sex_color))

variable_info <-
  extract_variable_info(temp_object2)

variable_info <-
  variable_info %>%
  dplyr::left_join(hmdb_ms1@spectra.info[, c("HMDB.ID", "Super_class", "Class")],
                   by = "HMDB.ID")

lipid_class = variable_info$Class
lipid_class[variable_info$Super_class != "Lipids and lipid-like molecules"] <-
  NA

h2 <-
  rowAnnotation(
    Super_class = variable_info$Super_class,
    lipid_class = lipid_class,
    col = list(Super_class = metabolite_super_class_color,
               lipid_class = lipid_class_color)
  )

plot <-
  ComplexHeatmap::Heatmap(
    matrix = temp_object2@expression_data,
    name = "z-score",
    show_column_names = FALSE,
    show_row_names = FALSE,
    row_names_gp = gpar(cex = 0.2),
    column_names_gp = gpar(cex = 0.2),
    top_annotation = h1,
    right_annotation = h2,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D",
    border = TRUE
    # km = 2,
    # column_split = 2
  )

plot <- ggplotify::as.ggplot(plot)
plot

# ggsave(plot,
#        filename = "marker_heatmap.pdf",
#        width = 10,
#        height = 7)




